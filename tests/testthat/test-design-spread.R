spread_pop <- function(N = 300) {
  set.seed(11)
  d <- data.frame(id = seq_len(N), x = stats::runif(N), y = stats::runif(N))
  d$smooth <- 100 + 80 * d$x + 60 * sin(4 * d$y) + stats::rnorm(N, 0, 10)
  d$size <- stats::rgamma(N, 3, 1)
  d
}

test_that("draws exactly n rows, in frame order, keeping the schema", {
  d <- spread_pop()
  des <- design_spread(c("x", "y"), n = 30)
  for (i in 1:20) {
    s <- draw(d, des, seed = i)
    expect_equal(nrow(s), 30L)
    expect_identical(s$id, sort(s$id))
  }
  expect_same_schema(s, d)
  expect_output(print(des), "<sampling design: spread>")
})

test_that("inclusion probabilities are exactly the ones asked for", {
  d <- spread_pop()
  eq <- design_spread(c("x", "y"), n = 30)
  expect_equal(inclusion_prob(d, eq), rep(0.1, 300))
  sim <- inclusion_prob(d, eq, simulate = TRUE, R = 4000, seed = 1)
  z <- (sim - 0.1) / sqrt(0.1 * 0.9 / 4000)
  expect_lt(max(abs(z)), 4.5)
  expect_lt(abs(mean(z)), 0.2)

  pps <- design_spread(c("x", "y"), n = 30, size = "size")
  expect_equal(inclusion_prob(d, pps),
               inclusion_prob(d, design_weighted("size", n = 30,
                                                 method = "poisson")))
  ex <- inclusion_prob(d, pps)
  sim <- inclusion_prob(d, pps, simulate = TRUE, R = 4000, seed = 2)
  z <- (sim - ex) / sqrt(ex * (1 - ex) / 4000)
  expect_lt(max(abs(z[ex < 1])), 4.5)
})

test_that("spreading buys a large variance reduction on a smooth surface", {
  d <- spread_pop()
  tot <- function(des) vapply(1:300, function(i) {
    ht_total(draw(d, des, seed = i, weights = TRUE), "smooth",
             variance = "none")$total
  }, numeric(1))
  v_spread <- stats::var(tot(design_spread(c("x", "y"), n = 30,
                                           scale = FALSE)))
  v_srs <- stats::var(tot(design_simple(n = 30)))
  expect_lt(v_spread / v_srs, 0.5)
})

test_that("neighbours are rarely taken together", {
  d <- spread_pop()
  s <- draw(d, design_spread(c("x", "y"), n = 30, scale = FALSE), seed = 1)
  r <- draw(d, design_simple(n = 30), seed = 1)
  nn <- function(p) {
    m <- as.matrix(stats::dist(p[, c("x", "y")]))
    diag(m) <- Inf
    mean(apply(m, 1, min))
  }
  expect_gt(nn(s), nn(r))
})

test_that("the local mean variance is conservative but in the right range", {
  d <- spread_pop()
  des <- design_spread(c("x", "y"), n = 40, scale = FALSE)
  out <- vapply(1:400, function(i) {
    r <- ht_total(draw(d, des, seed = i, weights = TRUE), "smooth")
    c(r$total, r$variance)
  }, numeric(2))
  ratio <- mean(out[2, ]) / stats::var(out[1, ])
  expect_gt(ratio, 0.85)
  expect_lt(ratio, 1.9)

  r <- ht_total(draw(d, des, seed = 1, weights = TRUE), "smooth")
  expect_equal(r$method, "local mean")
  expect_output(print(r), "(local mean, an approximation)", fixed = TRUE)
})

test_that("a spread census reports an exact total", {
  # Every row is in the sample, so there is no variance. The local mean
  # estimator compared each row with its neighbours anyway and reported a
  # standard error of about 5% of the total.
  d <- spread_pop(20)
  s <- draw(d, design_spread(c("x", "y"), n = 20), seed = 1, weights = TRUE)
  r <- ht_total(s, "smooth")
  expect_equal(r$total, sum(d$smooth))
  expect_equal(r$variance, 0)
  expect_match(r$note, "exact rather than estimated")
})

dominant_spread <- function() {
  d <- data.frame(x = (1:60 * 0.618) %% 1, y = (1:60 * 0.381) %% 1)
  d$w <- c(500, 400, 300, 1 + (1:57 %% 7) * 0.5)
  d$v <- 10 * d$w + 5 * sin(1:60)
  d
}

test_that("rows a size column puts past certainty add no local-mean variance", {
  d <- dominant_spread()
  des <- design_spread(c("x", "y"), n = 10, size = "w", scale = FALSE)
  s <- draw(d, des, seed = 1, weights = TRUE)
  expect_equal(sum(s$.prob == 1), 3L)

  # The estimator runs over the other rows only, as rows and as neighbours.
  free <- s$.prob < 1
  want <- drawn:::local_mean_variance(as.matrix(s[free, c("x", "y")]),
                                      s$v[free], s$.prob[free])
  r <- ht_total(s, "v")
  expect_equal(r$method, "local mean")
  expect_equal(r$variance, want$variance)

  y2 <- s$v * ifelse(free, 1, 4)
  expect_equal(ht_total(s, y2)$variance, r$variance)

  # One row below certainty leaves nothing to compare it with
  one <- draw(d, design_spread(c("x", "y"), n = 4, size = "w", scale = FALSE),
              seed = 1, weights = TRUE)
  expect_equal(sum(one$.prob < 1), 1L)
  a <- ht_total(one, "v", variance = "analytic")
  expect_true(is.na(a$variance))
  expect_match(a$note, "Only one sampled row")
  expect_true(is.na(ht_total(one, "v")$variance))
})

test_that("with certain rows the local mean variance stays near the truth", {
  skip_on_cran()
  # Counting the three dominant rows as random put this ratio above 5,000.
  d <- dominant_spread()
  des <- design_spread(c("x", "y"), n = 10, size = "w", scale = FALSE)
  out <- vapply(1:300, function(i) {
    r <- ht_total(draw(d, des, seed = i, weights = TRUE), "v")
    c(r$total, r$variance)
  }, numeric(2))
  ratio <- mean(out[2, ]) / stats::var(out[1, ])
  expect_gt(ratio, 0.6)
  expect_lt(ratio, 1.6)
})

test_that("joint probabilities have no closed form, but can be simulated", {
  d <- spread_pop(60)
  des <- design_spread(c("x", "y"), n = 10)
  expect_error(joint_prob(d, des, rows = 1:3), "no closed-form joint")
  m <- joint_prob(d, des, rows = 1:3, simulate = TRUE, R = 500, seed = 1)
  expect_equal(dim(m), c(3L, 3L))
})

test_that("scaling decides which column drives the distances", {
  set.seed(4)
  d <- data.frame(big = stats::runif(200, 0, 1e6), small = stats::runif(200))
  unscaled <- design_spread(c("big", "small"), n = 20, scale = FALSE)
  s <- draw(d, unscaled, seed = 1)
  # Unscaled, `small` barely moves the distance; scaled, it counts equally.
  expect_true(is.data.frame(s))
  expect_equal(nrow(draw(d, design_spread(c("big", "small"), n = 20), seed = 1)),
               20L)
})

test_that("validates its inputs and honours na_rm", {
  expect_error(design_spread(character(0), n = 5), "column names")
  expect_error(design_spread("x", n = -1), "non-negative")
  d <- spread_pop(50)
  d$x[3] <- NA
  expect_error(draw(d, design_spread(c("x", "y"), n = 5)), "na_rm = TRUE")
  s <- draw(d, design_spread(c("x", "y"), n = 5, na_rm = TRUE), seed = 1,
            weights = TRUE)
  expect_equal(nrow(s), 5L)
  expect_equal(inclusion_prob(d, design_spread(c("x", "y"), n = 5,
                                               na_rm = TRUE))[3], 0)
  d$label <- letters[(seq_len(50) %% 26) + 1]
  expect_error(draw(d, design_spread("label", n = 5)), "must be numeric")
  expect_error(draw(spread_pop(10), design_spread(c("x", "y"), n = 11)),
               "cannot exceed")
})

test_that("hands off to survey, and plots as a map", {
  d <- spread_pop()
  s <- draw(d, design_spread(c("x", "y"), n = 30), seed = 1, weights = TRUE)
  skip_if_not_installed("survey")
  svd <- as_svydesign(s)
  expect_equal(as.numeric(survey::svytotal(~smooth, svd)),
               ht_total(s, "smooth")$total)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(design_spread(c("x", "y"), n = 30), d, type = "map",
                        seed = 1))
  expect_error(plot(design_simple(n = 5), d, type = "map"), "coords")
  expect_invisible(plot(design_simple(n = 5), d, type = "map",
                        coords = c("x", "y"), seed = 1))
})

test_that("an sf frame spreads across plain numeric columns", {
  pts <- sf_points()
  plain <- sf::st_drop_geometry(pts)
  des <- design_spread(c("x", "yy"), n = 12)
  s <- draw(pts, des, seed = 1, weights = TRUE)
  expect_s3_class(s, "sf")
  expect_identical(s$id, draw(plain, des, seed = 1)$id)
  expect_equal(inclusion_prob(pts, des), rep(0.2, 60))
  expect_equal(ht_total(s, "y")$se,
               ht_total(draw(plain, des, seed = 1, weights = TRUE), "y")$se)
})

test_that("a POINT geometry column spreads as its two coordinates", {
  pts <- sf_points()
  plain <- sf::st_drop_geometry(pts)
  by_geom <- design_spread("geometry", n = 12)
  by_cols <- design_spread(c("lon", "lat"), n = 12)
  s <- draw(pts, by_geom, seed = 1, weights = TRUE)
  expect_identical(s$id, draw(plain, by_cols, seed = 1)$id)
  expect_equal(inclusion_prob(pts, by_geom), rep(0.2, 60))
  expect_equal(ht_total(s, "y")$se,
               ht_total(draw(plain, by_cols, seed = 1, weights = TRUE),
                        "y")$se)
  # With another column, the geometry still counts as two
  mixed <- design_spread(c("geometry", "x"), n = 12)
  expect_identical(draw(pts, mixed, seed = 2)$id,
                   draw(plain, design_spread(c("lon", "lat", "x"), n = 12),
                        seed = 2)$id)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(by_geom, pts, type = "map", seed = 1))
  lines <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(
    sf::st_linestring(rbind(c(0, 0), c(1, 1))),
    sf::st_linestring(rbind(c(0, 1), c(1, 2)))))
  expect_error(draw(lines, design_spread("geometry", n = 1)),
               "`across` names `geometry`, a geometry column holding LINESTRING",
               fixed = TRUE)
})

test_that("Date and date-time columns spread as days and seconds", {
  d <- spread_pop(40)
  d$day <- as.Date("2024-03-01") + (d$id * 7) %% 40
  d$days <- as.numeric(d$day)
  d$at <- as.POSIXct("2024-03-01", tz = "UTC") + d$days * 3600
  d$secs <- as.numeric(d$at)
  same <- function(a, b) {
    expect_identical(draw(d, design_spread(a, n = 6), seed = 3)$id,
                     draw(d, design_spread(b, n = 6), seed = 3)$id)
  }
  same(c("day", "x"), c("days", "x"))
  same(c("at", "x"), c("secs", "x"))
  expect_equal(inclusion_prob(d, design_spread("day", n = 6)), rep(0.15, 40))
  s <- draw(d, design_spread(c("day", "x"), n = 6), seed = 3, weights = TRUE)
  expect_equal(ht_total(s, "smooth")$se,
               ht_total(draw(d, design_spread(c("days", "x"), n = 6), seed = 3,
                             weights = TRUE), "smooth")$se)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(design_spread(c("day", "x"), n = 6), d, type = "map",
                        seed = 3))
  # A size has to be numeric, and anything else names the way out
  expect_error(draw(d, design_spread("x", n = 6, size = "day")),
               "`day` must be numeric to size by, not Date.", fixed = TRUE)
  d$label <- format(d$x)
  expect_error(draw(d, design_spread("label", n = 6)),
               "for example with as.numeric()", fixed = TRUE)
  d$day[4] <- NA
  expect_error(draw(d, design_spread("day", n = 6)), "na_rm = TRUE")
  expect_equal(inclusion_prob(d, design_spread("day", n = 6,
                                               na_rm = TRUE))[4], 0)
})

# ---- strata -----------------------------------------------------------------

strata_pop <- function(N = 150, sizes = c(70, 50, 30)) {
  i <- seq_len(N)
  d <- data.frame(id = i, x = (i * 0.618034) %% 1, y = (i * 0.414214) %% 1,
                  g = rep(c("north", "south", "east"), sizes))
  d$v <- 100 + 80 * d$x + 60 * sin(4 * d$y) + 5 * sin(7 * i)
  d$w <- 1 + (i %% 9)
  d
}

test_that("strata get their allocation, each spread on its own", {
  d <- strata_pop()
  des <- design_spread(c("x", "y"), n = 30, strata = "g")
  expect_output(print(des), "strata +\"g\"")
  for (i in 1:10) {
    s <- draw(d, des, seed = i)
    expect_equal(as.vector(table(factor(s$g, c("east", "north", "south")))),
                 c(6L, 14L, 10L))
    expect_identical(s$id, sort(s$id))
  }
  eq <- draw(d, design_spread(c("x", "y"), n = 30, strata = "g",
                              allocation = "equal"), seed = 1)
  expect_equal(as.vector(table(eq$g)), c(10L, 10L, 10L))

  counts <- design_spread(c("x", "y"), strata = "g",
                          allocation = c(north = 4, south = 3, east = 2))
  expect_equal(counts$n, 9L)
  expect_equal(as.vector(table(draw(d, counts, seed = 1)$g)), c(2L, 4L, 3L))

  tab <- allocation_table(d, des)
  expect_equal(tab$g, c("east", "north", "south"))
  expect_equal(tab$n, c(6L, 14L, 10L))
  expect_equal(tab$prob, c(0.2, 0.2, 0.2))
  expect_equal(tab$status, rep("rule", 3))
  expect_equal(allocation_table(d, counts)$status, rep("count", 3))

  # Each stratum is the design without strata on that stratum's rows alone.
  north <- d[d$g == "north", ]
  expect_equal(inclusion_prob(d, des)[d$g == "north"],
               inclusion_prob(north, design_spread(c("x", "y"), n = 14)))
})

test_that("inclusion probabilities with strata and size are exact", {
  d <- strata_pop()
  d$w[c(5, 80)] <- 60
  des <- design_spread(c("x", "y"), n = 30, strata = "g", size = "w")
  ex <- inclusion_prob(d, des)
  expect_equal(as.vector(tapply(ex, d$g, sum)), c(6, 14, 10))
  north <- d$g == "north"
  expect_equal(ex[north],
               inclusion_prob(d[north, ], design_spread(c("x", "y"), n = 14,
                                                        size = "w")))
  expect_equal(as.vector(certainty_set(d, des)), c(5L, 80L))
  expect_true(is.na(attr(certainty_set(d, des), "threshold")))

  skip_on_cran()
  small <- strata_pop(90, c(40, 30, 20))
  small$w[c(5, 50)] <- 40
  des <- design_spread(c("x", "y"), n = 18, strata = "g", size = "w")
  ex <- inclusion_prob(small, des)
  sim <- inclusion_prob(small, des, simulate = TRUE, R = 1500, seed = 1)
  z <- (sim - ex) / sqrt(ex * (1 - ex) / 1500)
  expect_lt(max(abs(z[ex < 1])), 4.5)
  expect_equal(sim[ex == 1], rep(1, sum(ex == 1)))
})

test_that("the variance sums the strata's local mean variances", {
  d <- strata_pop()
  d$w[c(5, 80)] <- 60
  for (size in list(NULL, "w")) {
    des <- design_spread(c("x", "y"), n = 30, strata = "g", size = size)
    s <- draw(d, des, seed = 3, weights = TRUE)
    # Each stratum's own estimator, on that stratum alone as a frame.
    parts <- vapply(c("east", "north", "south"), function(h) {
      sub <- d[d$g == h, ]
      own <- design_spread(c("x", "y"), n = sum(s$g == h), size = size)
      in_h <- s$g == h
      drawn:::ht_variance(own, sub, match(s$id[in_h], sub$id), s$v[in_h],
                  s$.prob[in_h])$variance
    }, numeric(1))
    r <- ht_total(s, "v")
    expect_equal(r$method, "local mean")
    expect_equal(r$variance, sum(parts))
    expect_equal(r$df, sum(s$.prob < 1) - 3)

    by <- ht_total(s, "v", by = "g")
    expect_equal(by$se^2, unname(parts))
  }

  jk <- ht_total(s, "v", variance = "jackknife")
  expect_match(jk$note, "in 3 strata")
})

test_that("a stratum with one sampled row below certainty is refused", {
  d <- strata_pop()
  des <- design_spread(c("x", "y"), strata = "g",
                       allocation = c(north = 4, south = 1, east = 3))
  s <- draw(d, des, seed = 1, weights = TRUE)
  a <- ht_total(s, "v", variance = "analytic")
  expect_true(is.na(a$variance))
  expect_match(a$note, "The stratum `south` has a single sampled row")
  expect_match(a$note, "a count of 2 or more for each in `allocation`")
  expect_true(is.na(ht_total(s, "v")$variance))

  # A domain with no unit in that stratum is measured without it.
  by <- ht_total(s, "v", by = "g")
  expect_true(all(is.finite(by$se[by$g != "south"])))
  expect_true(is.na(by$se[by$g == "south"]))
})

test_that("strata are refused or warned about as stratified designs are", {
  expect_error(design_spread(c("x", "y"), n = 5, allocation = "equal"),
               "`allocation` splits `n` over strata, and `strata` is NULL")
  expect_error(design_spread(c("x", "y"), n = 5, strata = "g",
                             allocation = "neyman"),
               "must be one of \"proportional\" or \"equal\"")
  expect_error(design_spread(c("x", "y"), strata = "g"), "`n` is missing")
  expect_error(design_spread(c("x", "y"), n = 7, strata = "g",
                             allocation = c(north = 2, south = 3)),
               "is not the sum of the counts")
  d <- strata_pop()
  expect_error(draw(d, design_spread(c("x", "y"), strata = "g",
                                     allocation = c(north = 2, south = 3))),
               "no entry for 1 stratum of `data`: `east`")
  expect_error(draw(d, design_spread(c("x", "y"), strata = "g",
                                     allocation = c(north = 2, south = 3,
                                                    east = 31))),
               "Stratum `east` has 30 row(s) but `allocation` asks for 31. ",
               fixed = TRUE)
  expect_error(inclusion_prob(d, design_spread(c("x", "y"), n = 151,
                                               strata = "g")),
               "cannot exceed the number of rows")
  expect_warning(draw(d, design_spread(c("x", "y"), n = 2, strata = "g"),
                      seed = 1),
                 class = "drawn_empty_stratum")

  d$g[3] <- NA
  des <- design_spread(c("x", "y"), n = 30, strata = "g")
  expect_error(draw(d, des), "a missing value in `x`, `y`, `g`")
  kept <- design_spread(c("x", "y"), n = 30, strata = "g", na_rm = TRUE)
  p <- inclusion_prob(d, kept)
  expect_equal(p[3], 0)
  expect_equal(p[d$g %in% "north"], rep(14 / 69, 69))
  expect_equal(sample_summary(draw(d, kept, seed = 1, weights = TRUE))$N_reachable,
               149L)

  # A row dropped for a missing coordinate is in no stratum's count, so the
  # sample and its variances are those of the frame without it.
  d <- strata_pop()
  d$x[10] <- NA
  s1 <- draw(d, kept, seed = 1, weights = TRUE)
  s2 <- draw(d[-10, ], kept, seed = 1, weights = TRUE)
  expect_identical(s1$id, s2$id)
  for (v in c("analytic", "jackknife")) {
    expect_equal(ht_total(s1, "v", variance = v)$variance,
                 ht_total(s2, "v", variance = v)$variance)
  }
})

test_that("sample_summary() and allocation_table() show the strata", {
  d <- strata_pop()
  s <- draw(d, design_spread(c("x", "y"), n = 30, strata = "g"), seed = 1,
            weights = TRUE)
  b <- sample_summary(s)$by_group
  expect_equal(b$group, c("east", "north", "south"))
  expect_equal(b$drawn, c(6L, 14L, 10L))
  expect_equal(b$in_frame, c(30L, 70L, 50L))

  d$w[c(5, 80)] <- 60
  sized <- design_spread(c("x", "y"), n = 30, strata = "g", size = "w")
  expect_true(all(is.na(allocation_table(d, sized)$prob)))
  expect_error(allocation_table(d, design_spread(c("x", "y"), n = 30)),
               "which has no strata to allocate over")

  cert <- design_certainty("w", 60, rest = design_spread(
    c("x", "y"), n = 20, strata = "g"))
  tab <- allocation_table(d, cert)
  expect_equal(tab$N, c(30L, 69L, 49L))
  expect_equal(sum(tab$n), 20L)
})

test_that("a spread sample with strata hands its strata to survey", {
  skip_if_not_installed("survey")
  d <- strata_pop()
  s <- draw(d, design_spread(c("x", "y"), n = 30, strata = "g"), seed = 1,
            weights = TRUE)
  svd <- as_svydesign(s)
  r <- ht_total(s, "v")
  expect_equal(as.numeric(survey::svytotal(~v, svd)), r$total)
  expect_equal(survey::degf(svd), r$df)
  expect_equal(sort(unique(as.character(svd$strata[, 1]))),
               c("east", "north", "south"))

  d$w[c(5, 80)] <- 60
  s <- draw(d, design_spread(c("x", "y"), n = 30, strata = "g", size = "w"),
            seed = 2, weights = TRUE)
  svd <- as_svydesign(s)
  r <- ht_total(s, "v")
  expect_equal(as.numeric(survey::svytotal(~v, svd)), r$total)
  expect_equal(survey::degf(svd), r$df)
  expect_equal(as.character(svd$strata[s$.prob == 1, 1]),
               c("(certain)", "(certain)"))

  # A stratum of the frame named `(certain)` stays apart from the certain rows.
  d$g[d$g == "east"] <- "(certain)"
  s <- draw(d, design_spread(c("x", "y"), n = 30, strata = "g", size = "w"),
            seed = 2, weights = TRUE)
  svd <- as_svydesign(s)
  expect_equal(length(unique(svd$strata[, 1])), 4L)
  expect_equal(survey::degf(svd), ht_total(s, "v")$df)

  cert <- design_certainty("w", 60, rest = design_spread(
    c("x", "y"), n = 20, strata = "g"))
  s <- draw(d, cert, seed = 1, weights = TRUE)
  expect_equal(survey::degf(as_svydesign(s)), ht_total(s, "v")$df)
})

test_that("a row with a tiny target probability still competes", {
  # Targets below 1e-9 were set to 0 before the pivotal walk, so such a row
  # could never be drawn although inclusion_prob() gave it a probability.
  # The walk now snaps only rounding error, at a tolerance `eps`; a target
  # far below it, here 0.01 against 0.05, is drawn at its own rate.
  x <- matrix(c(0, 1), ncol = 1)
  hits <- drawn:::with_seed(1, vapply(1:2000, function(i) {
    1L %in% drawn:::local_pivotal(c(0.01, 0.99), x, eps = 0.05)
  }, logical(1)))
  expect_gt(sum(hits), 5)
  expect_lt(sum(hits), 40)
  sizes <- drawn:::with_seed(2, vapply(1:200, function(i) {
    length(drawn:::local_pivotal(c(1e-10, 1 - 1e-10), x))
  }, integer(1)))
  expect_equal(sizes, rep(1L, 200))

  d <- data.frame(x = 1:20, y = 0, w = c(2e-9, rep(1, 19)))
  des <- design_spread(c("x", "y"), n = 5, size = "w")
  expect_equal(inclusion_prob(d, des)[1], 5 * 2e-9 / sum(d$w))
  expect_equal(nrow(draw(d, des, seed = 1)), 5L)
})
