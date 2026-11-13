dom_pop <- function() {
  local_seed(21)
  d <- data.frame(
    id = 1:400,
    site = rep(c("a", "b", "c", "d"), times = c(160, 120, 80, 40)),
    region = rep(c("east", "west"), 200),
    cl = rep(paste0("c", 1:40), each = 10)
  )
  d$y <- round(stats::runif(400, 10, 500))
  d
}

test_that("domain totals add up to the overall total", {
  d <- dom_pop()
  s <- draw(d, design_stratified("site", n = 80), seed = 1, weights = TRUE)
  by_region <- ht_total(s, "y", by = "region")
  expect_s3_class(by_region, "drawn_by")
  expect_equal(sum(by_region$total), ht_total(s, "y")$total)
  expect_equal(sort(by_region$region), c("east", "west"))
  expect_equal(sum(by_region$n), nrow(s))
})

test_that("domain estimates are unbiased across repeated samples", {
  d <- dom_pop()
  truth <- tapply(d$y, d$region, sum)
  des <- design_stratified("site", n = 80)
  est <- vapply(1:300, function(i) {
    r <- ht_total(draw(d, des, seed = i, weights = TRUE), "y", by = "region",
                  variance = "none")
    r$total[order(r$region)]
  }, numeric(2))
  expect_equal(rowMeans(est), as.numeric(truth), tolerance = 0.02)
})

test_that("domain standard errors match survey::svyby, design by design", {
  skip_if_not_installed("survey")
  d <- dom_pop()
  for (des in list(design_simple(n = 80), design_stratified("site", n = 80),
                   design_cluster("cl", n_clusters = 8))) {
    s <- draw(d, des, seed = 2, weights = TRUE)
    svd <- as_svydesign(s)

    ours <- ht_total(s, "y", by = "region")
    theirs <- survey::svyby(~y, ~region, svd, survey::svytotal)
    expect_equal(ours$total, unname(stats::coef(theirs)), tolerance = 1e-8)
    expect_equal(ours$se, unname(survey::SE(theirs)), tolerance = 1e-6)

    ours_m <- ht_mean(s, "y", by = "region")
    theirs_m <- survey::svyby(~y, ~region, svd, survey::svymean)
    expect_equal(ours_m$mean, unname(stats::coef(theirs_m)), tolerance = 1e-8)
    expect_equal(ours_m$se, unname(survey::SE(theirs_m)), tolerance = 1e-6)
  }
})

test_that("the Hajek domain variance tracks the real sampling variance", {
  d <- dom_pop()
  des <- design_simple(n = 60)
  out <- vapply(1:500, function(i) {
    r <- ht_mean(draw(d, des, seed = i, weights = TRUE), "y", by = "region")
    r <- r[r$region == "west", ]
    c(r$mean, r$se^2)
  }, numeric(2))
  ratio <- mean(out[2, ]) / stats::var(out[1, ])
  expect_gt(ratio, 0.8)
  expect_lt(ratio, 1.25)
})

test_that("the HT domain mean divides by the domain's known size", {
  d <- dom_pop()
  s <- draw(d, design_simple(n = 80), seed = 3, weights = TRUE)
  r <- ht_mean(s, "y", by = "site", estimator = "ht")
  N_a <- sum(d$site == "a")
  expect_equal(r$mean[r$site == "a"],
               sum(s$y[s$site == "a"] * s$.weight[s$site == "a"]) / N_a)
})

test_that("by accepts several columns and validates them", {
  d <- dom_pop()
  s <- draw(d, design_simple(n = 80), seed = 4, weights = TRUE)
  two <- ht_total(s, "y", by = c("site", "region"))
  expect_true(all(c("site", "region") %in% names(two)))
  expect_equal(sum(two$total), ht_total(s, "y")$total)
  expect_error(ht_total(s, "y", by = "nope"),
               "`by` names `nope`, which `sample` does not have.", fixed = TRUE)
  expect_error(ht_total(s, "y", by = 1), "`by` must name")
  expect_output(print(two), "by domain")
})

test_that("a recoded stratum column does not change domain variances", {
  d <- dom_pop()
  s <- draw(d, design_stratified("site", n = 80, min_per_stratum = 2), seed = 4,
            weights = TRUE)
  recoded <- s
  recoded$site <- "z"
  for (v in c("analytic", "jackknife")) {
    want <- ht_total(s, "y", by = "region", variance = v)
    got <- ht_total(recoded, "y", by = "region", variance = v)
    expect_equal(got$se, want$se, label = v)
    expect_equal(attr(got, "df"), attr(want, "df"), label = v)
  }
})

test_that("domains come out in key order, not sorted as text", {
  d <- dom_pop()
  d$month <- rep(1:12, length.out = 400)
  d$mfac <- factor(month.abb[d$month], levels = month.abb)
  d$case <- c("a", "B", "b")[(d$id %/% 12) %% 3 + 1]
  s <- draw(d, design_simple(n = 120), seed = 1, weights = TRUE)
  expect_equal(ht_total(s, "y", by = "month")$month, 1:12)
  expect_equal(as.character(ht_total(s, "y", by = "mfac")$mfac), month.abb)
  expect_equal(ht_total(s, "y", by = "case")$case, c("B", "a", "b"))
  with_other_collation(
    expect_equal(ht_total(s, "y", by = "case")$case, c("B", "a", "b"))
  )
  # Several columns: the first varies slowest.
  two <- ht_total(s, "y", by = c("mfac", "case"))
  expect_gt(nrow(two), 12)
  expect_identical(order(as.integer(two$mfac), match(two$case, c("B", "a", "b"))),
                   seq_len(nrow(two)))
})

test_that("domains and summaries of an sf sample carry no geometry", {
  pts <- sf_points()
  plain <- sf::st_drop_geometry(pts)
  des <- design_cluster("cl", n_clusters = 5)
  s <- draw(pts, des, seed = 1, weights = TRUE)
  by_site <- ht_total(s, "y", by = "site")
  expect_false(inherits(by_site, "sf"))
  expect_false("geometry" %in% names(by_site))
  want <- ht_total(draw(plain, des, seed = 1, weights = TRUE), "y",
                   by = "site")
  expect_equal(as.data.frame(by_site), as.data.frame(want))
  expect_equal(sample_summary(s)$by_group,
               sample_summary(draw(plain, des, seed = 1,
                                   weights = TRUE))$by_group)
  st <- draw(pts, design_stratified("site", n = 12), seed = 1, weights = TRUE)
  expect_equal(sample_summary(st)$by_group$in_frame, c(30L, 20L, 10L))
})

# ---- domains with too little to measure --------------------------------------

test_that("a Hajek domain mean inside one sampling unit has no se", {
  # Its linearised residuals sum to zero inside the one unit, so every
  # estimator gave a standard error of 0 and a zero-width "analytic" interval.
  d <- dom_pop()
  d$rare <- ifelse(d$id %% 40 == 0, "rare", "common")
  des <- design_simple(n = 40)
  hits <- vapply(1:100, function(i) {
    sum(draw(d, des, seed = i)$rare == "rare")
  }, integer(1))
  s <- draw(d, des, seed = which(hits == 1L)[1], weights = TRUE)
  expect_equal(sum(s$rare == "rare"), 1L)
  m <- ht_mean(s, "y", by = "rare")
  expect_true(is.na(m$se[m$rare == "rare"]))
  expect_equal(m$method[m$rare == "rare"], "none")
  expect_true(is.finite(m$se[m$rare == "common"]))
  expect_match(attr(m, "notes"), "The domain `rare` has one sampled row",
               fixed = TRUE, all = FALSE)
  # A domain total from one row is still estimated.
  expect_true(all(is.finite(ht_total(s, "y", by = "rare")$se)))

  # Ten rows of one cluster are still one unit.
  d$part <- ifelse(d$cl == "c3" & d$id %% 2 == 0, "part", "rest")
  cs <- draw(d, design_cluster("cl", n_clusters = 8), seed = 1, weights = TRUE)
  if (!"c3" %in% cs$cl) {
    seeds <- vapply(1:200, function(i) {
      "c3" %in% draw(d, design_cluster("cl", n_clusters = 8), seed = i)$cl
    }, logical(1))
    cs <- draw(d, design_cluster("cl", n_clusters = 8),
               seed = which(seeds)[1], weights = TRUE)
  }
  cm <- ht_mean(cs, "y", by = "part")
  expect_equal(cm$n[cm$part == "part"], 5)
  expect_true(is.na(cm$se[cm$part == "part"]))
  expect_match(attr(cm, "notes"), "in one cluster only", all = FALSE)
  expect_true(is.finite(cm$se[cm$part == "rest"]))
  for (v in c("analytic", "jackknife")) {
    expect_true(is.na(ht_mean(cs, "y", by = "part", variance = v)$se[1]))
  }
})

lonely_pop <- function() {
  # Neyman allocation gives stratum `b`, nearly constant in `aux`, one row of
  # 40 and takes all three rows of `c`.
  p <- data.frame(st = rep(c("a", "b", "c"), c(40, 40, 3)),
                  aux = c(seq(-1, 1, length.out = 40),
                          seq(-0.05, 0.05, length.out = 40), c(-100, 0, 100)),
                  y = 50 + (seq_len(83) * 37) %% 23)
  p$half <- p$st == "c" | seq_len(83) %% 2 == 0
  p
}

test_that("a domain clear of a lonely stratum is estimated, one across not", {
  p <- lonely_pop()
  des <- design_stratified("st", n = 20, allocation = "neyman",
                           allocation_by = "aux")
  s <- draw(p, des, seed = 1, weights = TRUE)
  expect_equal(as.vector(table(s$st)), c(16L, 1L, 3L))
  expect_true(all(s$.prob[s$st == "c"] == 1))

  for (v in c("auto", "analytic", "jackknife")) {
    r <- ht_total(s, "y", by = "st", variance = v)
    # Stratum `c` was taken whole, so its total is exact.
    expect_equal(r$se[r$st == "c"], 0, label = v)
    expect_equal(r$ci_lower[r$st == "c"], sum(p$y[p$st == "c"]), label = v)
    expect_equal(r$ci_upper[r$st == "c"], sum(p$y[p$st == "c"]), label = v)
    # Stratum `a` is unaffected by `b`, and `b` itself cannot be measured.
    expect_true(is.finite(r$se[r$st == "a"]), label = v)
    expect_true(is.na(r$se[r$st == "b"]), label = v)
  }
  # The analytic and jackknife figures for `a` agree.
  expect_equal(ht_total(s, "y", by = "st", variance = "jackknife")$se[1],
               ht_total(s, "y", by = "st", variance = "analytic")$se[1])

  # `half` cuts across `b`: its one sampled row is outside the TRUE domain,
  # but other rows of `b` are in it, so neither domain can be measured.
  expect_false(s$half[s$st == "b"])
  h <- ht_total(s, "y", by = "half")
  expect_true(all(is.na(h$se)))
  expect_match(attr(h, "notes"), "would understate", all = FALSE)

  # A `by` column the frame does not hold cannot be checked, so it is refused
  # and the note says the frame cannot tell.
  s$late <- s$st == "c"
  u <- ht_total(s, "y", by = "late")
  expect_true(all(is.na(u$se)))
  expect_match(attr(u, "notes"), "cannot be told from the frame",
               all = FALSE)
})

test_that("a domain of intervals taken whole is exact at 0 df", {
  # One row a day on days 1 to 5, taken whole at per_interval = 1, and five
  # rows a day on days 6 to 10, of which one is drawn.
  day <- c(1:5, rep(6:10, each = 5))
  p <- data.frame(ts = as.POSIXct("2024-01-01", tz = "UTC") +
                    (day - 1) * 86400 + 3600 * c(rep(0, 5), rep(0:4, 5)),
                  y = 10 * day + (seq_along(day) * 7) %% 5)
  p$early <- day <= 5
  des <- design_temporal("ts", from = "2024-01-01", to = "2024-01-11",
                         interval = 1, per_interval = 1, unit = "days")
  s <- draw(p, des, seed = 1, weights = TRUE)
  r <- ht_total(s, "y", by = "early")
  expect_equal(attr(r, "df"), 0)
  expect_equal(r$se[r$early], 0)
  expect_equal(r$ci_lower[r$early], sum(p$y[p$early]))
  expect_equal(r$ci_upper[r$early], sum(p$y[p$early]))
  expect_true(is.na(r$se[!r$early]))
  m <- ht_mean(s, "y", by = "early")
  expect_equal(m$ci_lower[m$early], mean(p$y[p$early]))
})

test_that("domain means decide a proportion from the frame, as ht_mean() does", {
  pop <- data.frame(id = 1:100, site = rep(c("a", "b"), 50),
                    kids = c(rep(0, 60), rep(1, 36), rep(3, 4)))
  s <- draw(pop, design_simple(n = 10), seed = 1, weights = TRUE)
  expect_true(all(s$kids %in% 0:1))
  expect_equal(attr(ht_mean(s, "kids", by = "site"), "ci_scale"), "linear")
  r <- ht_mean(s, "kids", by = "site", proportion = TRUE)
  expect_equal(attr(r, "ci_scale"), "logit")
  expect_output(print(r), "95% CI, logit, ")
})

test_that("a domain proportion's interval is the logit interval, in numbers", {
  d <- dom_pop()
  d$big <- d$y > 300
  s <- draw(d, design_stratified("site", n = 80), seed = 3, weights = TRUE)
  r <- ht_mean(s, "big", by = "region")
  df <- attr(r, "df")
  expect_equal(df, 76)
  p <- r$mean
  expect_true(all(p > 0 & p < 1))
  # The delta-method interval on the logit scale, mapped back.
  half <- stats::qt(0.975, df) * r$se / (p * (1 - p))
  expect_equal(r$ci_lower, stats::plogis(stats::qlogis(p) - half))
  expect_equal(r$ci_upper, stats::plogis(stats::qlogis(p) + half))
  # Both proportions are below a half, where the logit interval sits above
  # the symmetric one.
  expect_true(all(r$ci_lower > p - stats::qt(0.975, df) * r$se))

  skip_if_not_installed("survey")
  sv <- as_svydesign(s)
  theirs <- survey::svyby(~big, ~region, sv, survey::svyciprop,
                          method = "xlogit", vartype = "ci",
                          df = survey::degf(sv))
  expect_equal(r$ci_lower, unname(theirs$ci_l), tolerance = 1e-8)
  expect_equal(r$ci_upper, unname(theirs$ci_u), tolerance = 1e-8)
})

test_that("a domain proportion of 0 gets a beta interval and a note", {
  pop <- data.frame(id = 1:300, site = rep(c("a", "b", "c"), each = 100))
  pop$late <- pop$site == "c" & pop$id %% 4 == 0
  s <- draw(pop, design_stratified("site", n = 60), seed = 1, weights = TRUE)
  r <- ht_mean(s, "late", by = "site")
  # Twenty sampled rows in each domain, against the design's 57 degrees of
  # freedom.
  n_eff <- 20 * (stats::qt(0.025, 19) / stats::qt(0.025, 57))^2
  expect_equal(r$ci_lower[1:2], c(0, 0))
  expect_equal(r$ci_upper[1:2], rep(stats::qbeta(0.975, 1, n_eff), 2))
  expect_true(r$ci_lower[3] > 0 && r$ci_upper[3] < 1)
  expect_match(attr(r, "notes"), "The domains `a`, `b` have a proportion of 0",
               all = FALSE)
  expect_output(print(r), "The domains `a`, `b` have a proportion of 0")
})

test_that("a y constant over a domain gives its mean a variance of exactly 0", {
  pop <- data.frame(id = 1:200, site = rep(c("a", "b"), 100),
                    w = 1 + (1:200 * 7) %% 13)
  pop$y <- ifelse(pop$site == "a", 0.3, (1:200 * 11) %% 17)
  s <- draw(pop, design_weighted("w", n = 40, method = "poisson"), seed = 2,
            weights = TRUE)
  r <- ht_mean(s, "y", by = "site")
  expect_identical(r$mean[1], 0.3)
  expect_identical(r$se[1], 0)
  expect_true(r$se[2] > 0)
  expect_match(attr(r, "notes"), "every sampled row of the domain `a`",
               all = FALSE)
})

test_that("domain tables carry one-sided bounds and the unreachable count", {
  skew <- data.frame(id = 1:1000, g = rep(c("big", "rare"), c(990, 10)),
                     region = rep(c("east", "west"), 500),
                     y = 1 + (1:1000 * 7) %% 31)
  expect_warning(
    s <- draw(skew, design_stratified("g", n = 40), seed = 1, weights = TRUE),
    class = "drawn_empty_stratum")
  r <- ht_total(s, "y", by = "region", alternative = "less")
  two <- ht_total(s, "y", by = "region", level = 0.90)
  expect_equal(r$ci_upper, two$ci_upper)
  expect_equal(r$ci_lower, c(-Inf, -Inf))
  expect_equal(attr(r, "alternative"), "less")
  expect_equal(attr(r, "unreachable"), 10)
  expect_output(print(r), "95% upper bound, t, 39 df")
  expect_output(print(r), "unreachable  10 frame rows at probability 0")
  g <- ht_mean(s, "y", by = "region", alternative = "greater")
  expect_equal(g$ci_upper, c(Inf, Inf))
})

# ---- all domains at once, and contrasts --------------------------------------

contrast_pop <- function() {
  id <- 1:400
  d <- data.frame(
    id = id,
    site = rep(c("a", "b", "c", "d"), times = c(160, 120, 80, 40)),
    region = rep(c("east", "west"), 200),
    cl = rep(paste0("c", 1:40), each = 10),
    w = 1 + (id * 7) %% 19,
    ts = as.POSIXct("2024-01-01", tz = "UTC") + ((id - 1) %% 10) * 86400 +
      id * 60
  )
  d$y <- 10 + (id * 37) %% 491
  d$hit <- (id * 13) %% 7 < 3
  d
}

# Each domain's variance from its own run of the estimator, with nothing
# shared between domains: the way `by` worked before every domain was
# estimated against one context.
domain_loop <- function(s, y, by, what = "total", estimator = "hajek",
                        variance = "auto") {
  parts <- drawn:::ht_prepare(s, y, 0.95)
  key <- s[[by]]
  vapply(sort(unique(key)), function(d) {
    ind <- as.numeric(key == d)
    z <- if (what == "total") {
      ind * parts$y
    } else if (estimator == "ht") {
      ind * parts$y / sum(parts$pop[[by]] == d)
    } else {
      n_hat <- sum(ind / parts$pi)
      ind * (parts$y - sum(ind * parts$y / parts$pi) / n_hat) / n_hat
    }
    drawn:::ht_variance_dispatch(parts$design, parts$units, parts$pop,
                                 parts$rows, z, parts$pi, variance,
                                 parts$pop[[by]] == d)$variance
  }, numeric(1), USE.NAMES = FALSE)
}

test_that("domain variances worked out together equal one run per domain", {
  d <- contrast_pop()
  designs <- list(
    simple = design_simple(n = 60),
    reservoir = design_reservoir(n = 60),
    stratified = design_stratified("site", n = 60),
    temporal = design_temporal("ts", from = "2024-01-01", to = "2024-01-11",
                               interval = 1, per_interval = 6, unit = "days"),
    certainty = design_certainty("w", 18, design_stratified("site", n = 40)),
    cluster = design_cluster("cl", n_clusters = 8),
    multistage = design_multistage("cl", n_clusters = 8, n = 40),
    systematic = design_systematic(interval = 7),
    pps = design_weighted("w", n = 40, method = "systematic"),
    spread = design_spread("w", n = 40)
  )
  for (nm in names(designs)) {
    s <- draw(d, designs[[nm]], seed = 3, weights = TRUE)
    for (v in c("auto", "jackknife")) {
      lab <- paste(nm, v)
      tot <- ht_total(s, "y", by = "region", variance = v)
      expect_equal(tot$se^2, domain_loop(s, "y", "region", variance = v),
                   tolerance = 1e-10, label = lab)
      hajek <- ht_mean(s, "y", by = "region", variance = v)
      expect_equal(hajek$se^2,
                   domain_loop(s, "y", "region", "mean", variance = v),
                   tolerance = 1e-10, label = lab)
      ht <- ht_mean(s, "y", by = "site", estimator = "ht", variance = v)
      expect_equal(ht$se^2,
                   domain_loop(s, "y", "site", "mean", "ht", variance = v),
                   tolerance = 1e-10, label = lab)
    }
  }
})

test_that("domains with units in a lonely stratum fall back one by one", {
  p <- lonely_pop()
  p$third <- c("u", "v", "w")[seq_len(83) %% 3 + 1]
  des <- design_stratified("st", n = 20, allocation = "neyman",
                           allocation_by = "aux")
  s <- draw(p, des, seed = 1, weights = TRUE)
  s$late <- s$st == "c"
  for (by in c("st", "half", "third", "late")) {
    for (v in c("auto", "analytic")) {
      r <- ht_total(s, "y", by = by, variance = v)
      want <- domain_loop(s, "y", by, variance = v)
      expect_equal(is.na(r$se), is.na(want), label = paste(by, v))
      expect_equal(r$se^2, want, tolerance = 1e-10, label = paste(by, v))
    }
  }
})

test_that("200 domains of 2,000 sampled rows take well under 2 seconds", {
  skip_on_cran()
  N <- 40000
  pop <- data.frame(id = seq_len(N), st = sprintf("s%d", seq_len(N) %% 5),
                    dom = sprintf("d%03d", (seq_len(N) * 7919) %% 200))
  pop$y <- 100 + (pop$id * 37) %% 101
  s <- draw(pop, design_stratified("st", n = 2000), seed = 1, weights = TRUE)
  tot <- NULL
  mean_by <- NULL
  expect_lt(system.time(tot <- ht_total(s, "y", by = "dom"))[["elapsed"]], 2)
  expect_lt(system.time(mean_by <- ht_mean(s, "y", by = "dom"))[["elapsed"]],
            2)
  expect_equal(nrow(tot), 200L)
  expect_equal(tot$se^2, domain_loop(s, "y", "dom"), tolerance = 1e-10)
  expect_equal(mean_by$se^2, domain_loop(s, "y", "dom", "mean"),
               tolerance = 1e-10)
})

test_that("a contrast of domain totals is the total of the weighted variable", {
  d <- contrast_pop()
  w <- c(a = 1, b = -1, d = 0.5)
  for (des in list(design_stratified("site", n = 60),
                   design_cluster("cl", n_clusters = 8),
                   design_systematic(interval = 7))) {
    s <- draw(d, des, seed = 2, weights = TRUE)
    r <- ht_total(s, "y", by = "site", contrast = w)
    weighted <- ht_total(s, s$y * unname(c(w, c = 0)[s$site]))
    expect_s3_class(r, "drawn_contrast")
    expect_equal(r$contrast, weighted$total)
    expect_equal(r$se, weighted$se)
    expect_equal(r$ci, weighted$ci)
    expect_equal(r$method, weighted$method)
    by_site <- ht_total(s, "y", by = "site")
    expect_equal(r$estimates, stats::setNames(by_site$total[c(1, 2, 4)],
                                              c("a", "b", "d")))
    expect_equal(r$weights, w)
    expect_equal(r$df, attr(by_site, "df"))
  }

  # A Horvitz-Thompson mean divides each domain's total by its frame size.
  s <- draw(d, design_stratified("site", n = 60), seed = 2, weights = TRUE)
  m <- ht_mean(s, "y", by = "site", estimator = "ht",
               contrast = c(a = 1, b = -1))
  z <- s$y * ((s$site == "a") / 160 - (s$site == "b") / 120)
  expect_equal(m$contrast, ht_total(s, z)$total)
  expect_equal(m$se, ht_total(s, z)$se)
})

test_that("contrast standard errors match survey::svycontrast", {
  skip_if_not_installed("survey")
  d <- contrast_pop()
  for (des in list(design_stratified("site", n = 80),
                   design_cluster("cl", n_clusters = 8))) {
    s <- draw(d, des, seed = 2, weights = TRUE)
    svd <- as_svydesign(s)

    w <- c(east = 1, west = -1)
    ours <- ht_total(s, "y", by = "region", contrast = w)
    theirs <- survey::svycontrast(
      survey::svyby(~y, ~region, svd, survey::svytotal, covmat = TRUE), w)
    expect_equal(ours$contrast, as.numeric(stats::coef(theirs)))
    expect_equal(ours$se, as.numeric(survey::SE(theirs)), tolerance = 1e-8)

    w <- c(a = 1, b = -0.5, c = -0.5)
    ours <- ht_mean(s, "y", by = "site", contrast = w)
    theirs <- survey::svycontrast(
      survey::svyby(~y, ~site, svd, survey::svymean, covmat = TRUE), w)
    expect_equal(ours$contrast, as.numeric(stats::coef(theirs)))
    expect_equal(ours$se, as.numeric(survey::SE(theirs)), tolerance = 1e-8)
  }
})

test_that("a contrast of proportions has a symmetric interval", {
  d <- contrast_pop()
  s <- draw(d, design_stratified("site", n = 80), seed = 2, weights = TRUE)
  r <- ht_mean(s, "hit", by = "region", contrast = c(east = 1, west = -1))
  expect_equal(r$ci_scale, "linear")
  q <- stats::qt(0.975, r$df)
  expect_equal(r$ci, r$contrast + c(-1, 1) * q * r$se)
  up <- ht_mean(s, "hit", by = "region", contrast = c(east = 1, west = -1),
                alternative = "less")
  expect_equal(up$ci, c(-Inf, r$contrast + stats::qt(0.95, r$df) * r$se))
})

test_that("a contrast weighing a mean over one sampling unit has no se", {
  d <- contrast_pop()
  d$rare <- ifelse(d$id %% 40 == 0, "rare", "common")
  des <- design_simple(n = 40)
  hits <- vapply(1:100, function(i) {
    sum(draw(d, des, seed = i)$rare == "rare")
  }, integer(1))
  s <- draw(d, des, seed = which(hits == 1L)[1], weights = TRUE)
  r <- ht_mean(s, "y", by = "rare", contrast = c(rare = 1, common = -1))
  expect_true(is.na(r$se))
  expect_equal(r$method, "none")
  expect_match(r$note, "The domain `rare` has one sampled row", fixed = TRUE)
  expect_match(r$note, "so it has no standard error either", fixed = TRUE)
  # A contrast of totals is still estimated.
  expect_true(is.finite(ht_total(s, "y", by = "rare",
                                 contrast = c(rare = 1, common = -1))$se))
})

test_that("contrasts name domains as by does, and odd weights are refused", {
  d <- contrast_pop()
  d$empty <- ifelse(d$id == 400, "never", d$site)
  s <- draw(d, design_stratified("site", n = 60), seed = 2, weights = TRUE)
  expect_false("never" %in% s$empty)
  two <- ht_total(s, "y", by = c("site", "region"),
                  contrast = c("a . east" = 1, "a . west" = -1))
  expect_equal(names(two$weights), c("a . east", "a . west"))

  expect_error(ht_total(s, "y", contrast = c(a = 1, b = -1)),
               "`contrast` weighs the domains of `by`, and `by` is NULL",
               fixed = TRUE)
  expect_error(ht_total(s, "y", by = "site", contrast = c(1, -1)),
               "`contrast` must be a named numeric vector")
  expect_error(ht_total(s, "y", by = "site", contrast = c(a = "1")),
               "`contrast` must be a named numeric vector")
  expect_error(ht_total(s, "y", by = "site",
                        contrast = stats::setNames(c(1, -1), c("a", NA))),
               "`contrast` has a weight with no name")
  expect_error(ht_total(s, "y", by = "site", contrast = c(a = 1, a = -1)),
               "`contrast` names `a` more than once", fixed = TRUE)
  expect_error(ht_total(s, "y", by = "site", contrast = c(a = 1, b = NA)),
               "`contrast` gives `b` a weight of NA", fixed = TRUE)
  expect_error(ht_total(s, "y", by = "site", contrast = c(a = 1, z = -1)),
               paste("`contrast` names `z`, which is not a domain of `by`.",
                     "Name the weights after the domains, which are `a`,",
                     "`b`, `c`, `d`."), fixed = TRUE)
  expect_error(ht_total(s, "y", by = "empty", contrast = c(a = 1, never = -1)),
               "`contrast` names `never`, which has no sampled row",
               fixed = TRUE)
  expect_error(ht_total(s, "y", by = "site", contrast = c(a = 0, b = 0)),
               "`contrast` gives every domain a weight of 0", fixed = TRUE)
})

test_that("a contrast prints its weights, domains and interval", {
  d <- contrast_pop()
  s <- draw(d, design_stratified("site", n = 60), seed = 2, weights = TRUE)
  r <- ht_mean(s, "y", by = "site", contrast = c(a = 0.5, b = 0.5, c = -1))
  expect_output(print(r), "Contrast of Hajek means  (stratified design, n = 60)",
                fixed = TRUE)
  expect_output(print(r), "contrast 0.5 `a` + 0.5 `b` - `c`", fixed = TRUE)
  expect_output(print(r), "domains  `a` ", fixed = TRUE)
  expect_output(print(r), "95% CI", fixed = TRUE)
  t2 <- ht_total(s, "y", by = "site", contrast = c(b = -1, a = 1))
  expect_output(print(t2), "Contrast of Horvitz-Thompson totals")
  expect_output(print(t2), "contrast `a` - `b`", fixed = TRUE)
  expect_output(print(ht_mean(s, "y", by = "site", estimator = "ht",
                              contrast = c(a = -1, b = 1))),
                "contrast -`a` + `b`", fixed = TRUE)
})

# ---- the by table as numbers ------------------------------------------------

table_pop <- function() {
  id <- 1:300
  data.frame(id = id, site = rep(c("a", "b", "c"), times = c(150, 100, 50)),
             half = rep(c("east", "west"), 150),
             y = 20 + (id * 37) %% 101,
             flag = id <= 150 & id %% 3 == 0)
}

test_that("a by table holds each domain's variance and an integer count", {
  s <- draw(table_pop(), design_stratified("site", n = 60), seed = 1,
            weights = TRUE)
  b <- ht_total(s, "y", by = "site")
  expect_type(b$n, "integer")
  expect_equal(b$variance, b$se^2)
  expect_equal(names(b), c("site", "n", "total", "variance", "se",
                           "ci_lower", "ci_upper", "method"))
  expect_equal(attr(b, "variable"), "y")
  out <- capture.output(print(b))
  expect_false(any(grepl("variance", out)))
  expect_error(deff(b), "no design effect")
})

test_that("a by column named like a column of the table is refused", {
  # The count overwrote a `by` column called `n`, which then held the counts.
  d <- data.frame(id = 1:100, n = rep(c("x", "y"), 50), y = 1:100)
  s <- draw(d, design_simple(n = 20), seed = 1, weights = TRUE)
  expect_error(ht_total(s, "y", by = "n"), "`by` names `n`, which the table")
  s$variance <- s$n
  expect_error(ht_mean(s, "y", by = c("variance", "n")),
               paste0("`variance`, `n`, which the table of domain estimates ",
                      "uses for its own columns"), fixed = TRUE)
  # A contrast builds no table, so the name does not matter there.
  expect_s3_class(ht_total(s, "y", by = "n", contrast = c(x = 1, y = -1)),
                  "drawn_contrast")
})

test_that("a by table gives coef(), vcov(), confint() and tidy() by domain", {
  s <- draw(table_pop(), design_stratified("site", n = 60), seed = 1,
            weights = TRUE)
  b <- ht_mean(s, "y", by = "site")
  expect_equal(coef(b), c(a = b$mean[1], b = b$mean[2], c = b$mean[3]))
  v <- vcov(b)
  expect_equal(diag(v), c(a = b$variance[1], b = b$variance[2],
                          c = b$variance[3]))
  expect_equal(v[upper.tri(v)], c(0, 0, 0))
  expect_equal(unname(confint(b)), cbind(b$ci_lower, b$ci_upper))
  b8 <- ht_mean(s, "y", by = "site", level = 0.8)
  expect_equal(unname(confint(b, level = 0.8)),
               cbind(b8$ci_lower, b8$ci_upper))

  # A proportion keeps its logit interval, and the domains at 0 their beta
  # interval, at another level.
  p <- ht_mean(s, "flag", by = "site")
  p8 <- ht_mean(s, "flag", by = "site", level = 0.8)
  expect_equal(unname(confint(p, level = 0.8)),
               cbind(p8$ci_lower, p8$ci_upper))
  q <- ht_mean(s, "flag", by = "site", ci_method = "beta",
               alternative = "less")
  q8 <- ht_mean(s, "flag", by = "site", ci_method = "beta",
                alternative = "less", level = 0.8)
  expect_equal(unname(confint(q, level = 0.8)),
               cbind(q8$ci_lower, q8$ci_upper))

  # Domains of several columns are named as `contrast` names them.
  two <- ht_total(s, "y", by = c("site", "half"))
  expect_equal(names(coef(two))[1:2], c("a . east", "a . west"))

  skip_if_not_installed("survey")
  expect_equal(unname(survey::SE(b)), b$se)
  skip_if_not_installed("generics")
  t <- generics::tidy(b, conf.int = TRUE)
  expect_equal(names(t), c("site", "term", "estimate", "std.error",
                           "conf.low", "conf.high"))
  expect_equal(t$term, rep("y", 3))
  expect_equal(t$conf.high, b$ci_upper)
})

test_that("selecting rows of a by table keeps its methods; columns drop them", {
  s <- draw(table_pop(), design_stratified("site", n = 60), seed = 1,
            weights = TRUE)
  b <- ht_total(s, "y", by = "site")
  r <- b[c(3, 1), ]
  expect_equal(coef(r), coef(b)[c(3, 1)])
  expect_equal(confint(r, level = 0.8), confint(b, level = 0.8)[c(3, 1), ])
  expect_output(print(r), "by domain")

  # Selecting columns keeps the class but not the attributes the heading
  # reads, and print() stopped with "argument is of length zero".
  cols <- b[, c("site", "total")]
  out <- capture.output(print(cols))
  expect_false(any(grepl("by domain", out)))
  expect_true(any(grepl("site", out)) && any(grepl("total", out)))
  expect_output(print(subset(b, total > 0, select = c(site, se))), "se")
  expect_error(coef(cols), "attributes of a `by` table")
  expect_error(confint(cols), "attributes of a `by` table")

  plain <- as.data.frame(b)
  expect_identical(class(plain), "data.frame")
  expect_null(attr(plain, "what"))
  expect_equal(plain$total, b$total)
})

test_that("a contrast gives the same methods, named by its weights", {
  s <- draw(table_pop(), design_stratified("site", n = 60), seed = 1,
            weights = TRUE)
  k <- ht_total(s, "y", by = "site", contrast = c(a = 1, b = -1))
  expect_equal(coef(k), c("a - b" = k$contrast))
  expect_equal(vcov(k)[1, 1], k$variance)
  k8 <- ht_total(s, "y", by = "site", contrast = c(a = 1, b = -1),
                 level = 0.8)
  expect_equal(unname(confint(k, level = 0.8)[1, ]), k8$ci)
  df <- as.data.frame(k)
  expect_equal(df$domains, "a - b")
  expect_equal(df$contrast, k$contrast)
  expect_true(is.na(deff(k)))
})

test_that("several y by domain hold the covariances within each domain", {
  s <- draw(table_pop(), design_stratified("site", n = 60), seed = 1,
            weights = TRUE)
  s$v <- s$y %% 7
  b <- ht_total(s, c("y", "v"), by = "half")
  expect_equal(names(b)[1:3], c("variable", "half", "n"))
  one <- ht_total(s, "v", by = "half")
  expect_equal(b$total[b$variable == "v"], one$total)
  expect_equal(b$se[b$variable == "v"], one$se)
  V <- vcov(b)
  expect_equal(rownames(V), c("y . east", "y . west", "v . east",
                              "v . west"))
  expect_equal(unname(diag(V)), b$variance)
  expect_equal(V["y . east", "v . west"], 0)
  # Selecting rows keeps the matching part of the matrix.
  expect_equal(vcov(b[c(4, 1), ]), V[c(4, 1), c(4, 1)])
  expect_output(print(b), "by variable and domain")
  s$variable <- s$half
  expect_error(ht_total(s, c("y", "v"), by = "variable"),
               "`by` names `variable`")

  skip_if_not_installed("survey")
  svd <- as_svydesign(s)
  sb <- survey::svyby(~y + v, ~half, svd, survey::svytotal, covmat = TRUE)
  expect_equal(V["y . east", "v . east"], vcov(sb)["east:y", "east:v"])
  expect_equal(V["y . west", "v . west"], vcov(sb)["west:y", "west:v"])
  m <- ht_mean(s, c("y", "v"), by = "half")
  sm <- survey::svyby(~y + v, ~half, svd, survey::svymean, covmat = TRUE)
  expect_equal(vcov(m)["y . west", "v . west"],
               vcov(sm)["west:y", "west:v"])
})

test_that("domain ratios match survey::svyby(svyratio)", {
  d <- dom_pop()
  d$x <- 1 + d$id %% 9
  s <- draw(d, design_stratified("site", n = 80), seed = 2, weights = TRUE)
  b <- ht_ratio(s, "y", "x", by = "region")
  expect_equal(names(b), c("region", "n", "ratio", "variance", "se",
                           "ci_lower", "ci_upper", "method"))
  for (r in c("east", "west")) {
    i <- s$region == r
    expect_equal(b$ratio[b$region == r],
                 sum(s$y[i] * s$.weight[i]) / sum(s$x[i] * s$.weight[i]))
  }
  expect_output(print(b), "Ratio `y` / `x` by domain", fixed = TRUE)
  s$x0 <- ifelse(s$region == "east", 0, s$x)
  expect_error(ht_ratio(s, "y", "x0", by = "region"),
               "The estimated total of `x0` is 0 in the domain `east`")
  skip_if_not_installed("generics")
  expect_equal(generics::tidy(b)$term, c("y/x", "y/x"))

  skip_if_not_installed("survey")
  for (des in list(design_stratified("site", n = 80),
                   design_cluster("cl", n_clusters = 8))) {
    s <- draw(d, des, seed = 2, weights = TRUE)
    ours <- ht_ratio(s, "y", "x", by = "region")
    theirs <- survey::svyby(~y, ~region, as_svydesign(s), survey::svyratio,
                            denominator = ~x)
    expect_equal(ours$ratio, unname(coef(theirs)))
    expect_equal(ours$se, unname(survey::SE(theirs)))
  }
})
