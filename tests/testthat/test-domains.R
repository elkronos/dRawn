dom_pop <- function() {
  set.seed(21)
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
