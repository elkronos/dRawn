# The variance estimators the design calls for, and the cases where the naive
# one silently understated or overstated. Every claim here was a wrong number
# before it was a test.

strata_pop <- function() {
  set.seed(71)
  d <- data.frame(id = 1:400, g = rep(1:20, each = 20))
  d$y <- stats::rnorm(400, d$g * 3, 5)
  d
}

# ---- lonely strata -----------------------------------------------------------

test_that("one sampled row per stratum is refused, not reported as zero", {
  d <- strata_pop()
  s <- draw(d, design_stratified("g", n = 20), seed = 1, weights = TRUE)
  expect_true(all(table(s$g) == 1L))

  # Sen-Yates-Grundy has no within-stratum pair to sum over, and used to
  # report a variance of exactly 0 -- a zero-width interval -- here.
  r <- ht_total(s, "y", variance = "analytic")
  expect_true(is.na(r$variance))
  expect_match(r$note, "single sampled row")
  expect_match(r$note, "min_per_stratum = 2")

  # The stratified jackknife cannot delete from a one-row stratum either.
  auto <- ht_total(s, "y")
  expect_true(is.na(auto$variance))
  expect_equal(auto$method, "none")
  expect_true(all(is.na(auto$ci)))
})

test_that("a stratum taken whole is not lonely", {
  d <- data.frame(g = c("solo", rep("big", 50)), y = c(1000, 1:50))
  s <- draw(d, design_stratified("g", n = 11, min_per_stratum = 1), seed = 2,
            weights = TRUE)
  # `solo` holds one row and it was taken: nothing random, nothing to estimate
  r <- ht_total(s, "y")
  expect_equal(r$method, "analytic")
  expect_true(is.finite(r$se))
})

test_that("one row per interval is refused the same way", {
  d <- data.frame(
    id = 1:120,
    ts = rep(seq(as.POSIXct("2024-01-01", tz = "UTC"), by = "day",
                 length.out = 12), each = 10),
    y = stats::rnorm(120, 50, 10)
  )
  des <- design_temporal("ts", from = "2024-01-01", to = "2024-01-13",
                         interval = 1, per_interval = 1, unit = "days")
  r <- ht_total(draw(d, des, seed = 1, weights = TRUE), "y")
  expect_true(is.na(r$variance))
  expect_match(r$note, "per_interval = 2")
})

# ---- stratified jackknife ----------------------------------------------------

test_that("the jackknife deletes within strata, and matches the analytic form", {
  set.seed(8)
  d <- data.frame(id = 1:400, g = rep(c("a", "b", "c", "d"), each = 100))
  d$y <- rep(c(10, 100, 300, 900), each = 100) + stats::rnorm(400, 0, 5)
  s <- draw(d, design_stratified("g", n = 40), seed = 2, weights = TRUE)

  a <- ht_total(s, "y", variance = "analytic")
  j <- ht_total(s, "y", variance = "jackknife")
  # Deleting across strata counted the gaps BETWEEN them as noise: the
  # standard error came out sixty times too large.
  expect_equal(j$variance, a$variance, tolerance = 1e-8)
  expect_match(j$note, "in 4 strata")

  # Same for the linearised Hajek mean
  expect_equal(ht_mean(s, "y", variance = "jackknife")$variance,
               ht_mean(s, "y", variance = "analytic")$variance,
               tolerance = 1e-8)
})

test_that("temporal intervals are jackknife strata too", {
  set.seed(9)
  d <- data.frame(
    ts = rep(seq(as.POSIXct("2024-01-01", tz = "UTC"), by = "day",
                 length.out = 6), each = 30),
    y = rep(c(5, 50, 500, 5, 50, 500), each = 30) + stats::rnorm(180)
  )
  des <- design_temporal("ts", from = "2024-01-01", to = "2024-01-07",
                         interval = 1, per_interval = 5, unit = "days")
  s <- draw(d, des, seed = 1, weights = TRUE)
  expect_equal(ht_total(s, "y", variance = "jackknife")$variance,
               ht_total(s, "y", variance = "analytic")$variance,
               tolerance = 1e-8)
})

# ---- Neyman allocation with capping ------------------------------------------

test_that("an over-allocated stratum is taken whole and the rest re-split", {
  # a would get 60 * 1000 / 1600 = 37.5 of its 10 rows. Take it whole; the
  # remaining 50 split b:c in proportion 200*1 : 200*2, as Cochran sets out.
  got <- drawn:::allocate(60, c(a = 10, b = 200, c = 200), "neyman", 0,
                          cap = TRUE, spread = c(a = 100, b = 1, c = 2))
  expect_equal(unname(got), c(10L, 17L, 33L))

  # Equal allocation caps the same way
  got <- drawn:::allocate(30, c(a = 2, b = 50, c = 50), "equal", 0, cap = TRUE)
  expect_equal(unname(got), c(2L, 14L, 14L))
})

test_that("a constant stratum under Neyman is flagged, not silently dropped", {
  d <- data.frame(g = rep(c("a", "b", "c"), each = 100),
                  x = c(rep(5, 100), stats::rnorm(200, 0, 10)))
  des <- design_stratified("g", n = 30, allocation = "neyman",
                           allocation_by = "x")
  expect_warning(draw(d, des, seed = 1), "`a`.*allocated no rows")
  # ...and the floor fixes it
  fixed <- design_stratified("g", n = 30, allocation = "neyman",
                             allocation_by = "x", min_per_stratum = 2)
  expect_true(all(tapply(inclusion_prob(d, fixed), d$g, min) > 0))
})

test_that("an empty stratum biases a total, which is why it warns", {
  d <- data.frame(g = rep(c("n", "s", "e", "w"), times = c(300, 180, 90, 30)))
  d$y <- ifelse(d$g == "w", 100, 1)
  des <- design_stratified("g", n = 10)
  expect_equal(unname(tapply(inclusion_prob(d, des), d$g, max)["w"]), 0)
  expect_warning(s <- draw(d, des, seed = 1, weights = TRUE),
                 "30 frame row")
  # The 3,000 the west stratum holds can never enter the estimate
  expect_lt(ht_total(s, "y", variance = "none")$total, 1000)
})

test_that("simulation loops do not repeat the empty-stratum warning", {
  d <- data.frame(g = c(rep("big", 99), "rare"))
  expect_silent(inclusion_prob(d, design_stratified("g", n = 5),
                               simulate = TRUE, R = 50, seed = 1))
})

# ---- degrees of freedom and t intervals -------------------------------------

test_that("the design's degrees of freedom are PSUs minus strata", {
  set.seed(3)
  d <- data.frame(id = 1:240, site = rep(c("a", "b", "c"), each = 80),
                  cl = rep(paste0("c", 1:24), each = 10),
                  y = stats::runif(240, 5, 200))
  df_of <- function(des) {
    ht_total(draw(d, des, seed = 1, weights = TRUE), "y")$df
  }
  expect_equal(df_of(design_simple(n = 30)), 29)
  expect_equal(df_of(design_stratified("site", n = 30)), 27)
  expect_equal(df_of(design_cluster("cl", n_clusters = 5)), 4)
  expect_equal(df_of(design_multistage("cl", n_clusters = 6, n = 18)), 5)
})

test_that("the interval uses t on those degrees of freedom, or normal on request", {
  set.seed(4)
  d <- data.frame(cl = rep(1:24, each = 10), y = stats::runif(240, 5, 200))
  s <- draw(d, design_cluster("cl", n_clusters = 4), seed = 1, weights = TRUE)
  r <- ht_total(s, "y")
  expect_equal(r$df, 3)
  expect_equal(diff(r$ci) / 2, stats::qt(0.975, 3) * r$se)

  nrm <- ht_total(s, "y", df = Inf)
  expect_equal(diff(nrm$ci) / 2, stats::qnorm(0.975) * nrm$se)
  expect_output(print(r), "t, 3 df")
  expect_error(ht_total(s, "y", df = -1), "`df`")
})

test_that("t intervals restore coverage with few clusters", {
  set.seed(5)
  d <- data.frame(id = 1:240, cl = rep(1:24, each = 10))
  d$y <- rep(stats::rnorm(24, 50, 20), each = 10) + stats::rnorm(240, 0, 5)
  truth <- sum(d$y)
  hits <- vapply(1:600, function(i) {
    r <- ht_total(draw(d, design_cluster("cl", n_clusters = 4), seed = i,
                       weights = TRUE), "y")
    c(t = r$ci[1] <= truth && truth <= r$ci[2],
      z = abs(r$total - truth) <= stats::qnorm(0.975) * r$se)
  }, logical(2))
  cover <- rowMeans(hits)
  # The normal interval covers about 87% at four clusters
  expect_lt(cover[["z"]], 0.92)
  expect_gt(cover[["t"]], 0.92)
})

test_that("degrees of freedom agree with survey::degf()", {
  skip_if_not_installed("survey")
  set.seed(6)
  d <- data.frame(id = 1:240, site = rep(c("a", "b", "c"), each = 80),
                  cl = rep(paste0("c", 1:24), each = 10),
                  y = stats::runif(240, 5, 200))
  for (des in list(design_simple(n = 30), design_stratified("site", n = 30),
                   design_cluster("cl", n_clusters = 5),
                   design_multistage("cl", n_clusters = 6, n = 18))) {
    s <- draw(d, des, seed = 1, weights = TRUE)
    expect_equal(ht_total(s, "y")$df, survey::degf(as_svydesign(s)))
  }
})

# ---- Deville for systematic PPS ----------------------------------------------

test_that("systematic PPS uses Deville's approximation, exactly as published", {
  set.seed(5)
  pop <- data.frame(x = stats::rgamma(300, 2, 0.1))
  pop$y <- 3 * pop$x + stats::rnorm(300, 0, 15)
  s <- draw(pop, design_weighted("x", n = 60, method = "systematic"), seed = 1,
            weights = TRUE)
  r <- ht_total(s, "y")
  expect_equal(r$method, "deville")
  expect_output(print(r), "Deville approximation")

  p <- s$.prob; u <- s$y / p; c_k <- 1 - p; a <- c_k / sum(c_k)
  expect_equal(r$variance,
               sum(c_k * (u - sum(a * u))^2) / (1 - sum(a^2)))
  skip_if_not_installed("sampling")
  expect_equal(r$variance, sampling::varest(s$y, pik = s$.prob))
})

test_that("Deville tracks the real sampling variance of systematic PPS", {
  set.seed(5)
  pop <- data.frame(x = stats::rgamma(300, 2, 0.1))
  pop$y <- 3 * pop$x + stats::rnorm(300, 0, 15)
  d <- design_weighted("x", n = 60, method = "systematic")
  out <- vapply(1:600, function(i) {
    r <- ht_total(draw(pop, d, seed = i, weights = TRUE), "y")
    c(r$total, r$variance)
  }, numeric(2))
  ratio <- mean(out[2, ]) / stats::var(out[1, ])
  expect_gt(ratio, 0.85)
  expect_lt(ratio, 1.15)
})

test_that("certainty rows add nothing to Deville's approximation", {
  pop <- data.frame(x = c(1000, 900, stats::runif(98, 1, 10)))
  pop$y <- pop$x * 2 + stats::rnorm(100)
  s <- draw(pop, design_weighted("x", n = 10, method = "systematic"), seed = 3,
            weights = TRUE)
  expect_true(any(s$.prob == 1))
  r <- ht_total(s, "y")
  free <- s$.prob < 1
  u <- s$y[free] / s$.prob[free]; c_k <- 1 - s$.prob[free]; a <- c_k / sum(c_k)
  expect_equal(r$variance, sum(c_k * (u - sum(a * u))^2) / (1 - sum(a^2)))
})

# ---- successive differences for systematic ---------------------------------

test_that("successive differences follow the walk order, not frame order", {
  # A frame sorted on something related to y: the setting systematic sampling
  # is chosen for. Each population's true variance is computed exactly, by
  # enumerating every start, and the ratio is pooled over many populations.
  set.seed(12)
  N <- 400; k <- 10
  ratios <- vapply(1:40, function(r) {
    d <- data.frame(x = stats::runif(N, 0, 100))
    d$y <- d$x + stats::rnorm(N, 0, 20)
    out <- vapply(seq_len(k), function(st) {
      s <- draw(d, design_systematic(interval = k, order_by = "x"),
                seed = st * 1000 + r, weights = TRUE)
      c(ht_total(s, "y")$total, ht_total(s, "y")$variance,
        N^2 * (1 - 1 / k) * stats::var(s$y) / nrow(s))
    }, numeric(3))
    starts <- vapply(seq_len(k), function(st) {
      s <- draw(d, design_systematic(interval = k, order_by = "x", start = st))
      k * sum(s$y)
    }, numeric(1))
    truth <- mean((starts - mean(starts))^2)
    c(sd = mean(out[2, ]), srs = mean(out[3, ]), truth = truth)
  }, numeric(3))
  pooled <- rowMeans(ratios)
  ratio_sd <- pooled[["sd"]] / pooled[["truth"]]
  ratio_srs <- pooled[["srs"]] / pooled[["truth"]]
  # The simple random formula cannot see the ordering and overstates; the
  # successive-difference form can.
  expect_gt(ratio_srs, 1.8)
  expect_gt(ratio_sd, 0.7)
  expect_lt(ratio_sd, 1.4)
})

test_that("successive differences approximate the SRS figure on a shuffled frame", {
  set.seed(13)
  d <- data.frame(y = stats::rnorm(500, 100, 20))
  out <- vapply(1:500, function(i) {
    r <- ht_total(draw(d, design_systematic(interval = 10), seed = i,
                       weights = TRUE), "y")
    c(r$total, r$variance)
  }, numeric(2))
  ratio <- mean(out[2, ]) / stats::var(out[1, ])
  expect_gt(ratio, 0.75)
  expect_lt(ratio, 1.35)
})

# ---- one row per cluster -----------------------------------------------------

test_that("one row per cluster falls back to the ultimate-cluster jackknife", {
  set.seed(3)
  d <- data.frame(cl = rep(paste0("c", 1:24), each = 10))
  d$y <- rep(stats::rnorm(24, 50, 5), each = 10) + stats::rnorm(240, 0, 20)
  des <- design_multistage("cl", n_clusters = 8, n = 8)

  expect_true(is.na(ht_total(draw(d, des, seed = 1, weights = TRUE), "y",
                             variance = "analytic")$variance))
  out <- vapply(1:800, function(i) {
    r <- ht_total(draw(d, des, seed = i, weights = TRUE), "y")
    c(r$total, r$variance, r$method == "jackknife")
  }, numeric(3))
  expect_true(all(out[3, ] == 1))
  # The exact two-stage form ran at 0.70 here; the jackknife errs the safe way
  ratio <- mean(out[2, ]) / stats::var(out[1, ])
  expect_gt(ratio, 0.9)
  expect_lt(ratio, 1.6)
})

# ---- printing ----------------------------------------------------------------

test_that("the HT mean is labelled as such, not as Hajek", {
  d <- data.frame(y = stats::runif(100))
  s <- draw(d, design_simple(n = 20), seed = 1, weights = TRUE)
  expect_output(print(ht_mean(s, "y", estimator = "ht")),
                "^Horvitz-Thompson mean")
  expect_output(print(ht_mean(s, "y")), "^Hajek mean")
  expect_equal(ht_mean(s, "y", estimator = "ht")$estimator, "ht")
})

test_that("a logical y estimates a proportion", {
  d <- data.frame(y = rep(c(TRUE, FALSE), times = c(30, 70)))
  s <- draw(d, design_simple(n = 40), seed = 1, weights = TRUE)
  expect_equal(ht_mean(s, "y")$mean, mean(s$y))
})
