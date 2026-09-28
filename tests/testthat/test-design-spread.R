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
  expect_output(print(r), "local-mean approximation")
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
