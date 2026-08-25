est_pop <- function() {
  set.seed(404)
  N <- 300
  d <- data.frame(
    id   = 1:N,
    site = rep(c("a", "b", "c"), times = c(150, 100, 50)),
    cl   = rep(paste0("c", 1:30), each = 10),
    y    = round(stats::runif(N, 20, 400))
  )
  d$w <- d$y + 100
  d
}

test_that("ht_mean recovers the population mean under a range of designs", {
  d <- est_pop()
  truth <- mean(d$y)

  for (des in list(
    design_simple(n = 80),
    design_stratified("site", n = 80),
    design_weighted("w", n = 80, method = "poisson")
  )) {
    est <- replicate(400, {
      s <- draw(d, des, weights = TRUE)
      ht_mean(s, "y", variance = "none")$mean
    })
    se <- stats::sd(est) / sqrt(length(est))
    expect_lt(abs(mean(est) - truth), 4 * se)
  }
})

test_that("the two estimators differ, and Hajek is the steadier one", {
  d <- est_pop()
  # Unequal probabilities are what separate them; with equal ones they agree.
  eq <- draw(d, design_simple(n = 60), seed = 3, weights = TRUE)
  expect_equal(ht_mean(eq, "y", "hajek", variance = "none")$mean,
               ht_mean(eq, "y", "ht", variance = "none")$mean)

  des <- design_weighted("w", n = 60, method = "poisson")
  pair <- vapply(1:200, function(i) {
    s <- draw(d, des, seed = i, weights = TRUE)
    c(hajek = ht_mean(s, "y", "hajek", variance = "none")$mean,
      ht    = ht_mean(s, "y", "ht", variance = "none")$mean)
  }, numeric(2))
  expect_lt(stats::sd(pair["hajek", ]), stats::sd(pair["ht", ]))
})

test_that("the Hajek variance tracks the actual sampling variance", {
  d <- est_pop()
  des <- design_stratified("site", n = 90)
  out <- vapply(1:400, function(i) {
    s <- draw(d, des, seed = i, weights = TRUE)
    e <- ht_mean(s, "y")
    c(e$mean, e$variance)
  }, numeric(2))
  empirical <- stats::var(out[1, ])
  expect_gt(mean(out[2, ]) / empirical, 0.75)
  expect_lt(mean(out[2, ]) / empirical, 1.35)
})

test_that("the interval widens with the level and is centred on the estimate", {
  d <- est_pop()
  s <- draw(d, design_stratified("site", n = 60), seed = 11, weights = TRUE)
  a <- ht_mean(s, "y", level = 0.90)
  b <- ht_mean(s, "y", level = 0.99)
  expect_equal(a$mean, b$mean)
  expect_lt(diff(a$ci), diff(b$ci))
  expect_equal(mean(a$ci), a$mean)
})

test_that("ht_mean accepts a vector as well as a column name", {
  d <- est_pop()
  s <- draw(d, design_simple(n = 50), seed = 2, weights = TRUE)
  expect_equal(ht_mean(s, "y", variance = "none")$mean,
               ht_mean(s, s$y, variance = "none")$mean)
})

test_that("ht_mean reports why it has no variance rather than inventing one", {
  d <- est_pop()
  s <- draw(d, design_systematic(interval = 5), seed = 1, weights = TRUE)
  e <- ht_mean(s, "y")
  expect_true(is.finite(e$mean))
  expect_true(is.na(e$se))
  expect_true(nzchar(e$note))
  expect_output(print(e), "se\\s+NA")
})

test_that("ht_mean refuses samples it cannot estimate from", {
  d <- est_pop()
  bare <- draw(d, design_simple(n = 20), seed = 1)
  expect_error(ht_mean(bare, "y"), "weights = TRUE")

  s <- draw(d, design_simple(n = 20), seed = 1, weights = TRUE)
  expect_error(ht_mean(s, "nope"), "no column")
  expect_error(ht_mean(s, "site"), "numeric")
  expect_error(ht_mean(s, "y", level = 1), "between 0 and 1")

  na_y <- s
  na_y$y[1] <- NA
  expect_error(ht_mean(na_y, "y"), "missing")
})

test_that("deff falls below 1 when stratification helps and rises when it hurts", {
  set.seed(88)
  N <- 400
  d <- data.frame(
    id   = 1:N,
    site = rep(c("a", "b", "c", "d"), each = 100),
    cl   = rep(paste0("c", 1:40), each = 10)
  )
  # The grouping variable is the whole story of y, so stratifying on it removes
  # nearly all the variance and clustering on it keeps nearly all of it.
  d$y <- rep(c(20, 80, 160, 260), each = 100) + round(stats::rnorm(N, 0, 6))

  st <- deff(ht_total(draw(d, design_stratified("site", n = 40), seed = 1,
                           weights = TRUE), "y"))
  cl <- deff(ht_total(draw(d, design_cluster("cl", n_clusters = 4), seed = 1,
                           weights = TRUE), "y"))
  expect_lt(st, 0.5)
  expect_gt(cl, 2)
})

test_that("deff is available from both estimators and NA without a variance", {
  d <- est_pop()
  s <- draw(d, design_stratified("site", n = 60), seed = 5, weights = TRUE)
  expect_true(is.finite(deff(ht_total(s, "y"))))
  expect_true(is.finite(deff(ht_mean(s, "y"))))
  expect_true(is.na(deff(ht_mean(s, "y", variance = "none"))))
  expect_error(deff(1:10), "ht_total")
})

test_that("the jackknife deletes clusters, not rows", {
  # The primary sampling unit of a cluster design is the cluster. Deleting one
  # row at a time instead treats 80 correlated rows as 80 independent ones and
  # inflates the variance; the note names the count, so it is checkable.
  d <- est_pop()
  s <- draw(d, design_cluster("cl", n_clusters = 8), seed = 1, weights = TRUE)
  expect_equal(nrow(s), 80)

  for (r in list(ht_total(s, "y", variance = "jackknife"),
                 ht_mean(s, "y", variance = "jackknife"))) {
    expect_match(r$note, "Jackknife over 8 primary")
    expect_false(grepl("over 80 primary", r$note))
  }

  # And it lands on the analytic answer, which is the point of grouping right.
  expect_equal(ht_total(s, "y", variance = "jackknife")$variance,
               ht_total(s, "y", variance = "analytic")$variance,
               tolerance = 1e-6)
})

test_that("ht_mean falls back to the jackknife and says it did", {
  d <- est_pop()
  s <- draw(d, design_weighted("w", n = 60, method = "systematic"), seed = 1,
            weights = TRUE)
  e <- ht_mean(s, "y")
  expect_equal(e$method, "jackknife")
  expect_true(is.finite(e$se))
  expect_match(e$note, "jackknife was used")
})

test_that("printing an estimate shows the pieces someone would report", {
  d <- est_pop()
  s <- draw(d, design_stratified("site", n = 60), seed = 5, weights = TRUE)
  out <- utils::capture.output(print(ht_mean(s, "y")))
  expect_true(any(grepl("mean", out)))
  expect_true(any(grepl("estimate", out)))
  expect_true(any(grepl("95% CI", out, fixed = TRUE)))
  expect_true(any(grepl("deff", out)))
})

test_that("deff is exactly 1 for a simple random sample", {
  # The reference variance divided an n-term sum by N - 1, leaving it short by
  # n(N-1)/(N(n-1)) -- so deff came out at 1 + 1/n and a plain SRS of 10
  # printed as "worse than simple random sampling".
  set.seed(31)
  pop <- data.frame(y = stats::rnorm(5000, 100, 20))
  for (n in c(10, 20, 50, 200)) {
    got <- vapply(1:200, function(i) {
      deff(ht_total(draw(pop, design_simple(n = n), seed = i, weights = TRUE),
                    "y"))
    }, numeric(1))
    expect_equal(mean(got), 1, tolerance = 1e-9,
                 label = paste("deff at n =", n))
  }
})

test_that("deff is NA rather than negative when there is no variance", {
  set.seed(5)
  pop <- data.frame(y = stats::rnorm(60, 10, 3))
  r <- ht_total(draw(pop, design_systematic(interval = 5), seed = 1,
                     weights = TRUE), "y")
  expect_true(is.na(deff(r)))
  expect_true(is.na(drawn:::deff_value(-5, 1:10, rep(0.5, 10), 100, 10,
                                       "total")))
})
