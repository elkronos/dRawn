var_pop <- function() {
  set.seed(7)
  N <- 240
  d <- data.frame(
    id   = 1:N,
    site = rep(c("a", "b", "c"), times = c(120, 80, 40)),
    cl   = rep(paste0("c", 1:24), each = 10),
    y    = round(stats::runif(N, 5, 200))
  )
  d$w <- d$y + 50
  d$ts <- rep(seq(as.POSIXct("2024-01-01", tz = "UTC"), by = "day",
                  length.out = 12), each = 20)
  d
}

test_that("joint probabilities match the textbook forms", {
  d <- data.frame(id = 1:10)
  m <- joint_prob(d, design_simple(n = 4))
  expect_equal(dim(m), c(10L, 10L))
  expect_equal(diag(m), rep(0.4, 10))
  expect_equal(m[1, 2], 4 * 3 / (10 * 9))

  # Poisson: rows are independent, so the joint is the product
  w <- data.frame(id = 1:8, w = c(1, 1, 2, 2, 3, 3, 4, 4))
  pd <- design_weighted("w", n = 3, method = "poisson")
  pj <- joint_prob(w, pd)
  pi <- inclusion_prob(w, pd)
  expect_equal(pj[1, 2], pi[1] * pi[2])
})

test_that("rows selects a submatrix rather than the whole square", {
  d <- data.frame(id = 1:1000)
  m <- joint_prob(d, design_simple(n = 50), rows = c(3, 17, 400))
  expect_equal(dim(m), c(3L, 3L))
  expect_equal(diag(m), rep(50 / 1000, 3))
})

test_that("cluster joints differ within and between clusters", {
  d <- data.frame(id = 1:20, cl = rep(letters[1:5], each = 4))
  m <- joint_prob(d, design_cluster("cl", n_clusters = 2))
  expect_equal(m[1, 2], 2 / 5)            # same cluster: that cluster selected
  expect_equal(m[1, 5], 2 * 1 / (5 * 4))  # different clusters: both selected
})

test_that("systematic joints are zero for rows that cannot co-occur", {
  d <- data.frame(id = 1:20)
  m <- joint_prob(d, design_systematic(interval = 5))
  expect_equal(m[1, 6], 1 / 5)   # same residue class
  expect_equal(m[1, 2], 0)       # different residue class: impossible together
})

test_that("designs without a closed form refuse", {
  d <- var_pop()
  expect_error(joint_prob(d, design_weighted("w", n = 10, method = "systematic")),
               "sampling::UPsystematicpi2")
  expect_error(joint_prob(d, design_weighted("w", n = 10)), "no closed-form")
  expect_error(joint_prob(d, design_cluster("cl", n_clusters = 3, balanced = TRUE)),
               "no closed-form")
  expect_error(joint_prob(d, design_bootstrap(n_replicates = 5)),
               "not a probability sample")
  expect_error(joint_prob(d, design_multistage("cl", n_clusters = 4, n = 10)),
               "does not divide|not constant")
})

test_that("ht_total recovers the total and needs weights = TRUE", {
  d <- var_pop()
  s <- draw(d, design_stratified("site", n = 40), seed = 1, weights = TRUE)
  r <- ht_total(s, "y")
  expect_s3_class(r, "drawn_ht")
  expect_equal(r$n, 40)
  expect_true(is.finite(r$se))
  expect_true(r$ci[1] < r$total && r$total < r$ci[2])

  plain <- draw(d, design_stratified("site", n = 40), seed = 1)
  expect_error(ht_total(plain, "y"), "does not carry its design")
})

test_that("ht_total validates y and level", {
  d <- var_pop()
  s <- draw(d, design_simple(n = 30), seed = 1, weights = TRUE)
  expect_error(ht_total(s, "nope"), "no column")
  expect_error(ht_total(s, letters[1:30]), "numeric column name")
  expect_error(ht_total(s, "y", level = 1.5), "between 0 and 1")
})

test_that("the estimated variance matches the estimator's true variance", {
  d <- var_pop()
  designs <- list(
    simple     = design_simple(n = 30),
    stratified = design_stratified("site", n = 30),
    cluster    = design_cluster("cl", n_clusters = 6),
    multistage = design_multistage("cl", n_clusters = 8, n = 32),
    poisson    = design_weighted("w", n = 30, method = "poisson")
  )
  R <- 700
  for (nm in names(designs)) {
    est <- numeric(R); vr <- numeric(R)
    for (i in seq_len(R)) {
      s <- draw(d, designs[[nm]], seed = i, weights = TRUE)
      if (nrow(s) < 2) { est[i] <- NA; vr[i] <- NA; next }
      r <- ht_total(s, "y")
      est[i] <- r$total; vr[i] <- r$variance
    }
    ok <- !is.na(est) & !is.na(vr)
    ratio <- mean(vr[ok]) / stats::var(est[ok])
    expect_gt(ratio, 0.75, label = paste(nm, "variance ratio"))
    expect_lt(ratio, 1.35, label = paste(nm, "variance ratio"))
  }
})

test_that("systematic sampling reports no variance and says why", {
  d <- var_pop()
  s <- draw(d, design_systematic(interval = 8), seed = 1, weights = TRUE)
  r <- ht_total(s, "y")
  expect_true(is.finite(r$total))
  expect_true(is.na(r$variance))
  expect_match(r$note, "no design-unbiased variance")
  expect_true(all(is.na(r$ci)))
})

test_that("the print method shows the interval, or the reason there isn't one", {
  d <- var_pop()
  out <- capture.output(print(ht_total(draw(d, design_simple(n = 30), seed = 1,
                                            weights = TRUE), "y")))
  expect_true(any(grepl("Horvitz-Thompson total", out)))
  expect_true(any(grepl("95% CI", out)))

  out2 <- capture.output(print(ht_total(draw(d, design_systematic(interval = 8),
                                             seed = 1, weights = TRUE), "y")))
  expect_true(any(grepl("se\\s+NA", out2)))
})

test_that("the jackknife supplies a variance where the analytic form cannot", {
  set.seed(3); N <- 240
  pop <- data.frame(id = 1:N, cl = rep(paste0("c", 1:24), each = 10),
                    y = round(stats::runif(N, 5, 200)))
  pop$w <- pop$y + 50

  # Systematic PPS has first-order probabilities but no closed-form joint ones
  r <- ht_total(draw(pop, design_weighted("w", n = 30, method = "systematic"),
                     seed = 1, weights = TRUE), "y")
  expect_equal(r$method, "jackknife")
  expect_true(is.finite(r$se))
  expect_match(r$note, "jackknife was used")
})

test_that("jackknife and analytic agree where both exist", {
  set.seed(3); N <- 240
  pop <- data.frame(id = 1:N, cl = rep(paste0("c", 1:24), each = 10),
                    y = round(stats::runif(N, 5, 200)))
  for (d in list(design_simple(n = 30), design_cluster("cl", n_clusters = 6))) {
    s <- draw(pop, d, seed = 1, weights = TRUE)
    a <- ht_total(s, "y", variance = "analytic")$variance
    j <- ht_total(s, "y", variance = "jackknife")$variance
    expect_equal(j, a, tolerance = 1e-6)
  }
})

test_that("systematic sampling refuses the jackknife too", {
  set.seed(3)
  pop <- data.frame(id = 1:240, y = round(stats::runif(240, 5, 200)))
  r <- ht_total(draw(pop, design_systematic(interval = 8), seed = 1,
                     weights = TRUE), "y", variance = "jackknife")
  expect_true(is.na(r$variance))
  expect_match(r$note, "one primary sampling unit")
})

test_that("variance = 'none' and 'analytic' behave as documented", {
  pop <- data.frame(id = 1:100, y = 1:100)
  s <- draw(pop, design_simple(n = 20), seed = 1, weights = TRUE)
  expect_true(is.na(ht_total(s, "y", variance = "none")$variance))
  expect_equal(ht_total(s, "y", variance = "analytic")$method, "analytic")
  expect_error(ht_total(s, "y", variance = "nope"))
})

# ---- simulated joint probabilities -----------------------------------------

test_that("simulated joints converge on the closed form where one exists", {
  d <- data.frame(id = 1:20, site = rep(c("a", "b"), each = 10))
  for (des in list(design_simple(n = 6),
                   design_stratified("site", n = 6),
                   design_systematic(interval = 4))) {
    exact <- joint_prob(d, des)
    sim <- joint_prob(d, des, simulate = TRUE, R = 4000, seed = 1)
    expect_equal(dim(sim), dim(exact))
    expect_lt(max(abs(sim - exact)), 0.035)
  }
})

test_that("simulation supplies joints for designs with no closed form", {
  set.seed(5)
  d <- data.frame(id = 1:30, w = stats::runif(30, 1, 10))
  des <- design_weighted("w", n = 8, method = "systematic")
  expect_error(joint_prob(d, des), "UPsystematicpi2")

  m <- joint_prob(d, des, simulate = TRUE, R = 3000, seed = 2)
  expect_equal(dim(m), c(30L, 30L))
  expect_equal(m, t(m))
  # The diagonal is the first-order probability, which this design does have.
  expect_lt(max(abs(diag(m) - inclusion_prob(d, des))), 0.04)
  # No pair can be likelier than either row on its own.
  expect_true(all(m <= outer(diag(m), diag(m), pmin) + 1e-9))
})

test_that("simulation respects `rows` and is reproducible", {
  d <- data.frame(id = 1:40)
  des <- design_simple(n = 10)
  a <- joint_prob(d, des, rows = c(2, 9, 30), simulate = TRUE, R = 500,
                  seed = 42)
  b <- joint_prob(d, des, rows = c(2, 9, 30), simulate = TRUE, R = 500,
                  seed = 42)
  expect_equal(dim(a), c(3L, 3L))
  expect_identical(a, b)
  expect_false(identical(
    a, joint_prob(d, des, rows = c(2, 9, 30), simulate = TRUE, R = 500,
                  seed = 43)))
})

test_that("simulation validates R and refuses a reserved column name", {
  d <- data.frame(id = 1:20)
  expect_error(joint_prob(d, design_simple(n = 5), simulate = TRUE, R = 0),
               "R")
  clash <- data.frame(id = 1:20, .drawn_row_id = 1:20)
  expect_error(joint_prob(clash, design_simple(n = 5), simulate = TRUE,
                          R = 10),
               "already has a column")
})

# ---- the sampling unit has to be the one the design actually uses ----------

test_that("cluster variance is honest when clusters differ in size", {
  # Whole clusters are taken, so the row count is random unless every cluster
  # is the same size. Sen-Yates-Grundy assumes a fixed size; applied row by row
  # here it came out at 0.39 of the truth, went negative on 57% of samples, and
  # produced zero-width intervals.
  set.seed(9)
  sizes <- c(3, 4, 5, 6, 3, 4, 5, 6, 7, 8, 4, 5)
  pop <- data.frame(cl = rep(paste0("c", seq_along(sizes)), times = sizes))
  pop$y <- round(stats::rnorm(nrow(pop), 100, 30) +
                   rep(stats::rnorm(length(sizes), 0, 40), times = sizes))
  des <- design_cluster("cl", n_clusters = 4)

  out <- vapply(1:1500, function(i) {
    e <- ht_total(draw(pop, des, seed = i, weights = TRUE), "y")
    c(e$total, e$variance)
  }, numeric(2))

  expect_true(all(out[2, ] >= 0))
  ratio <- mean(out[2, ]) / stats::var(out[1, ])
  expect_gt(ratio, 0.85)
  expect_lt(ratio, 1.2)

  # The cluster-level form is algebraically the delete-a-cluster jackknife.
  s <- draw(pop, des, seed = 1, weights = TRUE)
  expect_equal(ht_total(s, "y", variance = "analytic")$variance,
               ht_total(s, "y", variance = "jackknife")$variance,
               tolerance = 1e-8)

  # A design where every cluster total is identical has no variance at all,
  # and one where they differ must not report zero.
  flat <- data.frame(cl = rep(c("a", "b", "c", "d"), each = 5), y = 1)
  expect_equal(ht_total(draw(flat, design_cluster("cl", n_clusters = 2),
                             seed = 1, weights = TRUE), "y")$variance, 0)
})

test_that("the jackknife declines where deleting rows misrepresents the design", {
  set.seed(4)
  pop <- data.frame(w = stats::runif(60, 1, 20), y = round(stats::rnorm(60, 100, 25)))

  # Poisson: rows are independent and the exact form is right there.
  s <- draw(pop, design_weighted("w", n = 12, method = "poisson"), seed = 1,
            weights = TRUE)
  jk <- ht_total(s, "y", variance = "jackknife")
  expect_true(is.na(jk$variance))
  expect_match(jk$note, "variance = \"analytic\"")
  expect_true(is.finite(ht_total(s, "y", variance = "analytic")$variance))

  # Systematic PPS is still a fixed-size draw without replacement, so the
  # finite population correction belongs and is applied. Dropping it overstates
  # the variance badly once the sampling fraction is large -- 1.77 of the truth
  # at n/N = 0.375, against 1.11 with the correction.
  ps <- draw(pop, design_weighted("w", n = 12, method = "systematic"), seed = 1,
             weights = TRUE)
  r <- ht_total(ps, "y")
  expect_equal(r$method, "jackknife")
  expect_match(r$note, "finite population correction")
})

test_that("the systematic-PPS jackknife tracks the real sampling variance", {
  set.seed(19)
  N <- 300
  pop <- data.frame(w = stats::rlnorm(N, 5, 0.9))
  # A response variable distinct from the size measure. Totalling the size
  # measure itself makes y/pi constant and every variance exactly zero, which
  # would let any estimator pass.
  pop$y <- round(pop$w * 0.3 + stats::rnorm(N, 0, 40))

  for (n in c(30, 90)) {
    des <- design_weighted("w", n = n, method = "systematic")
    o <- vapply(1:1200, function(i) {
      e <- ht_total(draw(pop, des, seed = i, weights = TRUE), "y")
      c(e$total, e$variance)
    }, numeric(2))
    ratio <- mean(o[2, ]) / stats::var(o[1, ])
    expect_gt(ratio, 0.75, label = paste("pps-systematic ratio at n =", n))
    expect_lt(ratio, 1.35, label = paste("pps-systematic ratio at n =", n))
  }
})

test_that("the jackknife treats certainty rows as fixed, not as sampling units", {
  # Certainty rows are in every possible sample. Deleting one and inflating the
  # rest invents variance the design does not have -- 4.3x too much.
  set.seed(12)
  pop <- data.frame(w = c(rep(100, 6), rep(1, 54)),
                    y = round(stats::rnorm(60, 100, 25), 3))
  des <- design_certainty("w", 50, design_simple(n = 10))
  s <- draw(pop, des, seed = 1, weights = TRUE)

  jk <- ht_total(s, "y", variance = "jackknife")
  expect_match(jk$note, "Jackknife over 10 primary")
  expect_false(grepl("over 16 primary", jk$note))
  expect_equal(jk$variance, ht_total(s, "y", variance = "analytic")$variance,
               tolerance = 1e-8)
})

test_that("a negative analytic variance is reported, not passed on", {
  # Sen-Yates-Grundy can go negative on an unlucky sample. sqrt() was already
  # guarded, but the reason was dropped, `auto` did not fall back, and deff()
  # divided by it and returned a negative design effect.
  fake <- list(variance = -1234, method = "analytic", note = NULL)
  expect_true(is.na(drawn:::deff_value(-1234, 1:10, rep(0.5, 10), 100, 10,
                                       "total")))

  pop <- data.frame(id = 1:60, y = 1:60)
  s <- draw(pop, design_simple(n = 20), seed = 1, weights = TRUE)
  bad <- drawn:::ht_variance_dispatch(
    structure(list(), class = c("drawn_design_nope", "drawn_design")),
    s, pop, attr(s, "drawn_rows"), s$y, s$.prob, "analytic")
  expect_true(is.na(bad$variance))
  expect_true(nzchar(bad$note))
})

test_that("a design neither estimator can handle says so once, truthfully", {
  set.seed(3)
  pop <- data.frame(y = stats::rnorm(60, 10, 3))
  r <- ht_total(draw(pop, design_systematic(interval = 5), seed = 1,
                     weights = TRUE), "y")
  # Nothing computed a variance, so nothing may claim to have been used.
  expect_true(is.na(r$variance))
  expect_equal(r$method, "none")
  expect_match(r$note, "no design-unbiased variance")
  expect_match(r$note, "one primary sampling unit")
})
