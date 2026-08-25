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
