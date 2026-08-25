test_that("closed forms match the textbook values", {
  d1 <- data.frame(id = 1:10)
  expect_equal(inclusion_prob(d1, design_simple(n = 3)), rep(0.3, 10))

  d2 <- data.frame(id = 1:20, g = rep(c("a", "b"), times = c(15, 5)))
  # n = 8 proportional over 15/5 allocates 6 and 2
  expect_equal(inclusion_prob(d2, design_stratified("g", n = 8)),
               c(rep(6/15, 15), rep(2/5, 5)))

  d4 <- data.frame(id = 1:20, cl = rep(letters[1:5], each = 4))
  expect_equal(inclusion_prob(d4, design_cluster("cl", n_clusters = 2)),
               rep(0.4, 20))

  expect_equal(inclusion_prob(data.frame(id = 1:12), design_reservoir(n = 4)),
               rep(1/3, 12))
})

test_that("systematic inclusion is 1/interval for every row, at any N", {
  for (N in c(20, 23, 29)) {
    p <- inclusion_prob(data.frame(id = seq_len(N)),
                        design_systematic(interval = 7))
    expect_equal(p, rep(1/7, N), info = paste("N =", N))
  }
})

test_that("a fixed systematic start is refused as a probability sample", {
  expect_error(
    inclusion_prob(data.frame(id = 1:20), design_systematic(interval = 5, start = 2)),
    "not a probability sample"
  )
})

test_that("sampling_weight is the reciprocal, and 0 becomes NA not Inf", {
  d <- data.frame(id = 1:10)
  expect_equal(sampling_weight(d, design_simple(n = 5)), rep(2, 10))

  # A row outside the window has probability 0 and no finite weight.
  ts <- data.frame(id = 1:6,
                   t = seq(as.POSIXct("2020-01-01", tz = "UTC"), by = "hour",
                           length.out = 6))
  w <- sampling_weight(ts, design_temporal("t", "2020-01-01", "2020-01-01 03:00:00",
                                         interval = 1, per_interval = 1,
                                         unit = "hours"))
  expect_true(all(is.na(w[4:6])))
  expect_true(all(is.finite(w[1:3])))
})

test_that("rows the design can never reach get probability 0", {
  ts <- data.frame(id = 1:6,
                   t = seq(as.POSIXct("2020-01-01", tz = "UTC"), by = "hour",
                           length.out = 6))
  p <- inclusion_prob(ts, design_temporal("t", "2020-01-01",
                                          "2020-01-01 03:00:00", interval = 1,
                                          per_interval = 1, unit = "hours"))
  expect_equal(p, c(1, 1, 1, 0, 0, 0))
})

test_that("designs with no closed form refuse rather than invent one", {
  df <- data.frame(id = 1:20, cl = rep(letters[1:4], each = 5), w = 1:20)

  expect_error(inclusion_prob(df, design_weighted("w", n = 3)),
               "no closed-form inclusion probability")
  expect_error(
    inclusion_prob(df, design_cluster("cl", n_clusters = 2, balanced = TRUE)),
    "no closed-form inclusion probability"
  )
  expect_error(
    inclusion_prob(df, design_multistage("cl", n_clusters = 2, n = 4,
                                         allocation = "proportional")),
    "no closed-form inclusion probability"
  )
  expect_error(inclusion_prob(df, design_bootstrap(n_replicates = 5)),
               "not a probability sample of a finite population")
})

test_that("simulate = TRUE works where the closed form refuses", {
  df <- data.frame(id = 1:12, cl = rep(letters[1:4], each = 3))
  p <- inclusion_prob(df, design_cluster("cl", n_clusters = 2, balanced = TRUE),
                      simulate = TRUE, R = 400, seed = 1)
  expect_length(p, 12)
  expect_true(all(p >= 0 & p <= 1))
  # Whole clusters of equal size: every row shares its cluster's rate.
  expect_equal(mean(p), 0.5, tolerance = 0.1)
})

test_that("simulation agrees with the closed form where both exist", {
  d <- data.frame(id = 1:20, g = rep(c("a", "b"), times = c(15, 5)))
  design <- design_stratified("g", n = 8)
  exact <- inclusion_prob(d, design)
  sim <- inclusion_prob(d, design, simulate = TRUE, R = 3000, seed = 1)
  expect_equal(sim, exact, tolerance = 0.06)
})

test_that("inclusion_prob validates its arguments", {
  df <- data.frame(id = 1:10)
  expect_error(inclusion_prob(df, "simple"), "must come from one of the design")
  expect_error(inclusion_prob(df, design_simple(n = 2), simulate = "yes"),
               "TRUE or FALSE")
  expect_error(inclusion_prob(df, design_simple(n = 2), simulate = TRUE, R = 0),
               "at least 1")
})
