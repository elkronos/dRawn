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

# ---- probabilities must describe the draw that actually happens ------------

test_that("multistage refuses where the per-cluster take is not constant", {
  # allocate() runs at draw time over the SELECTED clusters, so a remainder
  # goes to the largest of them and an undersized cluster is capped. Averaging
  # over the possibilities was out by 13% and could exceed 1.
  uneven <- data.frame(cl = rep(c("A", "B", "C", "D"), c(20, 10, 8, 6)))
  expect_error(inclusion_prob(uneven, design_multistage("cl", n_clusters = 2,
                                                        n = 7)),
               "no closed-form")
  tiny <- data.frame(cl = rep(c("x", "b", "c", "e"), c(1, 10, 10, 10)))
  expect_error(inclusion_prob(tiny, design_multistage("cl", n_clusters = 2,
                                                      n = 8)),
               "fewer than")

  # Where the take really is constant, it is exact.
  even <- data.frame(cl = rep(paste0("c", 1:8), each = 10))
  des <- design_multistage("cl", n_clusters = 4, n = 20)
  p <- inclusion_prob(even, des)
  expect_equal(unique(p), 0.25)
  sim <- inclusion_prob(even, des, simulate = TRUE, R = 4000, seed = 1)
  expect_lt(max(abs(sim - p)), 0.03)
})

test_that("max_items bounds what a reservoir design can reach", {
  d <- data.frame(y = 1:20)
  des <- design_reservoir(n = 4, max_items = 10)
  p <- inclusion_prob(d, des)
  expect_equal(p, c(rep(0.4, 10), rep(0, 10)))
  expect_equal(sum(p), 4)

  m <- joint_prob(d, des, rows = c(1, 2, 15))
  expect_equal(diag(m), c(0.4, 0.4, 0))
  expect_equal(m[1, 3], 0)               # can never co-occur
  expect_equal(m[1, 2], 4 * 3 / (10 * 9))
})

test_that("systematic probabilities follow `order_by`, not frame order", {
  # draw_design() sorts first. Reading residue classes off the unsorted frame
  # reported 0 for pairs that always co-occur and 1/k for pairs that never can.
  pop <- data.frame(id = 1:6, v = c(2, 5, 1, 6, 3, 4))
  des <- design_systematic(interval = 3, order_by = "v")
  m <- joint_prob(pop, des)
  sim <- joint_prob(pop, des, simulate = TRUE, R = 6000, seed = 1)
  expect_lt(max(abs(m - sim)), 0.04)
  expect_equal(m[1, 4], 0)

  # And a dropped row cannot be selected, so its probability is 0.
  na_pop <- data.frame(v = c(1:12, rep(NA, 8)))
  na_des <- design_systematic(interval = 3, order_by = "v", na_rm = TRUE)
  p <- inclusion_prob(na_pop, na_des)
  expect_equal(p, c(rep(1 / 3, 12), rep(0, 8)))
  expect_equal(sum(p), 4)
})

test_that("weights the draw would reject are rejected here too", {
  # pps_pi()'s rescaling loop runs off the rails on a non-positive weight and
  # returns a vector of 0s and 1s summing to the wrong total.
  bad <- data.frame(w = c(1, 2, 3, -4, 5))
  for (m in c("systematic", "poisson")) {
    expect_error(inclusion_prob(bad, design_weighted("w", n = 2, method = m)),
                 "must be positive")
  }
  expect_error(inclusion_prob(data.frame(w = c(1, 2, Inf)),
                              design_weighted("w", n = 2,
                                              method = "systematic")),
               "non-finite")
})

test_that("na_rm means the same thing on the probability path as on the draw", {
  cases <- list(
    list(d = data.frame(g = c(rep("a", 4), rep("b", 4), NA, NA)),
         des = function(na) design_stratified("g", n = 4, na_rm = na)),
    list(d = data.frame(g = c(rep("a", 4), rep("b", 4), NA, NA)),
         des = function(na) design_cluster("g", n_clusters = 1, na_rm = na)),
    list(d = data.frame(w = c(1:8, NA, NA)),
         des = function(na) design_weighted("w", n = 3, method = "poisson",
                                            na_rm = na)),
    list(d = data.frame(v = c(1:8, NA, NA)),
         des = function(na) design_systematic(interval = 3, order_by = "v",
                                              na_rm = na))
  )
  for (case in cases) {
    expect_error(draw(case$d, case$des(FALSE), seed = 1), "na_rm = TRUE")
    expect_error(inclusion_prob(case$d, case$des(FALSE)), "na_rm = TRUE")
    expect_equal(length(inclusion_prob(case$d, case$des(TRUE))), nrow(case$d))
  }
})

test_that("Neyman allocation reports probabilities for a design it can draw", {
  # stratum_spread() checked the whole frame while the draw checked only the
  # rows it kept, so a drawable design had no weights.
  d <- data.frame(g = c(rep("a", 5), rep("b", 5), NA, NA),
                  v = c(seq(1, 5), seq(2, 10, length.out = 5), NA, NA))
  des <- design_stratified("g", n = 4, allocation = "neyman",
                           allocation_by = "v", na_rm = TRUE)
  expect_equal(nrow(draw(d, des, seed = 1)), 4)
  p <- inclusion_prob(d, des)
  expect_equal(length(p), nrow(d))
  expect_equal(p[11:12], c(0, 0))
  expect_equal(nrow(draw(d, des, seed = 1, weights = TRUE)), 4)
})
