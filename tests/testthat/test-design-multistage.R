test_that("n is the total across selected clusters, split evenly by default", {
  df <- make_df()
  res <- draw(df, design_multistage("site", n_clusters = 2, n = 6), seed = 123)
  expect_equal(nrow(res), 6)
  expect_equal(length(unique(res$site)), 2)
  expect_equal(as.vector(table(res$site)), c(3L, 3L))
})

test_that("proportional allocation also totals n exactly", {
  df <- make_df()
  res <- draw(df, design_multistage("site", n_clusters = 4, n = 10,
                                    allocation = "proportional"), seed = 123)
  expect_equal(nrow(res), 10)
})

test_that("over-drawing the selected clusters is refused, with the real total", {
  small <- data.frame(x = 1:3, site = c("a", "b", "c"))
  expect_error(
    draw(small, design_multistage("site", n_clusters = 3, n = 6), seed = 123),
    "exceeds the 3 row\\(s\\) available across the 3 selected cluster"
  )
  # The constraint is the SELECTED clusters, not the whole frame.
  df <- make_df()
  expect_error(
    draw(df, design_multistage("site", n_clusters = 1, n = 40), seed = 1),
    "exceeds the 25 row\\(s\\) available across the 1 selected cluster"
  )
  expect_equal(
    nrow(draw(df, design_multistage("site", n_clusters = 1, n = 40,
                                    replace = TRUE), seed = 1)),
    40
  )
})

test_that("NA cluster labels never fabricate rows", {
  df <- data.frame(x = 1:30, cl = c(rep("a", 10), rep("b", 10), rep(NA, 10)))
  expect_error(draw(df, design_multistage("cl", n_clusters = 1, n = 5), seed = 1),
               "missing cluster label")

  res <- draw(df, design_multistage("cl", n_clusters = 1, n = 5, na_rm = TRUE),
              seed = 1)
  expect_false(anyNA(res$x))
  expect_equal(nrow(res), 5)
})

test_that("a single numeric cluster id is not read as a range", {
  d <- data.frame(x = 1:10, cl = rep(5, 10))
  expect_equal(nrow(draw(d, design_multistage("cl", n_clusters = 1, n = 2),
                         seed = 1)), 2)
})

test_that("n below n_clusters, or below what min_per_cluster needs, is refused", {
  expect_error(design_multistage("cl", n_clusters = 3, n = 0),
               "`n` \\(0\\) is smaller than `n_clusters` \\(3\\)")
  expect_error(design_multistage("cl", n_clusters = 5, n = 3),
               "Raise `n` to at least 5")
  expect_error(design_multistage("cl", n_clusters = 4, n = 6,
                                 min_per_cluster = 2),
               "needs at least 8 rows, but `n` is 6")
  expect_s3_class(design_multistage("cl", n_clusters = 4, n = 8,
                                    min_per_cluster = 2),
                  "drawn_design_multistage")
})

test_that("a selected cluster the allocation leaves empty is named", {
  d <- data.frame(id = 1:31, cl = c("tiny", rep("big", 30)))
  des <- design_multistage("cl", n_clusters = 2, n = 10,
                           allocation = "proportional")
  expect_warning(s <- draw(d, des, seed = 1),
                 "selected cluster `tiny` none of the 10 rows",
                 class = "drawn_empty_cluster")
  expect_equal(nrow(s), 10)
  # A cluster smaller than the minimum is taken whole, which is no shortfall.
  expect_no_warning(
    s2 <- draw(d, design_multistage("cl", n_clusters = 2, n = 10,
                                    allocation = "proportional",
                                    min_per_cluster = 2), seed = 1)
  )
  expect_equal(as.vector(table(s2$cl)), c(9L, 1L))
  # Simulation draws many times on purpose and stays quiet.
  expect_no_warning(inclusion_prob(d, des, simulate = TRUE, R = 20, seed = 1))
})

test_that("with replacement, a row's probability is that of being drawn at all", {
  # m / N_h is the expected number of draws of a row, not its probability.
  # Cluster `a` holds 4 rows and takes 5 draws, which replacement allows.
  sizes <- c(4, 8, 12, 16)
  d <- data.frame(id = 1:40, cl = rep(c("a", "b", "c", "d"), times = sizes))
  des <- design_multistage("cl", n_clusters = 2, n = 10, replace = TRUE)
  p <- inclusion_prob(d, des)
  expect_equal(p, 0.5 * (1 - (1 - 1 / rep(sizes, sizes))^5))
  expect_equal(sampling_weight(d, des), 1 / p)
  sim <- inclusion_prob(d, des, simulate = TRUE, R = 2000, seed = 1)
  expect_lt(max(abs(tapply(sim, d$cl, mean) - tapply(p, d$cl, mean))), 0.03)
  expect_error(joint_prob(d, des), "no closed-form joint")
})

test_that("a row with no cluster co-occurs with nothing", {
  d <- data.frame(id = 1:30, cl = c(rep("a", 10), rep("b", 10), NA,
                                    rep("c", 9)))
  des <- design_multistage("cl", n_clusters = 2, n = 4, na_rm = TRUE)
  m <- joint_prob(d, des, rows = c(1, 21, 22))
  expect_equal(m[2, ], c(0, 0, 0))
  expect_equal(m[, 2], c(0, 0, 0))
  expect_false(anyNA(joint_prob(d, des)))
  sim <- joint_prob(d, des, rows = c(1, 21, 22), simulate = TRUE, R = 2000,
                    seed = 1)
  expect_lt(max(abs(sim - m)), 0.03)
})

test_that("a factor's unused levels are not clusters", {
  d <- data.frame(id = 1:40, y = 1:40,
                  cl = factor(rep(c("a", "b", "c", "d"), each = 10),
                              levels = c("a", "b", "c", "d", "zz")))
  des <- design_multistage("cl", n_clusters = 2, n = 4)
  expect_equal(inclusion_prob(d, des), rep(0.1, 40))
  expect_equal(diag(joint_prob(d, des, rows = 1:3)), rep(0.1, 3))
  s <- draw(d, des, seed = 1, weights = TRUE)
  expect_equal(s$.prob, rep(0.1, 4))
  expect_true(is.finite(ht_total(s, "y")$se))
  expect_equal(nrow(sample_summary(s)$by_group), 2)
})

test_that("two cluster ids that print alike stay two clusters", {
  d <- data.frame(id = 1:40, cl = rep(c(0.3, 0.1 + 0.2, 0.5, 0.7), each = 10),
                  y = rep(c(1, 100, 5, 9), each = 10) + 1:40)
  des <- design_multistage("cl", n_clusters = 2, n = 4)
  expect_equal(inclusion_prob(d, des), rep(0.1, 40))
  m <- joint_prob(d, des, rows = c(1, 2, 11))
  expect_equal(m[1, 2], 0.5 * (2 * 1) / (10 * 9))
  expect_equal(m[1, 3], (2 / 10)^2 * (2 * 1) / (4 * 3))

  s <- draw(d, des, seed = 3, weights = TRUE)
  expect_equal(sum(s$cl == 0.3), 2)
  expect_equal(sum(s$cl == 0.1 + 0.2), 2)
  r <- ht_total(s, "y")
  expect_equal(r$method, "analytic")
  expect_equal(r$df, 1)
  expect_equal(nrow(sample_summary(s)$by_group), 2)

  skip_if_not_installed("survey")
  des_sv <- as_svydesign(s)
  expect_equal(survey::degf(des_sv), 1)
  expect_equal(as.numeric(survey::SE(survey::svytotal(~y, des_sv))), r$se)
})

test_that("a sample from one cluster has no variance, and says why", {
  # Sen-Yates-Grundy is negative on every such sample, and was 0 for a
  # constant `y`.
  d <- data.frame(id = 1:40, cl = rep(c("a", "b", "c", "d"), each = 10),
                  y = rep(5, 40))
  s <- draw(d, design_multistage("cl", n_clusters = 1, n = 5), seed = 2,
            weights = TRUE)
  for (v in c("analytic", "auto")) {
    r <- ht_total(s, "y", variance = v)
    expect_true(is.na(r$variance))
    expect_match(r$note, "A variance needs at least two clusters; this sample has 1.",
                 fixed = TRUE)
    expect_false(grepl("unlucky", r$note))
  }
})

test_that("the jackknife fallback says when it inflates the between-cluster part", {
  d <- data.frame(id = 1:100, cl = rep(sprintf("k%02d", 1:10), each = 10),
                  y = rep(1:10 * 10, each = 10) + 1:100 %% 7)
  note_for <- function(a) {
    ht_total(draw(d, design_multistage("cl", n_clusters = a, n = a), seed = 1,
                  weights = TRUE), "y")$note
  }
  half <- note_for(5)
  expect_match(half, "With 5 of the 10 clusters selected", fixed = TRUE)
  expect_match(half, "between clusters at 2 times its true size", fixed = TRUE)
  expect_false(grepl("true size", note_for(2)))
  all <- note_for(10)
  expect_match(all, "Every cluster is selected", fixed = TRUE)
  expect_match(all, "design_stratified()", fixed = TRUE)
})

# ---- clusters in proportion to size -----------------------------------------

pps_clusters <- function() {
  sizes <- c(40, 6, 8, 10, 6, 8, 10, 12)
  d <- data.frame(cl = rep(sprintf("k%d", 1:8), times = sizes))
  d$rows <- rep(sizes, times = sizes)
  d$y <- rep(c(50, 20, 35, 40, 25, 30, 45, 60), times = sizes) +
    seq_len(nrow(d)) %% 5 * 3
  d
}

test_that("size selects clusters in proportion to it, certain ones always", {
  d <- pps_clusters()
  des <- design_multistage("cl", n_clusters = 3, n = 6, size = "rows")
  p <- inclusion_prob(d, des)
  # k1 holds 40 of 100 rows: 3 * 0.4 passes 1, so it is certain, and the
  # other two selections are spread over the 60 rows left.
  first <- c(1, 2 * c(6, 8, 10, 6, 8, 10, 12) / 60)
  expect_equal(p, rep(first * 2 / c(40, 6, 8, 10, 6, 8, 10, 12),
                      times = c(40, 6, 8, 10, 6, 8, 10, 12)))
  expect_equal(sum(p), 6)
  # Below certainty the design is self-weighting: (2 * N_i / 60) * (2 / N_i)
  # for every row.
  expect_equal(unique(round(p[d$cl != "k1"], 12)), round(4 / 60, 12))
  for (i in 1:10) {
    s <- draw(d, des, seed = i)
    expect_true("k1" %in% s$cl)
    expect_identical(as.vector(table(s$cl)), rep(2L, 3))
  }
})

test_that("first-stage selection matches its probabilities in simulation", {
  skip_on_cran()
  d <- pps_clusters()
  des <- design_multistage("cl", n_clusters = 3, n = 6, size = "rows")
  sim <- inclusion_prob(d, des, simulate = TRUE, R = 2000, seed = 1)
  claim <- inclusion_prob(d, des)
  # Compared cluster by cluster, where the counts are large enough to test.
  expect_equal(as.vector(tapply(sim, d$cl, sum)),
               as.vector(tapply(claim, d$cl, sum)), tolerance = 0.06)
})

test_that("with equal measures the variance is the two-stage estimator", {
  d <- pps_clusters()
  d$one <- 1
  s <- draw(d, design_multistage("cl", n_clusters = 4, n = 12), seed = 3,
            weights = TRUE)
  des <- design_multistage("cl", n_clusters = 4, n = 12, size = "one")
  expect_equal(inclusion_prob(d, des),
               inclusion_prob(d, design_multistage("cl", 4, 12)))
  expect_equal(pps_stage_variance(des, d, s$.row, s$y, s$.prob)$variance,
               ht_total(s, "y")$variance)
})

test_that("a certain cluster is a stratum of its sampled rows", {
  d <- pps_clusters()
  s <- draw(d, design_multistage("cl", n_clusters = 4, n = 12, size = "rows"),
            seed = 2, weights = TRUE)
  r <- ht_total(s, "y")
  expect_identical(r$method, "deville")
  # Three clusters below certainty less one, and three rows of k1 less one.
  expect_identical(r$df, 4)
  k1 <- s$cl == "k1"
  within <- 40^2 * (1 - 3 / 40) * stats::var(s$y[k1]) / 3
  rest <- !k1
  u <- rowsum(s$y[rest] / s$.prob[rest], s$cl[rest])[, 1]
  p1 <- 3 * c(k2 = 6, k3 = 8, k4 = 10, k5 = 6, k6 = 8, k7 = 10,
              k8 = 12)[names(u)] / 60
  w_i <- vapply(names(u), function(k) {
    i <- s$cl == k
    N <- d$rows[match(k, d$cl)]
    N^2 * (1 - 3 / N) * stats::var(s$y[i]) / 3
  }, numeric(1))
  expect_equal(r$variance, deville_variance(u * p1, p1)$variance + within +
                 sum(w_i / p1))
  jk <- ht_total(s, "y", variance = "jackknife")
  expect_match(jk$note, "Jackknife over 6 primary sampling unit(s) in 2 strata",
               fixed = TRUE)
})

test_that("one row per cluster, or one cluster below certainty, declines", {
  d <- pps_clusters()
  s <- draw(d, design_multistage("cl", n_clusters = 3, n = 3, size = "rows"),
            seed = 1, weights = TRUE)
  r <- ht_total(s, "y")
  expect_true(is.na(r$se))
  expect_match(r$note, "one row per selected cluster", fixed = TRUE)
  expect_match(r$note, "cluster taken with certainty `k1` has a single",
               fixed = TRUE)
  # Measures that need not be row counts: k2 and k3 are certain, and one
  # selection is left for the other six clusters.
  big <- d
  big$rows[big$cl %in% c("k2", "k3")] <- 500
  s2 <- draw(big, design_multistage("cl", n_clusters = 3, n = 6,
                                    size = "rows"), seed = 1, weights = TRUE)
  r2 <- ht_total(s2, "y", variance = "analytic")
  expect_match(r2$note, "this sample has one, besides the 2 taken with",
               fixed = TRUE)
  expect_match(ht_total(s2, "y")$note,
               "The jackknife needs at least two clusters below certainty",
               fixed = TRUE)
})

test_that("a size that is not one per cluster, or not usable, is refused", {
  d <- pps_clusters()
  bad <- d
  bad$rows[2] <- 7
  des <- design_multistage("cl", n_clusters = 3, n = 6, size = "rows")
  expect_error(draw(bad, des, seed = 1),
               "one size for each cluster, but cluster `k1` has 40 and 7",
               fixed = TRUE)
  expect_error(inclusion_prob(bad, des), "one size for each cluster")
  bad$rows <- as.character(d$rows)
  expect_error(inclusion_prob(bad, des), "must be numeric", fixed = TRUE)
  expect_error(design_multistage("cl", n_clusters = 3, n = 7, size = "rows"),
               "Use n = 6 or n = 9.", fixed = TRUE)
  expect_error(design_multistage("cl", n_clusters = 3, n = 6, size = "rows",
                                 allocation = "proportional"),
               "`allocation` is \"proportional\"", fixed = TRUE)
  expect_error(joint_prob(d, des, rows = 1:3), "simulate = TRUE",
               fixed = TRUE)

  gap <- d
  gap$rows[45] <- NA
  expect_error(inclusion_prob(gap, des),
               "1 row(s) have a missing size in `rows`", fixed = TRUE)
  kept <- design_multistage("cl", n_clusters = 3, n = 6, size = "rows",
                            na_rm = TRUE)
  p <- inclusion_prob(gap, kept)
  expect_identical(p[45], 0)
  expect_equal(sum(p), 6)
})

# ---- weights conditional on the clusters selected ----------------------------

test_that("a take that depends on the clusters selected is weighted on them", {
  sizes <- c(6, 6, 9, 4, 6)
  d <- data.frame(id = seq_len(sum(sizes)),
                  cl = rep(c("a", "b", "c", "d", "e"), sizes))
  d$y <- (d$id * 7) %% 13 + rep(c(5, 20, 9, 30, 14), sizes)
  N_i <- rep(sizes, sizes)
  designs <- list(
    # A remainder, dealt in the order the clusters were selected when their
    # sizes tie.
    design_multistage("cl", n_clusters = 2, n = 7),
    # Cluster d is smaller than the take of 5, and is taken whole.
    design_multistage("cl", n_clusters = 3, n = 15),
    design_multistage("cl", n_clusters = 3, n = 10,
                      allocation = "proportional", min_per_cluster = 1)
  )
  for (des in designs) {
    expect_error(inclusion_prob(d, des), "no closed-form")
    a <- des$n_clusters
    for (seed in 1:25) {
      s <- draw(d, des, seed = seed, weights = TRUE)
      expect_identical(s$id, draw(d, des, seed = seed)$id)
      m <- as.vector(table(s$cl)[s$cl])
      expect_equal(s$.prob, a / 5 * m / N_i[s$.row])
      expect_setequal(attr(s, "drawn_selection"),
                      match(unique(s$cl), unique(d$cl)))
    }
    r <- ht_total(s, "y")
    expect_equal(r$method, "analytic")
    expect_equal(r$df, a - 1)
    expect_equal(r$total, sum(s$y / s$.prob))
  }
})

test_that("conditional weights are unbiased, by enumeration", {
  # Every ordered selection of 2 of 4 clusters, then every subset of the
  # rows each selected cluster gives. Ties in the allocation are broken by
  # the order of selection, so both orders are counted.
  sizes <- c(5, 3, 2, 4)
  d <- data.frame(cl = rep(1:4, sizes))
  d$y <- c(4, 9, 1, 7, 12, 30, 22, 26, 3, 5, 16, 18, 11, 20)
  des <- design_multistage("cl", n_clusters = 2, n = 5,
                           allocation = "proportional", min_per_cluster = 1)
  est <- numeric(0)
  for (i in 1:4) for (j in setdiff(1:4, i)) {
    p <- drawn:::conditional_inclusion(des, d, c(i, j))
    expect_true(all(is.na(p[!d$cl %in% c(i, j)])))
    # `.prob` is a / A times the chance of a row given the selection.
    take <- round(p * 4 / 2 * sizes[d$cl])
    total <- 0
    for (k in c(i, j)) {
      rows <- which(d$cl == k)
      subsets <- utils::combn(rows, take[rows[1]], simplify = FALSE)
      total <- total + mean(vapply(subsets, function(r) sum(d$y[r] / p[r]),
                                   numeric(1)))
    }
    expect_equal(take[d$cl == i][1] + take[d$cl == j][1], 5)
    est <- c(est, total)
  }
  expect_equal(mean(est), sum(d$y))
})

test_that("conditionally weighted totals are unbiased over repeated draws", {
  skip_on_cran()
  sizes <- c(4, 9, 12, 5, 20, 7, 15, 6, 10, 8)
  d <- data.frame(cl = rep(seq_along(sizes), sizes))
  d$y <- 3 * d$cl + seq_len(nrow(d)) %% 7 + rep(sizes, sizes)
  des <- design_multistage("cl", n_clusters = 4, n = 15,
                           allocation = "proportional", min_per_cluster = 2)
  totals <- vapply(1:1000, function(i) {
    s <- draw(d, des, seed = i, weights = TRUE)
    sum(s$y * s$.weight)
  }, numeric(1))
  expect_lt(abs(mean(totals) - sum(d$y)), 3 * stats::sd(totals) / sqrt(1000))
})

test_that("proportional weights need a row from every selected cluster", {
  d <- data.frame(id = 1:31, cl = c("tiny", rep("big", 30)), y = 1:31)
  des <- design_multistage("cl", n_clusters = 2, n = 10,
                           allocation = "proportional")
  expect_error(draw(d, des, seed = 1, weights = TRUE),
               "`min_per_cluster` is 0, so proportional allocation can give",
               fixed = TRUE)
  s <- draw(d, design_multistage("cl", n_clusters = 2, n = 10,
                                 allocation = "proportional",
                                 min_per_cluster = 1), seed = 1,
            weights = TRUE)
  expect_equal(as.vector(table(s$cl)[c("big", "tiny")]), c(9L, 1L))
  # Both clusters selected: no variance between them, and the one taken
  # whole adds none within.
  r <- ht_total(s, "y")
  big <- s$y[s$cl == "big"]
  expect_equal(r$variance, 30^2 * (1 - 9 / 30) * stats::var(big) / 9)
})

test_that("one row from a selected cluster of several falls back", {
  sizes <- c(2, 20, 20, 20)
  d <- data.frame(cl = rep(c("s", "x", "y", "z"), sizes))
  d$y <- seq_len(nrow(d)) %% 11
  des <- design_multistage("cl", n_clusters = 2, n = 8,
                           allocation = "proportional", min_per_cluster = 1)
  seed <- which(vapply(1:40, function(i) {
    "s" %in% draw(d, des, seed = i)$cl
  }, logical(1)))[1]
  s <- draw(d, des, seed = seed, weights = TRUE)
  expect_equal(sum(s$cl == "s"), 1L)
  r <- ht_total(s, "y")
  expect_equal(r$method, "jackknife")
  expect_match(r$note, "one row from a selected cluster of several")
})

test_that("a conditional sample keeps the clusters it was drawn from", {
  sizes <- c(6, 6, 9, 4, 6)
  d <- data.frame(id = seq_len(sum(sizes)),
                  cl = rep(c("a", "b", "c", "d", "e"), sizes))
  d$y <- d$id %% 5
  s <- draw(d, design_multistage("cl", n_clusters = 3, n = 13), seed = 2,
            weights = TRUE)
  expect_s3_class(ht_total(s, "y"), "drawn_ht")
  expect_s3_class(ht_mean(s[rev(seq_len(nrow(s))), ], "y"), "drawn_mean")
  expect_equal(sample_summary(s)$unreachable, 0L)

  lost <- s
  attr(lost, "drawn_selection") <- NULL
  expect_error(ht_total(lost, "y"),
               "does not record which clusters were selected")
  bad <- s
  attr(bad, "drawn_selection") <- c(1L, 1L, 2L)
  expect_error(ht_total(bad, "y"), "not 3 distinct clusters of the 5")
  other <- s
  sel <- attr(s, "drawn_selection")
  attr(other, "drawn_selection") <- c(setdiff(1:5, sel)[1], sel[-1])
  expect_error(ht_total(other, "y"), "different inclusion probability")
  edited <- s
  edited$.prob[1] <- edited$.prob[1] / 2
  expect_error(ht_total(edited, "y"), "differ from the inclusion probabilities")
})
