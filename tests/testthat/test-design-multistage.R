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
