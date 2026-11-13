test_that("selects whole clusters", {
  df <- make_df()
  res <- draw(df, design_cluster("site", n_clusters = 2), seed = 123)
  expect_equal(length(unique(res$site)), 2)
  expect_equal(nrow(res), 50)
})

test_that("balanced = TRUE equalises cluster sizes", {
  df <- make_df()
  uneven <- df[c(1:25, 26:45, 51:65, 76:100), ]
  res <- draw(uneven, design_cluster("site", n_clusters = 4, balanced = TRUE),
              seed = 123)
  expect_equal(length(unique(as.vector(table(res$site)))), 1)
})

test_that("asking for more clusters than exist is an error", {
  df <- make_df()
  expect_error(draw(df, design_cluster("site", n_clusters = 10), seed = 1),
               "exceeds the 4 cluster")
})

test_that("NA cluster labels are reported rather than treated as a cluster", {
  df <- data.frame(x = 1:30, cl = c(rep("a", 10), rep("b", 10), rep(NA, 10)))
  expect_error(draw(df, design_cluster("cl", n_clusters = 2), seed = 1),
               "missing cluster label")

  res <- draw(df, design_cluster("cl", n_clusters = 2, na_rm = TRUE), seed = 1)
  expect_false(anyNA(res$cl))
  expect_equal(nrow(res), 20)
})

test_that("clusters must name a single column", {
  expect_error(design_cluster(c("a", "b"), n_clusters = 2), "a single column")
})

test_that("the clusters chosen for a seed are sample(unique(labels))'s", {
  d <- data.frame(id = 1:60, cl = rep(c(5, 2, 9, 1, 7, 3), each = 10))
  for (k in c(1, 3)) {
    got <- unique(draw(d, design_cluster("cl", n_clusters = k), seed = 7)$cl)
    want <- drawn:::with_seed(7, sample(unique(d$cl), k))
    expect_setequal(got, want)
  }
})

test_that("two cluster ids that print alike stay two clusters", {
  # 0.3 and 0.1 + 0.2 both print as 0.3. Compared as text they were one
  # cluster: a two-cluster sample got no variance, and two rows of different
  # clusters got the same-cluster joint probability.
  d <- data.frame(id = 1:20, cl = rep(c(0.3, 0.1 + 0.2, 0.5, 0.7), each = 5),
                  y = c(1:5, 101:105, 21:25, 61:65))
  des <- design_cluster("cl", n_clusters = 3)
  expect_equal(inclusion_prob(d, des), rep(0.75, 20))
  m <- joint_prob(d, des, rows = c(1, 2, 6))
  expect_equal(m[1, 2], 0.75)
  expect_equal(m[1, 3], 3 * 2 / (4 * 3))

  s <- draw(d, des, seed = 2, weights = TRUE)
  expect_equal(sum(s$cl == 0.3), 5)
  expect_equal(sum(s$cl == 0.1 + 0.2), 5)
  u <- tapply(s$y / s$.prob, match(s$cl, unique(s$cl)), sum)
  r <- ht_total(s, "y")
  expect_equal(r$method, "analytic")
  expect_equal(r$variance, (1 - 3 / 4) * 3 * stats::var(u))
  expect_equal(r$df, 2)
  expect_equal(ht_total(s, "y", variance = "jackknife")$variance, r$variance)
  expect_equal(nrow(sample_summary(s)$by_group), 3)

  skip_if_not_installed("survey")
  des_sv <- as_svydesign(s)
  expect_equal(survey::degf(des_sv), 2)
  expect_equal(as.numeric(survey::SE(survey::svytotal(~y, des_sv))), r$se)
})

test_that("a factor's unused levels are not clusters", {
  d <- data.frame(id = 1:40, y = 1:40,
                  cl = factor(rep(c("a", "b", "c", "d"), each = 10),
                              levels = c("a", "b", "c", "d", "zz")))
  des <- design_cluster("cl", n_clusters = 2)
  expect_equal(inclusion_prob(d, des), rep(0.5, 40))
  s <- draw(d, des, seed = 1, weights = TRUE)
  expect_equal(nrow(s), 20)
  r <- ht_total(s, "y")
  expect_equal(r$method, "analytic")
  expect_equal(r$df, 1)
  expect_equal(nrow(sample_summary(s)$by_group), 2)
})

test_that("a missing cluster label is called that everywhere", {
  d <- data.frame(id = 1:30, cl = c(NA, rep(c("a", "b", "c"), each = 10))[1:30])
  des <- design_cluster("cl", n_clusters = 2)
  expect_error(draw(d, des, seed = 1), "1 row\\(s\\) have a missing cluster label")
  expect_error(inclusion_prob(d, des),
               "1 row\\(s\\) have a missing cluster label")
  expect_error(inclusion_prob(d, design_multistage("cl", n_clusters = 2, n = 4)),
               "1 row\\(s\\) have a missing cluster label")
})

test_that("selecting many clusters takes one pass over the frame", {
  skip_on_cran()
  # Each chosen cluster scanned the whole frame: 16 seconds for 1,000
  # clusters of a million rows.
  N <- 1e6
  d <- data.frame(cl = sprintf("c%07d", (seq_len(N) - 1) %/% 50 + 1))
  took <- system.time({
    s <- draw(d, design_cluster("cl", n_clusters = 1000), seed = 1)
    m <- draw(d, design_multistage("cl", n_clusters = 1000, n = 5000), seed = 1)
  })[["elapsed"]]
  expect_equal(nrow(s), 50000)
  expect_equal(length(unique(s$cl)), 1000)
  expect_equal(nrow(m), 5000)
  expect_lt(took, 5)
})
