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
