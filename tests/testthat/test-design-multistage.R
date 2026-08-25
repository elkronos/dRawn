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
