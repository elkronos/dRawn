# A data frame column can itself be a list or a matrix. Neither can serve as a
# grouping or ordering key, and before this check two designs failed silently
# rather than loudly.

odd_df <- function() {
  d <- data.frame(id = 1:20)
  d$lst <- as.list(rep(c("a", "b"), 10))
  d$mat <- matrix(1:40, nrow = 20)
  d
}

test_that("a matrix column is refused rather than silently mis-grouped", {
  d <- odd_df()
  # unique() on a matrix column works element-wise, which silently yields one
  # row from a twenty-row frame rather than an error.
  expect_error(draw(d, design_cluster("mat", n_clusters = 1), seed = 1),
               "which is a 2-dimensional matrix column")
  expect_error(draw(d, design_stratified("mat", n = 4), seed = 1),
               "which is a 2-dimensional matrix column")
  expect_error(draw(d, design_multistage("mat", n_clusters = 1, n = 2), seed = 1),
               "which is a 2-dimensional matrix column")
})

test_that("a list column is refused rather than failing inside base R", {
  d <- odd_df()
  # order() on a list column raises "unimplemented type 'list'" from base R.
  expect_error(draw(d, design_stratified("lst", n = 4), seed = 1),
               "which is a list column")
  expect_error(draw(d, design_systematic(interval = 2, order_by = "lst"), seed = 1),
               "which is a list column")
  expect_error(draw(d, design_cluster("lst", n_clusters = 1), seed = 1),
               "which is a list column")
})

test_that("the error names the offending column and the argument", {
  d <- odd_df()
  expect_error(draw(d, design_stratified(c("id", "mat"), n = 4), seed = 1),
               "`strata` names `mat`")
  expect_error(draw(d, design_systematic(interval = 2, order_by = "lst"), seed = 1),
               "`order_by` names `lst`")
})

test_that("inclusion_prob refuses them too, not just draw()", {
  d <- odd_df()
  expect_error(inclusion_prob(d, design_cluster("mat", n_clusters = 1)),
               "2-dimensional matrix column")
  expect_error(inclusion_prob(d, design_stratified("lst", n = 4)),
               "list column")
})

test_that("factors and ordinary vectors are still accepted", {
  d <- data.frame(
    id = 1:20,
    f  = factor(rep(c("a", "b"), 10)),
    ch = rep(c("x", "y"), 10),
    nu = as.numeric(rep(1:2, 10)),
    dt = rep(as.Date("2020-01-01") + 0:1, 10)
  )
  for (col in c("f", "ch", "nu", "dt")) {
    expect_equal(nrow(draw(d, design_stratified(col, n = 4), seed = 1)), 4,
                 info = col)
    expect_equal(nrow(draw(d, design_cluster(col, n_clusters = 1), seed = 1)), 10,
                 info = col)
  }
  expect_equal(nrow(draw(d, design_systematic(interval = 2, order_by = "nu"),
                         seed = 1)), 10)
})
