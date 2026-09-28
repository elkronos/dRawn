test_that("walks the data at a fixed interval", {
  df <- make_df()
  expect_equal(draw(df, design_systematic(interval = 10, start = 1))$id,
               seq(1, 100, by = 10))
})

test_that("order_by sorts before walking, and rows come back in frame order", {
  df <- make_df()
  res <- draw(df, design_systematic(interval = 10, order_by = "value"), seed = 123)
  # Every selected row sits the same distance from the start of the sorted walk
  pos <- rank(df$value)[res$id]
  expect_length(unique((pos - 1) %% 10), 1L)
  expect_length(res$id, 10L)
  # ...and the contract in ?designs: frame order, not walk order
  expect_identical(res$id, sort(res$id))
})

test_that("a random start lands inside 1:interval", {
  df <- make_df()
  starts <- vapply(1:50, function(i) {
    draw(df, design_systematic(interval = 7), seed = i)$id[1]
  }, numeric(1))
  expect_true(all(starts >= 1 & starts <= 7))
})

test_that("interval and start are validated at construction", {
  expect_error(design_systematic(interval = 0), "at least 1")
  expect_error(design_systematic(interval = 2.5), "single whole number")
  expect_error(design_systematic(interval = 10, start = 11),
               "must lie in 1:interval")
})

test_that("a start past the end warns and returns an empty frame", {
  df <- make_df()
  expect_warning(res <- draw(df, design_systematic(interval = 200, start = 150)),
                 "past the last row")
  expect_equal(nrow(res), 0L)
  expect_same_schema(res, df)
})

test_that("NA sort keys are reported rather than silently sorted last", {
  df <- data.frame(v = c(5, NA, 4, 2, 3), z = letters[1:5])
  expect_error(draw(df, design_systematic(interval = 2, order_by = "v", start = 1)),
               "missing `v`")
  res <- draw(df, design_systematic(interval = 2, order_by = "v", start = 1,
                                    na_rm = TRUE))
  expect_false(anyNA(res$v))
})
