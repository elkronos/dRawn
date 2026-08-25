test_that("draws the requested number of rows", {
  df <- make_df()
  expect_equal(nrow(draw(df, design_simple(n = 10), seed = 123)), 10)
})

test_that("refuses to over-draw without replacement, but allows it with", {
  df <- make_df()
  expect_error(draw(df, design_simple(n = 101), seed = 1),
               "cannot exceed the number of rows")
  expect_equal(nrow(draw(df, design_simple(n = 150, replace = TRUE), seed = 1)),
               150)
})

test_that("n is validated at construction time", {
  expect_error(design_simple(n = -1), "must be non-negative")
  expect_error(design_simple(n = 2.7), "single whole number")
  expect_error(design_simple(n = NA), "single whole number")
  expect_error(design_simple(n = c(5, 10)), "single whole number")
  expect_equal(design_simple(n = 0)$n, 0L)
})

test_that("replace is validated", {
  expect_error(design_simple(n = 5, replace = "yes"), "must be TRUE or FALSE")
  expect_error(design_simple(n = 5, replace = NA), "must be TRUE or FALSE")
})

test_that("n = 0 draws nothing but keeps the schema", {
  df <- make_df()
  res <- draw(df, design_simple(n = 0), seed = 1)
  expect_equal(nrow(res), 0L)
  expect_same_schema(res, df)
})

test_that("a tibble comes back a tibble", {
  skip_if_not_installed("tibble")
  tb <- tibble::as_tibble(make_df())
  expect_s3_class(draw(tb, design_simple(n = 3), seed = 1), "tbl_df")
})
