test_that("simple resampling returns n_replicates replicates of n rows", {
  df <- make_df()
  res <- draw(df, design_bootstrap(n_replicates = 5, n = 10), seed = 123)
  expect_equal(nrow(res), 50)
  expect_equal(as.vector(table(res$.replicate)), rep(10L, 5))
  expect_identical(names(res), c(".replicate", names(df)))
})

test_that("n defaults to nrow(data)", {
  df <- make_df()
  res <- draw(df, design_bootstrap(n_replicates = 2), seed = 1)
  expect_equal(as.vector(table(res$.replicate)), rep(nrow(df), 2))
})

test_that("block replicates are stitched from blocks across the whole frame", {
  d <- data.frame(id = 1:1000)
  res <- draw(d, design_bootstrap(n_replicates = 20, n = 50, method = "block",
                                  block_length = 5), seed = 123)
  reps <- split(res, res$.replicate)

  # One contiguous window of 50 could span at most 50 positions.
  spans <- vapply(reps, function(z) diff(range(z$id)), numeric(1))
  expect_true(all(spans > 200))

  # Several separate runs per replicate, not one.
  runs <- vapply(reps, function(z) {
    u <- sort(unique(z$id))
    sum(diff(u) > 1L) + 1L
  }, numeric(1))
  expect_true(all(runs >= 5))

  expect_equal(as.vector(table(res$.replicate)), rep(50L, 20))
})

test_that("block_length is validated against the data", {
  df <- make_df()
  expect_error(
    draw(df, design_bootstrap(n_replicates = 2, method = "block",
                              block_length = 500), seed = 1),
    "exceeds the number of rows"
  )
})

test_that("an unknown method is rejected at construction", {
  expect_error(design_bootstrap(method = "invalid"))
})

test_that("replicates can be split apart and summarised", {
  df <- data.frame(id = 1:50, value = (1:50) / 10)
  res <- draw(df, design_bootstrap(n_replicates = 8, n = 20), seed = 1)
  means <- vapply(split(res, res$.replicate), function(r) mean(r$value),
                  numeric(1))
  expect_length(means, 8)
  expect_true(all(is.finite(means)))
})
