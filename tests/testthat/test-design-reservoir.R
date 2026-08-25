test_that("a data frame in gives a data frame out", {
  df <- make_df()
  res <- draw(df, design_reservoir(n = 10), seed = 123)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 10)
  expect_same_schema(res, df)
})

test_that("a list stream gives a list back", {
  res <- draw(as.list(1:100), design_reservoir(n = 5), seed = 123)
  expect_type(res, "list")
  expect_equal(length(res), 5)
  expect_true(all(unlist(res) %in% 1:100))
})

test_that("a generator function is accepted", {
  i <- 0L
  gen <- function() {
    i <<- i + 1L
    if (i > 50L) NULL else i
  }
  res <- draw(gen, design_reservoir(n = 5), seed = 1)
  expect_equal(length(res), 5)
  expect_true(all(unlist(res) %in% 1:50))
})

test_that("a short stream is trimmed, not padded with NULLs", {
  res <- draw(as.list(1:5), design_reservoir(n = 10), seed = 1)
  expect_equal(length(res), 5)
  expect_false(any(vapply(res, is.null, logical(1))))

  df_res <- draw(data.frame(x = 1:5), design_reservoir(n = 10), seed = 1)
  expect_equal(nrow(df_res), 5)
})

test_that("max_items warns when it truncates", {
  expect_warning(res <- draw(as.list(1:100), design_reservoir(n = 5, max_items = 20),
                             seed = 1),
                 "may have held more")
  expect_true(all(unlist(res) <= 20))

  expect_warning(draw(data.frame(x = 1:100), design_reservoir(n = 5, max_items = 20),
                      seed = 1),
                 "below the 100 available rows")
})

test_that("max_items is validated", {
  expect_error(design_reservoir(n = 5, max_items = 0), "Inf or a single whole")
  expect_error(design_reservoir(n = 5, max_items = 2.5), "Inf or a single whole")
})

test_that("an unsupported stream type is reported clearly", {
  expect_error(draw(1:10, design_reservoir(n = 3), seed = 1),
               "must be a data frame, list, connection")
})

test_that("the sample is uniform over the stream", {
  set.seed(4)
  counts <- integer(20)
  for (i in 1:4000) {
    drawn <- unlist(draw(as.list(1:20), design_reservoir(n = 5)))
    counts[drawn] <- counts[drawn] + 1L
  }
  expect_gt(suppressWarnings(stats::chisq.test(counts)$p.value), 0.001)
})
