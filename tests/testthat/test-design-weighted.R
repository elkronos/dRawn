weighted_df <- function() {
  set.seed(99)
  data.frame(id = 1:100, weight = stats::runif(100, 0.1, 10))
}

test_that("draws the requested number of rows", {
  df <- weighted_df()
  expect_equal(nrow(draw(df, design_weighted("weight", n = 10), seed = 123)), 10)
})

test_that("the caller's weights come back untouched", {
  d <- data.frame(id = 1:10, w = 1:10)
  out <- draw(d, design_weighted("w", n = 5), seed = 1)
  expect_equal(out$w, d$w[match(out$id, d$id)])
  expect_type(out$w, "integer")
})

test_that("heavier rows are drawn more often", {
  d <- data.frame(id = 1:4, w = c(1, 1, 1, 100))
  draws <- vapply(1:400, function(i) draw(d, design_weighted("w", n = 1),
                                          seed = i)$id, numeric(1))
  expect_gt(mean(draws == 4), 0.8)
})

test_that("even the lightest row remains selectable", {
  d <- data.frame(id = 1:5, w = c(1, 2, 3, 4, 5))
  drawn <- vapply(1:600, function(i) draw(d, design_weighted("w", n = 1),
                                          seed = i)$id, numeric(1))
  expect_true(1 %in% drawn)
})

test_that("non-positive, non-finite and non-numeric weights are rejected", {
  df <- weighted_df()

  bad <- df
  bad$weight <- -abs(bad$weight)
  expect_error(draw(bad, design_weighted("weight", n = 10), seed = 1),
               "must be positive")

  inf_df <- df
  inf_df$weight[1] <- Inf
  expect_error(draw(inf_df, design_weighted("weight", n = 10), seed = 1),
               "non-finite")

  chr <- df
  chr$weight <- as.character(chr$weight)
  expect_error(draw(chr, design_weighted("weight", n = 10), seed = 1),
               "must be numeric")
})

test_that("NA weights follow the same na_rm contract as other keys", {
  df <- weighted_df()
  df$weight[1] <- NA
  expect_error(draw(df, design_weighted("weight", n = 10), seed = 1),
               "missing weight")
  expect_equal(
    nrow(draw(df, design_weighted("weight", n = 10, na_rm = TRUE), seed = 1)),
    10
  )
})

test_that("over-drawing without replacement is refused", {
  d <- data.frame(id = 1:5, w = 1:5)
  expect_error(draw(d, design_weighted("w", n = 10), seed = 1),
               "cannot exceed the number of rows")
})
