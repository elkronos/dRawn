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

test_that("a text order_by is walked in byte order in every locale", {
  d <- data.frame(id = 1:12,
                  key = c("a", "B", "c", "D", "e", "F", "g", "H", "i", "J",
                          "k", "L"),
                  y = 1:12)
  des <- design_systematic(interval = 3, order_by = "key")
  walk <- order(d$key, method = "radix")
  pos <- integer(12)
  pos[walk] <- 1:12
  expect_identical(drawn:::systematic_positions(des, d), pos)
  ids <- draw(d, des, seed = 1)$id
  expect_length(ids, 4L)
  expect_length(unique((pos[ids] - 1L) %% 3L), 1L)
  s <- draw(d, des, seed = 1, weights = TRUE)
  want <- ht_total(s, "y")$se
  pij <- joint_prob(d, des)
  with_other_collation({
    expect_identical(drawn:::systematic_positions(des, d), pos)
    expect_identical(draw(d, des, seed = 1)$id, ids)
    expect_equal(ht_total(s, "y")$se, want)
    expect_equal(joint_prob(d, des), pij)
  })
})

# One sample for each start: every start has probability 1 / interval, so
# averages over these are exact expectations.
every_start <- function(d, des) {
  out <- list()
  i <- 0L
  while (length(out) < des$interval && i < 500L) {
    i <- i + 1L
    s <- suppressWarnings(draw(d, des, seed = i, weights = TRUE))
    out[[paste(s$.row, collapse = ",")]] <- s
  }
  unname(out)
}

test_that("a random sample size adds its own variance", {
  # 244 rows at interval 8 give samples of 30 or 31 rows. Successive
  # differences cannot see that, so this near-constant `y` got 0.016 of its
  # true variance, and a constant one a standard error of 0.
  d <- data.frame(y = 1000 + (seq_len(244) * 37) %% 41 - 20, one = 1)
  des <- design_systematic(interval = 8)
  samples <- every_start(d, des)
  expect_length(samples, 8L)
  est <- vapply(samples, function(s) {
    r <- ht_total(s, "y")
    c(r$total, r$variance)
  }, numeric(2))
  ratio <- mean(est[2, ]) / mean((est[1, ] - sum(d$y))^2)
  expect_gt(ratio, 0.8)
  expect_lt(ratio, 1.25)

  # For a constant `y` the term is exact: the total is 8 times the sample
  # size, which is 30 or 31 with equal chance.
  r <- ht_total(samples[[1]], "one")
  expect_equal(r$se, 8 * 0.5)
  expect_match(r$note, "random, 30 or 31 rows")

  # The mean over the known N carries the same term; the Hajek mean's
  # residuals average 0 and do not.
  m <- ht_mean(samples[[1]], "one", estimator = "ht")
  expect_equal(m$se, 8 * 0.5 / 244)
  expect_equal(ht_mean(samples[[1]], "y")$variance,
               ht_mean(samples[[1]], "y", variance = "analytic")$variance)
  h <- vapply(samples, function(s) {
    u <- (s$y - mean(s$y)) / nrow(s)
    k <- 8
    n <- nrow(s)
    c(ht_mean(s, "y")$variance,
      (1 - 1 / k) * n / (2 * (n - 1)) * sum(diff(u)^2))
  }, numeric(2))
  expect_equal(h[1, ], h[2, ])
})

test_that("an interval that divides the frame adds nothing", {
  d <- data.frame(y = 1000 + (seq_len(240) * 37) %% 41 - 20)
  s <- draw(d, design_systematic(interval = 8), seed = 1, weights = TRUE)
  r <- ht_total(s, "y")
  n <- nrow(s)
  expect_equal(r$variance,
               (1 - 1 / 8) * n / (2 * (n - 1)) * sum(diff(8 * s$y)^2))
  expect_false(grepl("random", r$note))
})

test_that("the rows an order_by drops do not count towards the remainder", {
  # 246 rows, 2 with no sort key: 244 walked, 30 or 31 drawn.
  d <- data.frame(key = c(NA, NA, seq_len(244)), one = 1)
  s <- draw(d, design_systematic(interval = 8, order_by = "key", na_rm = TRUE),
            seed = 1, weights = TRUE)
  r <- ht_total(s, "one")
  expect_equal(r$se, 4)
  expect_match(r$note, "The 244 rows walked")
})

test_that("an order_by the frame lacks is refused on every path", {
  # The probability paths gave every row 0 and every weight NA for a
  # misspelt or missing column, where draw() refused it.
  d <- data.frame(id = 1:24, v = (24:1) / 2, y = 1:24)
  des <- design_systematic(interval = 4, order_by = "V")
  want <- "`data` is missing column: `V`."
  expect_error(draw(d, des, seed = 1), want, fixed = TRUE)
  expect_error(inclusion_prob(d, des), want, fixed = TRUE)
  expect_error(sampling_weight(d, des), want, fixed = TRUE)
  expect_error(joint_prob(d, des, rows = 1:3), want, fixed = TRUE)
  expect_error(
    drawn:::ht_variance.drawn_design_systematic(des, d, 1:6, 1:6, rep(0.25, 6)),
    want, fixed = TRUE)
})

test_that("deleting order_by from the frame after the draw leaves the sample", {
  # The stored copy of the frame shared its names with the caller's, so
  # deleting the column in place took it out of the copy too, and ht_total()
  # refused the sample.
  skip_if_not_installed("data.table")
  d <- data.frame(id = 1:80, w = (80:1) / 4, y = sin(1:80) * 10 + 1:80)
  s <- draw(d, design_systematic(interval = 8, order_by = "w"), seed = 1,
            weights = TRUE)
  want <- ht_total(s, "y")
  data.table::set(d, j = "w", value = NULL)
  expect_identical(names(attr(s, "drawn_population")), c("id", "w", "y"))
  expect_equal(ht_total(s, "y"), want)
})
