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

test_that("design_systematic() takes `interval` or `n`, and `start` with `interval`", {
  expect_identical(design_systematic(7), design_systematic(interval = 7))
  expect_error(design_systematic(), "needs `interval` or `n`")
  expect_error(design_systematic(5, n = 3), "not both")
  expect_error(design_systematic(n = 3, start = 1),
               "`start` cannot be fixed with `n`")
  expect_error(design_systematic(interval = 2.5),
               "give the sample size as `n`")
  expect_error(design_systematic(n = 2.5), "single whole number")
  expect_error(design_systematic(n = -1), "non-negative")
  des <- design_systematic(n = 7)
  expect_output(print(des), "n +7")
  expect_identical(drawn:::design_label(des), "design_systematic(n = 7)")
})

test_that("n draws exactly n rows at a step that need not be whole", {
  d <- data.frame(id = 1:50, v = 50:1)
  for (seed in 1:20) {
    ids <- draw(d, design_systematic(n = 7), seed = seed)$id
    expect_length(ids, 7L)
    # A step of 50 / 7 rows leaves gaps of 7 or 8 along the walk.
    expect_true(all(diff(ids) %in% 7:8), info = paste("seed", seed))
    # With order_by the gaps of 12 or 13 are counted along the sorted walk.
    ids <- draw(d, design_systematic(n = 4, order_by = "v"), seed = seed)$id
    expect_length(ids, 4L)
    expect_true(all(diff(sort(d$v[ids])) %in% 12:13),
                info = paste("seed", seed))
  }
  expect_equal(nrow(draw(d, design_systematic(n = 0), seed = 1)), 0L)
  expect_identical(draw(d, design_systematic(n = 50), seed = 1)$id, 1:50)
})

test_that("n gives every row n / N and every pair its overlap", {
  d <- data.frame(id = 1:50)
  des <- design_systematic(n = 7)
  expect_equal(inclusion_prob(d, des), rep(7 / 50, 50))
  jp <- joint_prob(d, des)
  expect_equal(diag(jp), rep(7 / 50, 50))
  # The size is fixed, so each row's pairs add up to n times its own pi.
  expect_equal(rowSums(jp), rep(7 * 7 / 50, 50))
  expect_equal(jp[1, 2], 0)
  sim <- joint_prob(d, des, simulate = TRUE, R = 2000, seed = 1)
  expect_lt(max(abs(jp - sim)), 0.04)
  # A step that divides N is the interval design.
  d60 <- data.frame(id = 1:60)
  expect_equal(joint_prob(d60, design_systematic(n = 10)),
               joint_prob(d60, design_systematic(interval = 6)))
})

test_that("averaging over the start gives the joint probabilities and the total", {
  # The rows taken change only where `u` crosses frac(j * n / N), so the
  # expectation over `u` is a sum over the pieces between those points.
  N <- 13L
  n <- 5L
  d <- data.frame(y = c(4, 9, 1, 7, 3, 12, 6, 2, 8, 5, 11, 10, 15))
  cuts <- sort(unique(c(0, 1, (seq_len(N) * n / N) %% 1)))
  mid <- (cuts[-1] + cuts[-length(cuts)]) / 2
  len <- diff(cuts)
  pairs <- matrix(0, N, N)
  total <- 0
  for (i in seq_along(mid)) {
    s <- drawn:::fractional_walk(mid[i], n, N)
    expect_length(unique(s), n)
    pairs[s, s] <- pairs[s, s] + len[i]
    total <- total + len[i] * sum(d$y[s]) * N / n
  }
  expect_equal(joint_prob(d, design_systematic(n = n)), pairs)
  expect_equal(total, sum(d$y))
})

test_that("samples drawn with `interval` are the same for a seed as before", {
  d <- data.frame(id = 1:100, value = (1:100) / 10)
  des <- design_systematic(interval = 7)
  expect_identical(draw(d, des, seed = 1)$id, seq(1L, 100L, by = 7L))
  expect_identical(draw(d, des, seed = 2)$id, seq(5L, 100L, by = 7L))
  expect_identical(draw(d, des, seed = 4)$id, seq(3L, 100L, by = 7L))
  expect_identical(
    draw(d, design_systematic(interval = 6, order_by = "value"), seed = 9)$id,
    seq(3L, 100L, by = 6L))
})

test_that("a fixed-size sample has no term for a random size", {
  d <- data.frame(y = 1000 + (seq_len(244) * 37) %% 41 - 20, one = 1)
  s <- draw(d, design_systematic(n = 30), seed = 1, weights = TRUE)
  expect_equal(nrow(s), 30L)
  # Every sample gives the frame's row count exactly.
  r <- ht_total(s, "one")
  expect_equal(r$total, 244)
  expect_equal(r$se, 0)
  r <- ht_total(s, "y")
  u <- s$y * 244 / 30
  expect_equal(r$variance, (1 - 30 / 244) * 30 / (2 * 29) * sum(diff(u)^2))
  expect_false(grepl("random", r$note))
  expect_match(r$note, "the step of 8.133 rows", fixed = TRUE)
  expect_equal(r$df, 29)
})

test_that("an n above the rows walked is refused on every path", {
  d <- data.frame(v = c(NA, 1:9))
  des <- design_systematic(n = 10, order_by = "v", na_rm = TRUE)
  want <- "`n` (10) cannot exceed the number of rows walked (9)."
  expect_error(draw(d, des, seed = 1), want, fixed = TRUE)
  expect_error(inclusion_prob(d, des), want, fixed = TRUE)
  expect_error(joint_prob(d, des), want, fixed = TRUE)
  # The row with no sort key is never walked.
  ok <- design_systematic(n = 4, order_by = "v", na_rm = TRUE)
  expect_equal(inclusion_prob(d, ok), c(0, rep(4 / 9, 9)))
  jp <- joint_prob(d, ok)
  expect_equal(jp[1, ], rep(0, 10))
  expect_equal(rowSums(jp)[-1], rep(4 * 4 / 9, 9))
})

test_that("a fixed-size systematic sample goes to survey with the rows walked", {
  skip_if_not_installed("survey")
  d <- data.frame(id = 1:50, v = c(NA, 49:1), y = sin(1:50) * 10 + 1:50)
  s <- draw(d, design_systematic(n = 7, order_by = "v", na_rm = TRUE),
            seed = 3, weights = TRUE)
  sv <- as_svydesign(s)
  r <- ht_total(s, "y")
  expect_equal(unname(coef(survey::svytotal(~y, sv))), r$total)
  expect_equal(survey::degf(sv), r$df)
  expect_equal(unname(sv$fpc$popsize[, 1]), rep(49, 7))
})
