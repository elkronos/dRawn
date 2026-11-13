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
  # Without either it is built, for design_stratified(within =), and refused
  # wherever it would have to draw on its own.
  bare <- design_systematic(order_by = "id")
  d <- data.frame(id = 1:10)
  want <- "`design_systematic()` has no `interval` or `n`"
  expect_error(draw(d, bare, seed = 1), want, fixed = TRUE)
  expect_error(inclusion_prob(d, bare), want, fixed = TRUE)
  expect_error(joint_prob(d, bare), want, fixed = TRUE)
  expect_error(design_systematic(start = 2), "`start` needs `interval`")
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

# ---- several random starts ---------------------------------------------------

# Every set of `m` starts from 1:(m * k), each walked at a step of m * k:
# the rows each set takes, all equally likely.
every_start_set <- function(N, k, m) {
  sets <- utils::combn(m * k, m)
  lapply(seq_len(ncol(sets)), function(j) {
    st <- sets[, j]
    sort(unlist(lapply(st[st <= N], function(s) seq(s, N, by = m * k))))
  })
}

test_that("n_starts is checked when the design is built", {
  des <- design_systematic(interval = 10, n_starts = 3)
  expect_identical(des$n_starts, 3L)
  expect_identical(design_systematic(10)$n_starts, 1L)
  expect_output(print(des), "n_starts +3")
  expect_identical(drawn:::design_label(des),
                   "design_systematic(interval = 10, n_starts = 3)")
  expect_identical(drawn:::design_label(design_systematic(10)),
                   "design_systematic(interval = 10)")
  expect_error(design_systematic(10, n_starts = 0), "at least 1")
  expect_error(design_systematic(10, n_starts = 1.5), "single whole number")
  expect_error(design_systematic(n = 4, n_starts = 2),
               "Give `interval` in place of `n`")
  expect_error(design_systematic(n_starts = 2), "needs `interval`")
  expect_error(design_systematic(3, start = 1, n_starts = 2),
               "`start` cannot be fixed with `n_starts` (2)", fixed = TRUE)
  expect_error(design_systematic(1e9, n_starts = 3), "must be at most")
})

test_that("n_starts draws that many replicates, each at a step of n_starts * interval", {
  d <- data.frame(id = 1:120)
  for (seed in 1:10) {
    ids <- draw(d, design_systematic(interval = 10, n_starts = 4),
                seed = seed)$id
    # 120 rows at a step of 40: three rows from each of four starts.
    expect_length(ids, 12L)
    expect_equal(as.vector(table((ids - 1L) %% 40L)), rep(3L, 4))
  }
  # A seed draws the rows it drew before n_starts existed.
  expect_identical(draw(d, design_systematic(7, n_starts = 1), seed = 2)$id,
                   seq(5L, 120L, by = 7L))
})

test_that("several starts give every pair its exact joint probability", {
  # 23 rows, interval 3 and 2 starts from 1:6: averaging over all 15 sets
  # of starts gives the joint probabilities exactly.
  d <- data.frame(y = c(4, 9, 1, 7, 3, 12, 6, 2, 8, 5, 11, 10, 15, 3, 9, 14,
                        2, 6, 8, 1, 13, 7, 5))
  des <- design_systematic(interval = 3, n_starts = 2)
  sets <- every_start_set(23L, 3L, 2L)
  pairs <- matrix(0, 23, 23)
  for (rows in sets) pairs[rows, rows] <- pairs[rows, rows] + 1 / length(sets)
  jp <- joint_prob(d, des)
  expect_equal(jp, pairs)
  expect_equal(inclusion_prob(d, des), rep(1 / 3, 23))
  # Same residue class modulo 6, then two classes.
  expect_equal(jp[1, 7], 1 / 3)
  expect_equal(jp[1, 2], 2 * 1 / (6 * 5))
  sim <- joint_prob(d, des, simulate = TRUE, R = 2000, seed = 1)
  expect_lt(max(abs(jp - sim)), 0.04)
})

test_that("the variance between replicates is unbiased over every set of starts", {
  check <- function(y, k, m) {
    N <- length(y)
    des <- design_systematic(interval = k, n_starts = m)
    est <- vapply(every_start_set(N, k, m), function(rows) {
      c(k * sum(y[rows]),
        drawn:::replicate_variance(des, seq_len(N), rows, y[rows],
                                   rep(1 / k, length(rows)))$variance)
    }, numeric(2))
    expect_equal(mean(est[1, ]), sum(y))
    expect_equal(mean(est[2, ]), mean((est[1, ] - sum(y))^2))
  }
  # 23 rows at a step of 6, so the replicates differ in size.
  check(c(4, 9, 1, 7, 3, 12, 6, 2, 8, 5, 11, 10, 15, 3, 9, 14, 2, 6, 8, 1,
          13, 7, 5), k = 3L, m = 2L)
  # 10 rows at a step of 12: starts 11 and 12 draw nothing, and count with a
  # total of 0.
  check(c(5, 1, 8, 3, 9, 2, 7, 4, 6, 10), k = 4L, m = 3L)
})

test_that("ht_total() uses the replicates, on n_starts - 1 degrees of freedom", {
  d <- data.frame(y = sin(1:24) * 10 + 1:24, one = 1)
  des <- design_systematic(interval = 4, n_starts = 3)
  s <- draw(d, des, seed = 5, weights = TRUE)
  r <- ht_total(s, "y")
  t_r <- vapply(split(s$y / s$.prob, (s$.row - 1L) %% 12L), sum, 0) * 3
  expect_equal(r$variance, (1 - 1 / 4) * var(t_r) / 3)
  expect_identical(r$method, "analytic")
  expect_equal(r$df, 2)
  # With a fixed size, the Sen-Yates-Grundy sum over the exact joint
  # probabilities gives the same figure.
  pij <- joint_prob(d, des, rows = s$.row)
  expect_equal(r$variance,
               drawn:::syg_variance(s$y, s$.prob, pij)$variance)
  # Every replicate holds two rows here, so the count is exact.
  expect_equal(ht_total(s, "one")$se, 0)
  expect_match(ht_total(s, "y", variance = "jackknife")$note,
               "starts as its primary sampling units")
})

test_that("a Hajek mean over one replicate's rows has no standard error", {
  # Its residuals sum to 0 inside the replicate, so every replicate total
  # is 0 and the variance between them would be 0 too.
  d <- data.frame(y = (1:60 * 7) %% 23)
  s <- draw(d, design_systematic(interval = 6, n_starts = 2), seed = 1,
            weights = TRUE)
  s$dom <- ifelse((s$.row - 1L) %% 12L == (s$.row[1] - 1L) %% 12L, "x", "z")
  r <- ht_mean(s, "y", by = "dom")
  expect_true(all(is.na(r$se)))
  expect_match(attr(r, "notes")[1],
               "sampled rows in one replicate only", fixed = TRUE)
  expect_true(is.finite(ht_mean(s, "y")$se))
})

test_that("several starts go to survey as replicates with the same standard error", {
  skip_if_not_installed("survey")
  d <- data.frame(id = 1:50, v = c(NA, 49:1), y = sin(1:50) * 10 + 1:50,
                  g = rep(c("a", "b"), 25))
  des <- design_systematic(interval = 5, n_starts = 3, order_by = "v",
                           na_rm = TRUE)
  s <- draw(d, des, seed = 3, weights = TRUE)
  sv <- as_svydesign(s)
  r <- ht_total(s, "y")
  st <- survey::svytotal(~y, sv)
  expect_equal(unname(coef(st)), r$total)
  expect_equal(as.numeric(survey::SE(st)), r$se)
  expect_equal(survey::degf(sv), r$df)
  expect_equal(unname(sv$fpc$popsize[, 1]), rep(15, nrow(s)))
  by <- ht_mean(s, "y", by = "g")
  sb <- survey::svyby(~y, ~g, sv, survey::svymean)
  expect_equal(by$se, unname(sb$se))

  # A start past the last row drew nothing, and survey would count the
  # replicates with rows only.
  d2 <- data.frame(y = 1:10)
  des2 <- design_systematic(interval = 4, n_starts = 3)
  s2 <- draw(d2, des2, seed = 2, weights = TRUE)
  expect_length(unique((s2$.row - 1L) %% 12L), 2L)
  expect_error(as_svydesign(s2), "1 of its 3 starts fell past the last")
  expect_true(is.finite(ht_total(s2, "y")$se))

  # An interval of 1 takes every row, whatever the starts.
  s3 <- draw(data.frame(y = c(3, 8)), design_systematic(1, n_starts = 3),
             seed = 1, weights = TRUE)
  r3 <- ht_total(s3, "y")
  expect_equal(c(r3$total, r3$se), c(11, 0))
  expect_match(r3$note, "taken with certainty")
  expect_equal(as.numeric(survey::SE(survey::svytotal(~y, as_svydesign(s3)))),
               0)
})

test_that("a certainty design over several starts estimates, and goes to survey", {
  d <- data.frame(y = c(1000, 2000, 1:48))
  des <- design_certainty("y", 500, design_systematic(interval = 4,
                                                      n_starts = 2))
  s <- draw(d, des, seed = 2, weights = TRUE)
  r <- ht_total(s, "y")
  rest <- s[s$.prob < 1, ]
  t_r <- vapply(split(rest$y / rest$.prob, (rest$.row - 3L) %% 8L), sum, 0) * 2
  expect_equal(r$variance, (1 - 1 / 4) * var(t_r) / 2)
  expect_equal(r$df, 1)
  skip_if_not_installed("survey")
  # The certainty rows are one more unit beside the replicates.
  sv <- as_svydesign(s)
  expect_equal(as.numeric(survey::SE(survey::svytotal(~y, sv))), r$se,
               tolerance = 1e-8)
  expect_equal(survey::degf(sv), 1)
})

test_that("order_by takes one column, and a column of ranks walks by several", {
  expect_error(design_systematic(3, order_by = c("region", "size")),
               "order(order(data$region, data$size", fixed = TRUE)
  d <- data.frame(region = rep(c("b", "a"), each = 6), size = c(6:1, 1:6))
  d$walk <- order(order(d$region, d$size, method = "radix"))
  des <- design_systematic(interval = 4, order_by = "walk")
  ids <- draw(d, des, seed = 1)
  # The walk runs through region a, then b, each by size.
  walk <- order(d$region, d$size, method = "radix")
  pos <- match(seq_len(12), walk)
  expect_length(unique((pos[match(ids$walk, d$walk)] - 1L) %% 4L), 1L)
  expect_identical(drawn:::systematic_positions(des, d), pos)
})
