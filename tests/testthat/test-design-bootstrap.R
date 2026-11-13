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

# Each replicate drawn and subset on its own, then bound together, as
# design_bootstrap() built them before it subset the data once.
one_by_one <- function(data, R, n, method = "simple", block_length = NULL) {
  N <- nrow(data)
  parts <- lapply(seq_len(R), function(i) {
    idx <- if (method == "simple") {
      sample.int(N, size = n, replace = TRUE)
    } else {
      starts <- sample.int(N - block_length + 1L, ceiling(n / block_length),
                           replace = TRUE)
      unlist(lapply(starts, function(s) s:(s + block_length - 1L)))[seq_len(n)]
    }
    out <- data[idx, , drop = FALSE]
    rownames(out) <- NULL
    out[[".replicate"]] <- rep(i, nrow(out))
    out[, c(".replicate", names(data)), drop = FALSE]
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

test_that("the stacked rows are those of each replicate drawn on its own", {
  d <- data.frame(id = 1:40, v = (1:40) / 3,
                  s = rep(c("p", "q"), 20),
                  f = factor(rep(c("x", "y", "z", "x"), 10)),
                  t = as.Date("2024-01-01") + 0:39)
  expect_identical(draw(d, design_bootstrap(6, n = 15), seed = 4),
                   with_seed(4, one_by_one(d, 6, 15)))
  expect_identical(
    draw(d, design_bootstrap(5, n = 13, method = "block", block_length = 4),
         seed = 4),
    with_seed(4, one_by_one(d, 5, 13, "block", 4))
  )
  skip_if_not_installed("tibble")
  tb <- tibble::as_tibble(d)
  expect_identical(draw(tb, design_bootstrap(3), seed = 8),
                   with_seed(8, one_by_one(tb, 3, 40)))
})

test_that("indices and weights hold the replicates the rows hold", {
  d <- data.frame(id = 1:12, y = (1:12)^2)
  rows <- draw(d, design_bootstrap(5, n = 8), seed = 3)
  idx <- draw(d, design_bootstrap(5, n = 8, output = "indices"), seed = 3)
  w <- draw(d, design_bootstrap(5, n = 8, output = "weights"), seed = 3)

  expect_true(is.matrix(idx))
  expect_type(idx, "integer")
  expect_identical(dim(idx), c(8L, 5L))
  expect_identical(as.vector(idx), rows$id)

  expect_type(w, "integer")
  expect_identical(dim(w), c(12L, 5L))
  expect_identical(w, apply(idx, 2, tabulate, nbins = 12))
  expect_equal(colSums(w), rep(8, 5))

  # One row per replicate is still a matrix.
  expect_identical(dim(draw(d, design_bootstrap(3, n = 1, output = "indices"),
                            seed = 1)), c(1L, 3L))
})

test_that("bootstrap counts are replicate weights for survey", {
  skip_if_not_installed("survey")
  d <- data.frame(y = c(3, 8, 1, 9, 4, 7, 2, 6, 5, 10, 12, 0, 11, 15, 13))
  w <- draw(d, design_bootstrap(60, output = "weights"), seed = 2)
  reps <- survey::svrepdesign(data = d, repweights = w,
                              weights = rep(1, nrow(d)), type = "bootstrap",
                              combined.weights = TRUE)
  means <- colSums(w * d$y) / colSums(w)
  est <- survey::svymean(~y, reps)
  expect_equal(unname(coef(est)), mean(d$y))
  expect_equal(unname(survey::SE(est))[1], stats::sd(means))
})

test_that("circular runs wrap around and reach every row equally often", {
  d <- data.frame(id = 1:10)
  idx <- draw(d, design_bootstrap(200, n = 8, method = "circular",
                                  block_length = 4, output = "indices"),
              seed = 5)
  runs <- matrix(idx, nrow = 4)
  expect_true(all((diff(runs) %% 10L) == 1L))
  expect_true(any(diff(runs) == -9L))

  d30 <- data.frame(id = 1:30)
  counts <- function(method) {
    rowSums(draw(d30, design_bootstrap(2000, n = 30, method = method,
                                       block_length = 5, output = "weights"),
                 seed = 7))
  }
  circ <- counts("circular")
  expect_true(all(abs(circ / 2000 - 1) < 0.1))
  # The non-circular moving-block bootstrap reaches its end rows less often.
  block <- counts("block")
  expect_lt(block[1], 0.4 * 2000)
  expect_lt(block[30], 0.4 * 2000)
})

test_that("the default block length is floor(N^(1/3)) at perfect cubes too", {
  des <- design_bootstrap(1, method = "block")
  expect_identical(bootstrap_block_length(des, 64L), 4L)
  expect_identical(bootstrap_block_length(des, 63L), 3L)
  expect_identical(bootstrap_block_length(des, 1000L), 10L)
  expect_identical(bootstrap_block_length(des, 2L), 1L)
  expect_identical(
    bootstrap_block_length(design_bootstrap(1, method = "block",
                                            block_length = 7), 64L), 7L)

  idx <- draw(data.frame(id = 1:64),
              design_bootstrap(1, n = 64, method = "block",
                               output = "indices"), seed = 1)
  expect_true(all(diff(matrix(idx, nrow = 4)) == 1L))
})

test_that("duplicated column names come back as they went in", {
  dd <- data.frame(a = 1:5, a = 6:10, check.names = FALSE)
  out <- draw(dd, design_bootstrap(2, n = 3), seed = 1)
  expect_identical(names(out), c(".replicate", "a", "a"))
  expect_identical(out[[3]] - out[[2]], rep(5L, 6))

  skip_if_not_installed("tibble")
  tb <- tibble::as_tibble(dd, .name_repair = "minimal")
  out_tb <- draw(tb, design_bootstrap(2, n = 3), seed = 1)
  expect_s3_class(out_tb, "tbl_df")
  expect_identical(names(out_tb), c(".replicate", "a", "a"))
  expect_identical(out_tb[[3]] - out_tb[[2]], rep(5L, 6))
})

test_that("n = 0 gives no rows but keeps the .replicate column", {
  d <- data.frame(id = 1:5)
  out <- draw(d, design_bootstrap(3, n = 0), seed = 1)
  expect_identical(names(out), c(".replicate", "id"))
  expect_equal(nrow(out), 0)
  expect_type(out$.replicate, "integer")
  expect_length(split(out, factor(out$.replicate, 1:3)), 3)

  expect_identical(dim(draw(d, design_bootstrap(3, n = 0, output = "indices"),
                            seed = 1)), c(0L, 3L))
  w <- draw(d, design_bootstrap(3, n = 0, output = "weights"), seed = 1)
  expect_identical(w, matrix(0L, 5, 3))
})

test_that("a stacked result past 10 million rows warns with its own class", {
  d <- data.frame(id = 1:5)
  cond <- tryCatch(draw(d, design_bootstrap(2, n = 5000001), seed = 1),
                   drawn_bootstrap_size = function(w) w)
  expect_s3_class(cond, "drawn_bootstrap_size")
  expect_s3_class(cond, "drawn_warning")
  expect_match(conditionMessage(cond), "10,000,002 rows", fixed = TRUE)
  expect_match(conditionMessage(cond), "output = \"indices\"", fixed = TRUE)

  skip_on_cran()
  expect_no_warning(
    draw(d, design_bootstrap(2, n = 5000001, output = "indices"), seed = 1)
  )
})

test_that("the refusals give the expected count the method gives", {
  d <- data.frame(id = 1:20)
  msg <- function(des, ...) {
    tryCatch(inclusion_prob(d, des, ...), error = conditionMessage)
  }
  for (m in c("simple", "circular")) {
    expect_match(msg(design_bootstrap(5, method = m)),
                 "that is n_replicates * n / nrow(data)", fixed = TRUE)
  }
  block <- design_bootstrap(5, method = "block")
  for (got in list(msg(block), msg(block, simulate = TRUE, R = 10))) {
    expect_match(got, "first and last block_length - 1 rows appear less often",
                 fixed = TRUE)
    expect_no_match(got, "that is n_replicates", fixed = TRUE)
  }
})

test_that("method and output are checked when the design is built", {
  expect_error(design_bootstrap(5, output = "matrix"),
               "^`output` must be one of \"rows\", \"indices\" or \"weights\"")
  expect_error(design_bootstrap(5, block_length = 3),
               "only applies to method = \"block\" or \"circular\"")
  expect_silent(design_bootstrap(5, method = "circular", block_length = 3))
})

test_that("a selection plot draws the rows whatever the output", {
  d <- data.frame(id = 1:30)
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit({grDevices::dev.off(); unlink(f)}, add = TRUE)
  for (out in c("rows", "indices", "weights")) {
    expect_no_error(plot(design_bootstrap(2, n = 4, output = out), d,
                         type = "selection", seed = 1))
  }
})

test_that("reindex() subsets a plain data frame as `[` does", {
  d <- data.frame(i = 1:6, s = letters[1:6],
                  f = factor(c("a", "b", "a", "c", "b", "a"),
                             levels = c("c", "b", "a", "z")),
                  t = as.POSIXct("2024-01-01", tz = "America/New_York") +
                    3600 * 0:5)
  d$l <- I(list(1, "a", NULL, 1:3, NA, list(2)))
  d$m <- matrix(1:12, nrow = 6, dimnames = list(letters[1:6], c("p", "q")))
  d$df <- data.frame(u = 6:1)
  names(d)[2] <- "i"
  rownames(d) <- paste0("r", 1:6)
  attr(d, "note") <- "kept"
  for (idx in list(c(2L, 2L, 5L, 1L), integer(0), 6:1)) {
    ref <- d[idx, , drop = FALSE]
    rownames(ref) <- NULL
    expect_identical(reindex(d, idx), ref)
  }
})
