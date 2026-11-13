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

test_that("a connection that was not open is read before it is closed", {
  tf <- tempfile()
  on.exit(unlink(tf), add = TRUE)
  writeLines(as.character(1:20), tf)
  con <- file(tf)                       # created, not opened
  got <- unlist(draw(con, design_reservoir(n = 3), seed = 1))
  expect_length(got, 3)
  expect_true(all(got %in% as.character(1:20)))
})

test_that("a stream of exactly max_items is not reported as truncated", {
  make_gen <- function(k) {
    i <- 0L
    function() {
      i <<- i + 1L
      if (i > k) NULL else i
    }
  }
  expect_silent(draw(make_gen(5), design_reservoir(n = 3, max_items = 5),
                     seed = 1))
  expect_warning(draw(make_gen(9), design_reservoir(n = 3, max_items = 5),
                      seed = 1), "may have held more")
})

make_counter <- function(k) {
  i <- 0L
  function() {
    i <<- i + 1L
    if (i > k) NULL else i
  }
}

test_that("a NULL element of a list is an item, not the end of the stream", {
  expect_identical(draw(list(1, NULL, 3, 4, 5), design_reservoir(n = 10),
                        seed = 1),
                   list(1, NULL, 3, 4, 5))

  recs <- list(a = "a", b = NULL, c = "c", d = "d", e = "e", f = "f")
  seen <- unlist(lapply(1:200, function(s) {
    names(draw(recs, design_reservoir(n = 2), seed = s))
  }))
  expect_true(all(table(factor(seen, letters[1:6])) > 40))
})

test_that("a list is sampled by position, as a data frame is", {
  expect_identical(
    unlist(draw(as.list(1:50), design_reservoir(n = 5), seed = 9)),
    draw(data.frame(x = 1:50), design_reservoir(n = 5), seed = 9)$x
  )
})

test_that("the stream counter is a double, which does not overflow at 2^31", {
  out <- reservoir_stream(make_counter(10), n = 2, max_items = Inf)
  expect_type(out$count, "double")
  expect_equal(out$count, 10)
})

test_that("a stream's items come back in the order they arrived", {
  for (s in 1:5) {
    got <- unlist(draw(make_counter(500), design_reservoir(n = 20), seed = s))
    expect_length(got, 20)
    expect_false(is.unsorted(got, strictly = TRUE))
  }
})

test_that("growing the reservoir keeps the items Algorithm L keeps", {
  # The reservoir set aside all n slots at the start. Growing it as the stream
  # fills it, and returning items in arrival order, must not change which
  # items are kept for a seed.
  algorithm_l <- function(items, n) {
    res <- vector("list", n)
    w <- exp(log(stats::runif(1)) / n)
    nxt <- n + floor(log(stats::runif(1)) / log(1 - w)) + 1
    for (k in seq_along(items)) {
      if (k <= n) {
        res[[k]] <- items[[k]]
      } else if (k == nxt) {
        res[[sample.int(n, 1L)]] <- items[[k]]
        w <- w * exp(log(stats::runif(1)) / n)
        nxt <- nxt + floor(log(stats::runif(1)) / log(1 - w)) + 1
      }
    }
    unlist(res[seq_len(min(n, length(items)))])
  }
  for (case in list(c(n = 3, N = 100), c(n = 1500, N = 4000),
                    c(n = 2000, N = 1800))) {
    got <- unlist(draw(make_counter(case[["N"]]),
                       design_reservoir(n = case[["n"]]), seed = 2))
    ref <- with_seed(2, algorithm_l(seq_len(case[["N"]]), case[["n"]]))
    expect_identical(got, sort(ref))
  }
})

test_that("a large n costs no memory for a short stream", {
  skip_on_cran()
  # All n slots were set aside before the first item was read, so n = 2e7
  # took 150 MB for a stream of five items.
  before <- gc(reset = TRUE)[2, 6]
  got <- draw(make_counter(5), design_reservoir(n = 2e7), seed = 1)
  expect_lt(gc()[2, 6] - before, 50)
  expect_identical(unlist(got), 1:5)
})
