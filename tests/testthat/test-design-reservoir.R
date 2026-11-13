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
                 "below the 100 available rows", class = "drawn_stream_truncated")
})

test_that("a simulation over a truncated frame warns once, with the count", {
  d <- data.frame(x = 1:100)
  des <- design_reservoir(n = 5, max_items = 20)
  caught <- list()
  p <- withCallingHandlers(
    inclusion_prob(d, des, simulate = TRUE, R = 300, seed = 1),
    warning = function(w) {
      caught[[length(caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  expect_length(caught, 1L)
  expect_s3_class(caught[[1L]], "drawn_stream_truncated")
  expect_match(conditionMessage(caught[[1L]]),
               "only the first 20. This happened in 300 of the 300 simulated draws.",
               fixed = TRUE)
  expect_lt(max(abs(p[1:20] - 0.25)), 0.08)
  expect_true(all(p[21:100] == 0))
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

test_that("a header read off an open connection is left out of the sample", {
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf), add = TRUE)
  pop <- data.frame(id = 1:5, v = letters[1:5])
  utils::write.csv(pop, tf, row.names = FALSE)

  # Every line is an item, the header included.
  got <- unlist(draw(file(tf), design_reservoir(n = 10), seed = 1))
  expect_identical(got, readLines(tf))

  con <- file(tf, "r")
  on.exit(close(con), add = TRUE, after = FALSE)
  header <- readLines(con, n = 1)
  lines <- draw(con, design_reservoir(n = 10), seed = 1)
  expect_length(lines, 5)
  expect_identical(utils::read.csv(text = c(header, unlist(lines))), pop)
})

test_that("a connection read in blocks gives the line-by-line sample", {
  tf <- tempfile()
  on.exit(unlink(tf), add = TRUE)
  writeLines(as.character(1:25000), tf)
  line_by_line <- function(con) {
    function() {
      line <- readLines(con, n = 1L)
      if (length(line) == 0L) NULL else line
    }
  }
  # 25,000 lines span three blocks of 10,000.
  for (n in c(3, 1500)) {
    con <- file(tf, "r")
    one <- draw(line_by_line(con), design_reservoir(n = n), seed = 4)
    close(con)
    expect_identical(draw(file(tf), design_reservoir(n = n), seed = 4), one)
  }

  # Blocks of 7 lines put swaps at every position within a block, and a cap
  # partway through one.
  for (cap in c(Inf, 100)) {
    con <- file(tf, "r")
    got <- with_seed(8, reservoir_blocks(as_stream(con, block = 7L)$read,
                                         n = 20, max_items = cap))
    close(con)
    ref <- with_seed(8, reservoir_stream(make_counter(25000), n = 20,
                                         max_items = cap))
    expect_identical(as.integer(unlist(got$items)), unlist(ref$items))
    expect_identical(got$truncated, is.finite(cap))
  }
})

test_that("the line read past max_items goes back onto an open connection", {
  tf <- tempfile()
  on.exit(unlink(tf), add = TRUE)
  writeLines(as.character(1:20), tf)

  con <- file(tf, "r")
  on.exit(close(con), add = TRUE, after = FALSE)
  expect_warning(got <- draw(con, design_reservoir(n = 3, max_items = 5),
                             seed = 1), "may have held more")
  expect_true(all(unlist(got) %in% as.character(1:5)))
  # Line 6 was read to see that the stream held more, and pushed back.
  expect_identical(readLines(con, n = 1), "6")
  expect_silent(draw(con, design_reservoir(n = 3, max_items = 14), seed = 1))
  expect_identical(readLines(con), character(0))

  # A connection in binary mode cannot take a line back, so it loses it.
  bin <- file(tf, "rb")
  on.exit(close(bin), add = TRUE, after = FALSE)
  expect_warning(draw(bin, design_reservoir(n = 3, max_items = 5), seed = 1),
                 "may have held more")
  expect_identical(readLines(bin, n = 1), "7")
})

# A generator that returns its arguments one per call, then NULL.
calls <- function(...) {
  out <- list(...)
  i <- 0L
  function() {
    i <<- i + 1L
    if (i > length(out)) NULL else out[[i]]
  }
}

batches_of <- function(pop, size) {
  start <- 0L
  function() {
    if (start >= nrow(pop)) return(NULL)
    rows <- start + seq_len(min(size, nrow(pop) - start))
    start <<- start + size
    pop[rows, , drop = FALSE]
  }
}

test_that("batches = TRUE samples rows across a generator's data frames", {
  pop <- data.frame(id = 1:103, g = rep(c("a", "b"), length.out = 103))
  des <- design_reservoir(n = 8, batches = TRUE)
  got <- draw(batches_of(pop, 10), des, seed = 3)
  expect_identical(got, reindex(pop, got$id))
  expect_equal(nrow(got), 8)
  expect_false(is.unsorted(got$id, strictly = TRUE))

  # The rows are the items, so the sample is the one a stream of the row
  # numbers gives, however the rows are cut into batches.
  expect_identical(got$id, unlist(draw(make_counter(103),
                                       design_reservoir(n = 8), seed = 3)))
  for (size in c(1, 7, 103, 500)) {
    expect_identical(draw(batches_of(pop, size), des, seed = 3), got)
  }

  # Without it, each data frame is one item.
  whole <- draw(batches_of(pop, 10), design_reservoir(n = 3), seed = 3)
  expect_type(whole, "list")
  expect_true(all(vapply(whole, nrow, integer(1)) %in% c(10L, 3L)))
})

test_that("batches = TRUE gives every row the same chance", {
  skip_on_cran()
  pop <- data.frame(id = 1:30)
  des <- design_reservoir(n = 5, batches = TRUE)
  counts <- integer(30)
  for (s in 1:1000) {
    got <- draw(batches_of(pop, 7), des, seed = s)$id
    counts[got] <- counts[got] + 1L
  }
  expect_equal(sum(counts), 5000)
  expect_gt(stats::chisq.test(counts)$p.value, 0.001)
})

test_that("max_items counts rows under batches = TRUE", {
  pop <- data.frame(id = 1:103)
  expect_warning(
    got <- draw(batches_of(pop, 10),
                design_reservoir(n = 4, max_items = 25, batches = TRUE),
                seed = 2),
    "may have held more"
  )
  ref <- suppressWarnings(draw(make_counter(103),
                               design_reservoir(n = 4, max_items = 25),
                               seed = 2))
  expect_identical(got$id, unlist(ref))
  expect_silent(draw(batches_of(pop, 10),
                     design_reservoir(n = 4, max_items = 103, batches = TRUE),
                     seed = 2))
})

test_that("batches = TRUE keeps the batches' class and columns", {
  pop <- data.frame(id = 1:40, f = factor(rep(c("x", "y"), 20)))
  # Empty batches are skipped, and a stream of nothing else keeps the columns.
  empty <- calls(pop[0, ], pop[0, ])
  expect_identical(draw(empty, design_reservoir(n = 3, batches = TRUE)),
                   pop[0, ])
  expect_identical(draw(function() NULL,
                        design_reservoir(n = 3, batches = TRUE)),
                   data.frame())

  skip_if_not_installed("tibble")
  got <- draw(batches_of(tibble::as_tibble(pop), 6),
              design_reservoir(n = 5, batches = TRUE), seed = 1)
  expect_s3_class(got, "tbl_df")
  expect_identical(levels(got$f), c("x", "y"))

  skip_if_not_installed("data.table")
  env <- new.env(parent = globalenv())
  env$s <- draw(batches_of(data.table::as.data.table(pop), 6),
                design_reservoir(n = 5, batches = TRUE), seed = 1)
  expect_s3_class(env$s, "data.table")
  expect_no_warning(evalq(s[, twice := id * 2L], env))
})

test_that("batches = TRUE refuses what is not a batch of rows", {
  expect_error(design_reservoir(n = 3, batches = NA), "`batches` must be")
  des <- design_reservoir(n = 3, batches = TRUE)
  expect_error(draw(list(data.frame(id = 1)), des), "`data` is a list")

  tf <- tempfile()
  on.exit(unlink(tf), add = TRUE)
  writeLines(as.character(1:5), tf)
  con <- file(tf)
  on.exit(close(con), add = TRUE, after = FALSE)
  expect_error(draw(con, des), "`data` is a connection")

  pop <- data.frame(id = 1:3, v = 4:6)
  expect_error(draw(calls(pop, as.matrix(pop)), des),
               "returned a matrix on call 2")
  expect_error(draw(calls(pop, 1:3), des),
               "returned an integer on call 2", fixed = TRUE)
  expect_error(draw(calls(pop, cbind(pop, z = 1)), des),
               "Batch 2 from `data` adds `z`")
  expect_error(draw(calls(pop, pop[2:1]), des),
               "Batch 2 from `data` has the columns `v`, `id`")
})

test_that("a stream of the wrong kind is named with its article", {
  expect_error(draw(1:10, design_reservoir(n = 2)),
               "zero-argument function for a reservoir design, not an integer.",
               fixed = TRUE)
})
