test_that("draw() rejects anything that is not a design", {
  df <- make_df()
  expect_error(draw(df, "stratified"), "must come from one of the design")
  expect_error(draw(df, list(n = 5)), "must come from one of the design")
})

test_that("is_design() recognises designs and nothing else", {
  expect_true(is_design(design_simple(n = 1)))
  expect_true(is_design(design_bootstrap()))
  expect_false(is_design(mtcars))
  expect_false(is_design(NULL))
})

test_that("designs print their parameters", {
  out <- capture.output(print(design_stratified("site", n = 20)))
  expect_match(out[1], "sampling design: stratified")
  expect_true(any(grepl("strata", out)))
  expect_true(any(grepl("proportional", out)))
})

test_that("a design is reusable across data sets", {
  d <- design_stratified("site", n = 8)
  df <- make_df()
  expect_equal(nrow(draw(df, d, seed = 1)), 8)
  expect_equal(nrow(draw(df[1:60, ], d, seed = 1)), 8)
})

test_that("every design returns the input's class and column order", {
  df <- make_df()
  df$w <- abs(df$value) + 1
  results <- list(
    draw(df, design_simple(n = 4), seed = 1),
    draw(df, design_stratified("site", n = 4), seed = 1),
    draw(df, design_systematic(interval = 10), seed = 1),
    draw(df, design_cluster("site", n_clusters = 2), seed = 1),
    draw(df, design_multistage("site", n_clusters = 2, n = 4), seed = 1),
    draw(df, design_weighted("w", n = 3), seed = 1),
    draw(df, design_reservoir(n = 4), seed = 1),
    draw(df, design_temporal("time", "2020-01-01 00:00:00",
                             "2020-01-01 12:00:00", interval = 6,
                             per_interval = 1, unit = "hours"), seed = 1)
  )
  for (res in results) expect_same_schema(res, df)
})

test_that("seeding is local to the draw", {
  df <- make_df()
  local_seed(42)
  expected <- stats::runif(3)

  designs <- list(
    design_simple(n = 2),
    design_stratified("site", n = 4),
    design_systematic(interval = 10),
    design_cluster("site", n_clusters = 2),
    design_multistage("site", n_clusters = 2, n = 4),
    design_reservoir(n = 3),
    design_bootstrap(n_replicates = 2, n = 3),
    design_temporal("time", "2020-01-01 00:00:00", "2020-01-01 12:00:00",
                    interval = 6, per_interval = 1, unit = "hours")
  )
  for (d in designs) {
    local_seed(42)
    invisible(draw(df, d, seed = 7))
    expect_equal(stats::runif(3), expected, info = class(d)[1])
  }
})

test_that("an unseeded draw still advances the RNG normally", {
  df <- make_df()
  local_seed(1)
  a <- draw(df, design_simple(n = 3))$id
  b <- draw(df, design_simple(n = 3))$id
  expect_false(identical(a, b))
})

test_that("the same seed gives the same sample", {
  df <- make_df()
  d <- design_simple(n = 10)
  expect_equal(draw(df, d, seed = 123), draw(df, d, seed = 123))
})

test_that("an empty or non-data-frame input is reported clearly", {
  expect_error(draw(data.frame(x = integer(0)), design_simple(n = 1)),
               "has no rows")
  expect_error(draw(1:10, design_simple(n = 1)), "must be a data frame")
})

test_that("weights = TRUE stores a copy of the frame and its probabilities", {
  df <- make_df()
  des <- design_stratified("site", n = 8)
  s <- draw(df, des, seed = 1, weights = TRUE)
  expect_identical(attr(s, "drawn_population"), df)
  expect_equal(attr(s, "drawn_pi"), inclusion_prob(df, des))
  expect_identical(attr(s, "drawn_rows"), s$.row)
  expect_type(s$.row, "integer")
  expect_equal(s$.prob, attr(s, "drawn_pi")[s$.row])
})

test_that("the frame copy keeps a tibble's and an sf frame's class", {
  skip_if_not_installed("tibble")
  tb <- tibble::as_tibble(make_df())
  s <- draw(tb, design_simple(n = 5), seed = 1, weights = TRUE)
  expect_identical(attr(s, "drawn_population"), tb)

  skip_if_not_installed("sf")
  pts <- sf::st_as_sf(make_df(), coords = c("lon", "lat"), crs = 4326,
                      remove = FALSE)
  s <- draw(pts, design_simple(n = 5), seed = 1, weights = TRUE)
  expect_s3_class(s, "sf")
  expect_identical(attr(s, "drawn_population"), pts)
})

test_that("an in-place edit of a data.table frame does not reach the sample", {
  skip_if_not_installed("data.table")
  dt <- data.table::as.data.table(make_df())
  data.table::setkeyv(dt, "site")
  s <- draw(dt, design_simple(n = 5), seed = 1, weights = TRUE)
  pop <- attr(s, "drawn_population")
  expect_s3_class(pop, "data.table")
  expect_identical(data.table::key(pop), "site")
  before <- pop$value
  data.table::set(dt, j = "value", value = 0)
  data.table::setorderv(dt, "id", order = -1L)
  expect_identical(attr(s, "drawn_population")$value, before)
})

test_that("renaming, reordering or deleting frame columns does not reach the sample", {
  # The copy of the frame shared its names vector with the caller's, which
  # setnames(), setcolorder() and deleting a column rewrite in place.
  skip_if_not_installed("data.table")
  for (as_table in c(FALSE, TRUE)) {
    d <- make_df()
    if (as_table) d <- data.table::as.data.table(d)
    s <- draw(d, design_stratified("site", n = 8), seed = 1, weights = TRUE)
    want <- ht_total(s, "value")
    cols <- names(make_df())
    data.table::setnames(d, "site", "SITE")
    data.table::setcolorder(d, rev(names(d)))
    data.table::set(d, j = "lon", value = NULL)
    pop <- attr(s, "drawn_population")
    expect_identical(names(pop), cols)
    expect_identical(pop$site, make_df()$site)
    expect_equal(ht_total(s, "value"), want)
  }
})

test_that("a sample drawn from a data.table takes := without a warning", {
  skip_if_not_installed("data.table")
  dt <- data.table::as.data.table(make_df())
  for (w in c(FALSE, TRUE)) {
    # `:=` works only where data.table is attached or imported, which the
    # global environment counts as.
    env <- new.env(parent = globalenv())
    env$s <- draw(dt, design_stratified("site", n = 8), seed = 1, weights = w)
    expect_s3_class(env$s, "data.table")
    expect_no_warning(evalq(s[, obs := value * 2], env))
    expect_equal(env$s$obs, env$s$value * 2)
    if (w) expect_false(is.null(attr(env$s, "drawn_design")))
  }
})

test_that("a value outside an argument's choices is refused by name", {
  # match.arg() said "'arg' should be one of", which names no argument.
  msg <- function(expr) tryCatch(expr, error = conditionMessage)
  expect_equal(
    msg(design_stratified("g", n = 30, allocation = "neymann")),
    paste0("`allocation` must be one of \"proportional\", \"equal\" or ",
           "\"neyman\", not \"neymann\".")
  )
  expect_match(msg(design_multistage("g", n_clusters = 2, n = 4,
                                     allocation = "even")),
               "^`allocation` must be one of \"equal\" or \"proportional\"")
  expect_match(msg(design_weighted("w", n = 5, method = "pps")),
               "^`method` must be one of")
  expect_match(msg(design_bootstrap(5, method = "wild")),
               "^`method` must be one of")
  expect_match(msg(design_temporal("t", "2020-01-01", "2020-01-02", 1, 1,
                                   unit = "fortnights")),
               "^`unit` must be one of")
  expect_match(msg(plan_size(0.1, sd = 1, target = "median")),
               "^`target` must be one of")
  expect_match(msg(plot(design_simple(n = 2), data.frame(x = 1:5),
                        type = "pie")),
               "^`type` must be one of")
  s <- draw(data.frame(y = 1:10), design_simple(n = 5), seed = 1,
            weights = TRUE)
  expect_match(msg(ht_total(s, "y", variance = "boot")),
               "^`variance` must be one of")
  expect_match(msg(ht_mean(s, "y", estimator = "median")),
               "^`estimator` must be one of")
  expect_match(msg(ht_mean(s, "y", variance = NA)), "not NA\\.$")

  # The default and a unique abbreviation work as they did with match.arg().
  expect_equal(design_stratified("g", n = 30)$allocation, "proportional")
  expect_equal(design_stratified("g", n = 30, allocation = "ney",
                                 allocation_by = "w")$allocation, "neyman")
  expect_equal(ht_total(s, "y", variance = "jack")$method, "jackknife")
})

test_that("a printed design says what a NULL default means", {
  out <- capture.output(print(design_bootstrap(20, method = "block")))
  expect_true(any(grepl("^  n +nrow\\(data\\)$", out)))
  expect_true(any(grepl("^  block_length +nrow\\(data\\)\\^\\(1/3\\), rounded down$",
                        out)))
  out <- capture.output(print(design_bootstrap(20)))
  expect_false(any(grepl("block_length", out)))
  out <- capture.output(print(design_systematic(interval = 5)))
  expect_true(any(grepl("^  start +random$", out)))
  out <- capture.output(print(design_systematic(interval = 5, start = 2)))
  expect_true(any(grepl("^  start +2$", out)))
})

test_that("messages give a class its article", {
  expect_error(draw(data.frame(id = 1:5), 1L),
               "constructors, not an integer.", fixed = TRUE)
  expect_error(draw(1:5, design_simple(n = 2), weights = TRUE),
               "needs a data frame, not an integer.", fixed = TRUE)
  expect_error(design_certainty("w", 1, list()),
               "not a list.", fixed = TRUE)
})

# ---- seeds, generators and the record of a draw -----------------------------

# Run `code` under another generator, then put back the caller's random number
# state, or its absence, and generator.
under_kind <- function(kind, code) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv())
  old_kind <- RNGkind()
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else {
      do.call(RNGkind, as.list(old_kind))
      suppressWarnings(rm(".Random.seed", envir = globalenv()))
    }
  })
  do.call(RNGkind, as.list(kind))
  code
}

test_that("a seed draws the same sample under any generator", {
  df <- make_df()
  des <- design_stratified("site", n = 8)
  ref <- draw(df, des, seed = 11)
  for (kind in list("L'Ecuyer-CMRG", "Knuth-TAOCP-2002",
                    list("Mersenne-Twister", "Box-Muller", "Rounding"))) {
    got <- suppressWarnings(under_kind(kind, {
      before <- RNGkind()
      out <- draw(df, des, seed = 11)
      # The caller's generator is back in place afterwards.
      expect_identical(RNGkind(), before)
      out
    }))
    expect_identical(got, ref)
  }
  expect_equal(inclusion_prob(df, design_simple(n = 8), simulate = TRUE,
                              R = 50, seed = 3),
               under_kind("L'Ecuyer-CMRG",
                          inclusion_prob(df, design_simple(n = 8),
                                         simulate = TRUE, R = 50, seed = 3)))
})

test_that("a seeded draw leaves the caller's generator and stream alone", {
  df <- make_df()
  under_kind("L'Ecuyer-CMRG", {
    before <- get(".Random.seed", envir = globalenv())
    draw(df, design_simple(n = 3), seed = 1)
    expect_identical(get(".Random.seed", envir = globalenv()), before)
    expect_identical(RNGkind()[[1L]], "L'Ecuyer-CMRG")
  })
  # With no stream started, none is left behind, and the generator is the
  # caller's.
  under_kind("L'Ecuyer-CMRG", {
    rm(".Random.seed", envir = globalenv())
    draw(df, design_simple(n = 3), seed = 1)
    expect_false(exists(".Random.seed", envir = globalenv(),
                        inherits = FALSE))
    expect_identical(RNGkind()[[1L]], "L'Ecuyer-CMRG")
  })
})

test_that("a seed outside the integer range or given as text is refused", {
  df <- make_df()
  des <- design_simple(n = 3)
  for (s in list(2^31, -2^31, 1e15)) {
    expect_error(draw(df, des, seed = s), "must lie between -2147483647 and")
  }
  expect_no_error(draw(df, des, seed = .Machine$integer.max))
  expect_error(draw(df, des, seed = "1"), 'not the string "1"', fixed = TRUE)
  # The probability functions check it whether or not they simulate.
  expect_error(inclusion_prob(df, des, seed = "abc"), 'the string "abc"',
               fixed = TRUE)
  expect_error(joint_prob(df, des, rows = 1:2, seed = 2^31), "must lie")
})

test_that("a weighted sample records how it was drawn", {
  df <- make_df()
  des <- design_stratified("site", n = 8)
  s <- draw(df, des, seed = 5, weights = TRUE)
  meta <- attr(s, "drawn_meta")
  expect_identical(meta$seed, 5)
  expect_identical(meta$rng_kind,
                   c("Mersenne-Twister", "Inversion", "Rejection"))
  expect_identical(meta$version, as.character(utils::packageVersion("drawn")))
  expect_s3_class(meta$time, "POSIXct")
  expect_identical(meta$frame_rows, 100L)
  expect_identical(meta$made_by, "draw()")
  expect_identical(sample_summary(s)$drawn, meta)
  expect_match(capture.output(sample_summary(s)),
               "with seed 5, by drawn [0-9.]+ on [0-9-]+$", all = FALSE)

  # Without a seed, the session's generator is recorded.
  u <- under_kind("L'Ecuyer-CMRG", draw(df, des, weights = TRUE))
  expect_null(attr(u, "drawn_meta")$seed)
  expect_identical(attr(u, "drawn_meta")$rng_kind[[1L]], "L'Ecuyer-CMRG")
  expect_match(capture.output(sample_summary(u)),
               "without a seed (L'Ecuyer-CMRG", fixed = TRUE, all = FALSE)

  # A sample rebuilt from a file does not know its seed.
  back <- as_drawn(as.data.frame(s)[, c(".row", names(df))], des, df)
  expect_identical(attr(back, "drawn_meta")$made_by, "as_drawn()")
  expect_null(attr(back, "drawn_meta")$seed)
  expect_match(capture.output(sample_summary(back)), "rebuilt",
               all = FALSE)
})

test_that("use_reserves() keeps the record of the draw", {
  df <- make_df()
  l <- draw(df, design_simple(n = 4, reserve = 2), seed = 9, weights = TRUE)
  s <- use_reserves(l, c(TRUE, FALSE, FALSE, FALSE, FALSE, NA))
  expect_identical(attr(s, "drawn_meta"), attr(l, "drawn_meta"))
})

# ---- row names --------------------------------------------------------------

test_that("row names of the frame's own stay with their rows", {
  cars <- mtcars[, c("mpg", "cyl", "wt")]
  cars$id <- seq_len(nrow(cars))
  s <- draw(cars, design_simple(n = 5), seed = 1)
  expect_identical(rownames(s), rownames(cars)[s$id])
  w <- draw(cars, design_stratified("cyl", n = 6), seed = 1, weights = TRUE)
  expect_identical(rownames(w), rownames(cars)[w$.row])
  # A frame left with integer row names by an earlier subset keeps them.
  sub <- cars[cars$cyl != 6, ]
  w <- draw(sub, design_simple(n = 4), seed = 2, weights = TRUE)
  expect_identical(rownames(w), rownames(sub)[w$.row])

  # A row drawn twice gets a name made unique, as `[` gives it.
  r <- draw(cars, design_simple(n = 40, replace = TRUE), seed = 1)
  expect_false(anyDuplicated(rownames(r)) > 0)
  expect_identical(rownames(r), make.unique(rownames(cars)[r$id]))
  b <- draw(cars, design_bootstrap(n_replicates = 2, n = 5), seed = 1)
  expect_false(anyDuplicated(rownames(b)) > 0)
})

test_that("automatic row names are numbered afresh, also after rows drop", {
  df <- make_df()
  df$o <- df$id
  df$o[c(3, 50)] <- NA
  df$w <- df$id
  df$w[7] <- NA
  df$cl <- rep(1:10, each = 10)
  df$cl[c(1, 95)] <- NA
  df$big <- df$id
  df$big[4] <- NA
  designs <- list(
    design_systematic(interval = 7, order_by = "o", na_rm = TRUE),
    design_weighted("w", n = 6, na_rm = TRUE),
    design_cluster("cl", n_clusters = 2, na_rm = TRUE),
    design_certainty("big", 95, rest = design_simple(n = 3), na_rm = TRUE)
  )
  for (d in designs) {
    s <- draw(df, d, seed = 4)
    expect_identical(attr(s, "row.names"), seq_len(nrow(s)),
                     info = class(d)[1])
  }
})

# ---- arguments --------------------------------------------------------------

test_that("a design and data given the wrong way round are named as swapped", {
  df <- make_df()
  des <- design_simple(n = 3)
  expect_error(draw(des, df), "look swapped. Call draw(data, design)",
               fixed = TRUE)
  expect_error(des |> draw(df), "look swapped")
  expect_error(inclusion_prob(des, df), "Call inclusion_prob(data, design)",
               fixed = TRUE)
  expect_error(joint_prob(des, df), "Call joint_prob(data, design)",
               fixed = TRUE)
})

test_that("a design class of the user's own draws through draw_design()", {
  df <- data.frame(id = 1:10)
  first <- structure(list(n = 3), class = c("drawn_test_first", "drawn_design"))
  .S3method("draw_design", "drawn_test_first", function(design, data) {
    data[seq_len(design$n), , drop = FALSE]
  })
  expect_identical(draw(df, first)$id, 1:3)
  expect_equal(inclusion_prob(df, first, simulate = TRUE, R = 5, seed = 1),
               rep(c(1, 0), c(3, 7)))
  expect_error(draw(df, first, weights = TRUE),
               "simulate = TRUE) estimates the probabilities", fixed = TRUE)
  expect_error(joint_prob(df, first, rows = 1:2),
               "joint_prob(data, design, simulate = TRUE)", fixed = TRUE)
  # A method that drops columns is caught rather than counted as drawing
  # nothing.
  thin <- structure(list(), class = c("drawn_test_thin", "drawn_design"))
  .S3method("draw_design", "drawn_test_thin", function(design, data) {
    data.frame(id = data$id[1:2])
  })
  expect_error(inclusion_prob(df, thin, simulate = TRUE, R = 2, seed = 1),
               "returned rows without the `.drawn_row_id` column")
})
