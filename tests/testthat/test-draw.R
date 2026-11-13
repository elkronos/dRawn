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
  set.seed(42)
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
    set.seed(42)
    invisible(draw(df, d, seed = 7))
    expect_equal(stats::runif(3), expected, info = class(d)[1])
  }
})

test_that("an unseeded draw still advances the RNG normally", {
  df <- make_df()
  set.seed(1)
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
