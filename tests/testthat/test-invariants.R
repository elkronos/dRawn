# Behaviour that is easy to get wrong and easy to regress: base R's sharper
# edges, the cases where a design must refuse rather than guess, and the
# invariants every design is supposed to hold to.
# Each pins a specific case where a plausible implementation gives a wrong
# answer silently, rather than erroring.

test_that("a single numeric cluster id is not reinterpreted as a range", {
  # sample(5, 1) samples from 1:5, so the only real cluster was never selected.
  d <- data.frame(x = 1:10, cl = rep(5, 10))
  expect_equal(nrow(draw(d, design_cluster("cl", n_clusters = 1), seed = 1)), 10)
  expect_equal(
    nrow(draw(d, design_multistage("cl", n_clusters = 1, n = 2), seed = 1)), 2
  )
})

test_that("capping a result keeps every stratum represented", {
  # max_rows was head() of a stratum-sorted frame, so a cap of 10 returned
  # group "a" and nothing else. The argument is gone: n IS the cap, and it is
  # allocated across strata rather than applied afterwards.
  df <- make_df()
  res <- draw(df, design_stratified("site", n = 10), seed = 1)
  expect_equal(nrow(res), 10)
  expect_equal(length(unique(res$site)), 4)

  seen <- unique(unlist(lapply(1:20, function(i) {
    draw(df, design_stratified("site", n = 10), seed = i)$site
  })))
  expect_setequal(seen, letters[1:4])
})

test_that("a temporal sample covers the whole window (invariant, not a guard)", {
  # The temporal half of A2 lived entirely in max_rows, which is gone, so this
  # cannot regress the way it did. It passes against the original code too.
  # Kept as an invariant; the actual guard against A2 returning is the
  # "removed arguments stay removed" test at the end of this file.
  res <- draw(make_ts(), design_temporal("ts", "2020-01-01", "2020-01-02",
                                         interval = 6, per_interval = 1,
                                         unit = "hours"), seed = 1)
  expect_equal(nrow(res), 4)
  expect_gt(as.numeric(difftime(max(res$ts), min(res$ts), units = "hours")), 6)
})

test_that("multistage allocation totals the request, up or down", {
  df <- make_df()
  # Used to return 10 for a request of 20 (denominator was the whole population).
  expect_equal(
    nrow(draw(df, design_multistage("site", n_clusters = 2, n = 20,
                                    allocation = "proportional"), seed = 1)),
    20
  )
  # And 20 for a request of 5, because max(1, ...) floored every cluster.
  many <- data.frame(x = 1:1000, cl = as.character(rep(1:100, each = 10)))
  expect_equal(
    nrow(draw(many, design_multistage("cl", n_clusters = 20, n = 5,
                                      allocation = "proportional"), seed = 1)),
    5
  )
})

test_that("the caller's weights come back untouched", {
  d <- data.frame(id = 1:10, w = 1:10)
  out <- draw(d, design_weighted("w", n = 5), seed = 1)
  expect_equal(out$w, d$w[match(out$id, d$id)])
  expect_type(out$w, "integer")
})

test_that("the lightest row stays selectable", {
  # min-max and z-score both subtracted the minimum and added double.eps,
  # leaving the lightest row at ~2.2e-16. Rescaling is gone entirely: prob=
  # normalises by sum, so it never bought anything.
  d <- data.frame(id = 1:5, w = c(1, 2, 3, 4, 5))
  drawn <- vapply(1:600, function(i) draw(d, design_weighted("w", n = 1),
                                          seed = i)$id, numeric(1))
  expect_true(1 %in% drawn)
})

test_that("balanced clustering works with unused factor levels", {
  d <- data.frame(x = 1:30,
                  cl = factor(rep(c("a", "b", "c"), each = 10),
                              levels = c("a", "b", "c", "zz")))
  res <- draw(d, design_cluster("cl", n_clusters = 2, balanced = TRUE), seed = 1)
  expect_equal(nrow(res), 20)
})

test_that("block replicates are stitched from blocks across the frame", {
  d <- data.frame(id = 1:1000)
  res <- draw(d, design_bootstrap(n_replicates = 20, n = 50, method = "block",
                                  block_length = 5), seed = 1)
  reps <- split(res, res$.replicate)
  spans <- vapply(reps, function(z) diff(range(z$id)), numeric(1))
  expect_true(all(spans > 200))
})

test_that("bootstrap returns exactly n_replicates of equal size (invariant)", {
  # A8 was the max_rows branch rbinding, de-duplicating and re-splitting the
  # replicates. That branch is gone, so this passes against the original code
  # too. Kept as an invariant, not as a regression guard.
  df <- make_df()
  res <- draw(df, design_bootstrap(n_replicates = 5, n = 10), seed = 1)
  expect_equal(length(unique(res$.replicate)), 5)
  expect_equal(unique(as.vector(table(res$.replicate))), 10L)
})

test_that("POSIXct and Date columns are used as-is", {
  expect_silent(
    res <- draw(make_ts(), design_temporal("ts", "2020-01-01", "2020-01-02",
                                           interval = 6, per_interval = 2,
                                           unit = "hours"), seed = 1)
  )
  expect_equal(nrow(res), 8)

  dates <- data.frame(id = 1:10,
                      d = seq(as.Date("2020-01-01"), by = "day", length.out = 10))
  expect_equal(
    nrow(draw(dates, design_temporal("d", "2020-01-01", "2020-01-06",
                                     interval = 1, per_interval = 1,
                                     unit = "days"), seed = 1)),
    5
  )
})

test_that("a region list of three or more geometries is unioned", {
  skip_if_not_installed("sf")
  d <- data.frame(id = 1:3, lon = c(5, 25, 45), lat = c(5, 25, 45))
  expect_equal(
    nrow(draw(d, design_spatial(c("lon", "lat"), three_boxes(), n = 3), seed = 1)),
    3
  )
})

test_that("count arguments reject negative, non-numeric and length-2 values", {
  # max_rows accepted all three: head(x, -3) trimmed from the end, "5" was
  # ignored, and c(3, 2) silently dropped columns.
  expect_error(design_simple(n = -3), "non-negative")
  expect_error(design_simple(n = "5"), "whole number")
  expect_error(design_simple(n = c(3, 2)), "whole number")
  expect_error(design_stratified("g", n = -1), "non-negative")
  expect_error(design_bootstrap(n_replicates = 0), "at least 1")
})

test_that("single-column frames survive intact", {
  d <- data.frame(cl = rep(c("a", "b"), each = 5))
  res <- draw(d, design_cluster("cl", n_clusters = 1), seed = 1)
  expect_identical(names(res), "cl")
  expect_equal(nrow(res), 5)

  res_ms <- draw(d, design_multistage("cl", n_clusters = 1, n = 2), seed = 1)
  expect_identical(names(res_ms), "cl")

  one_col <- data.frame(v = c(3, 1, 2))
  res_sys <- draw(one_col, design_systematic(interval = 1, order_by = "v",
                                             start = 1))
  expect_identical(names(res_sys), "v")
  expect_equal(res_sys$v, c(1, 2, 3))
})

test_that("a short stream is trimmed, not padded with NULLs", {
  res <- draw(as.list(1:5), design_reservoir(n = 10), seed = 1)
  expect_equal(length(res), 5)
  expect_false(any(vapply(res, is.null, logical(1))))
})

test_that("empty results keep the input's columns", {
  df <- make_df()
  empty <- draw(df, design_temporal("time", "2021-01-01", "2021-01-02",
                                    interval = 6, per_interval = 2,
                                    unit = "hours"), seed = 1)
  expect_equal(ncol(empty), ncol(df))
  expect_identical(names(empty), names(df))
})

test_that("NA cluster labels never fabricate rows", {
  d <- data.frame(x = 1:30, cl = c(rep("a", 10), rep("b", 10), rep(NA, 10)))

  # `==` against an NA label yields NA, and indexing with NA manufactures
  # all-NA rows that were never in the data.
  expect_error(draw(d, design_multistage("cl", n_clusters = 1, n = 5), seed = 1),
               "missing cluster label")

  res <- draw(d, design_multistage("cl", n_clusters = 1, n = 5, na_rm = TRUE),
              seed = 1)
  expect_false(anyNA(res$x))
  expect_equal(nrow(res), 5)

  # A guard reading the NA-inflated count would let this through, returning 15
  # rows from a 10-row cluster with replace = FALSE.
  expect_error(
    draw(d, design_multistage("cl", n_clusters = 1, n = 15, na_rm = TRUE),
         seed = 1),
    "exceeds the 10 row\\(s\\) available"
  )
})

test_that("an antimeridian-spanning region warns instead of failing mutely", {
  skip_if_not_installed("sf")
  df <- make_df()
  expect_warning(
    try(draw(df, design_spatial(c("lon", "lat"), antimeridian_poly(), n = 5),
             seed = 1), silent = TRUE),
    "degrees of longitude"
  )
  expect_silent(
    res <- draw(df, design_spatial(c("lon", "lat"), wide_poly(), n = 5), seed = 1)
  )
  expect_equal(nrow(res), 5)
})

test_that("the warning does not fire on legitimate multi-part regions", {
  skip_if_not_installed("sf")
  d <- data.frame(id = 1:3, lon = c(5, 25, 45), lat = c(5, 25, 45))
  expect_silent(draw(d, design_spatial(c("lon", "lat"), three_boxes(), n = 3),
                     seed = 1))
})

test_that("allocation hits the total and does not floor rare strata", {
  uneven <- data.frame(
    x = 1:143,
    g = rep(letters[1:7], times = c(31, 29, 23, 19, 17, 13, 11))
  )
  expect_equal(nrow(draw(uneven, design_stratified("g", n = 50), seed = 1)), 50)

  skewed <- data.frame(x = 1:1000, g = c(rep("big", 999), "rare"))
  expect_equal(nrow(draw(skewed, design_stratified("g", n = 10), seed = 1)), 10)
})

test_that("month buckets follow the calendar, not 30.44-day durations", {
  df <- data.frame(
    id = 1:365,
    ts = seq(as.POSIXct("2020-01-01", tz = "UTC"), by = "day", length.out = 365)
  )
  res <- draw(df, design_temporal("ts", "2020-01-01", "2020-04-01", interval = 1,
                                  per_interval = 1, unit = "months"), seed = 1)
  expect_equal(sort(as.integer(format(res$ts, "%m"))), 1:3)
})

test_that("count arguments are validated rather than passed to sample()", {
  expect_error(design_simple(n = 2.7), "single whole number")
  expect_error(design_systematic(interval = 0), "at least 1")
  expect_error(design_systematic(interval = 2.5), "single whole number")
})

test_that("drawing leaves the caller's RNG stream alone", {
  df <- make_df()
  set.seed(42)
  expected <- stats::runif(3)
  set.seed(42)
  invisible(draw(df, design_simple(n = 2), seed = 99))
  expect_equal(stats::runif(3), expected)
})

test_that("column order is preserved by every design", {
  df <- make_df()
  designs <- list(
    design_simple(n = 4),
    design_stratified("site", n = 4),
    design_systematic(interval = 10),
    design_cluster("site", n_clusters = 2),
    design_multistage("site", n_clusters = 2, n = 4),
    design_reservoir(n = 4),
    design_temporal("time", "2020-01-01 00:00:00", "2020-01-01 12:00:00",
                    interval = 6, per_interval = 1, unit = "hours")
  )
  for (d in designs) {
    expect_identical(names(draw(df, d, seed = 1)), names(df), info = class(d)[1])
  }
})

test_that("removed arguments stay removed", {
  # max_rows caused both halves of A2 and all of A8. normalization caused A5.
  # Reintroducing either should be a deliberate act, not an accident, so pin
  # them here: this is the real guard for the three cases whose failure mode is
  # now structurally impossible rather than merely fixed.
  expect_error(design_simple(n = 5, max_rows = 3), "unused argument")
  expect_error(design_stratified("g", n = 5, max_rows = 3), "unused argument")
  expect_error(design_cluster("g", n_clusters = 1, max_rows = 3),
               "unused argument")
  expect_error(design_multistage("g", n_clusters = 1, n = 5, max_rows = 3),
               "unused argument")
  expect_error(design_weighted("w", n = 5, max_rows = 3), "unused argument")
  expect_error(design_reservoir(n = 5, max_rows = 3), "unused argument")
  expect_error(design_systematic(interval = 5, max_rows = 3), "unused argument")
  expect_error(design_bootstrap(n_replicates = 5, max_rows = 3),
               "unused argument")
  expect_error(design_weighted("w", n = 5, normalization = "min-max"),
               "unused argument")
  expect_error(design_spatial(c("x", "y"), region = NULL, n = 5,
                              complex_region = TRUE), "unused argument")
})

test_that("weights = TRUE keeps the input's class, not just its columns", {
  skip_if_not_installed("tibble")
  # cbind() on a data.frame silently demotes a tibble, and the whole contract
  # is that what comes back matches what went in.
  tb <- tibble::tibble(id = 1:100, site = rep(c("a", "b"), each = 50))
  for (w in c(FALSE, TRUE)) {
    out <- draw(tb, design_stratified("site", n = 20), seed = 1, weights = w)
    expect_s3_class(out, "tbl_df")
  }
  wt <- draw(tb, design_simple(n = 10), seed = 1, weights = TRUE)
  expect_identical(names(wt), c(".prob", ".weight", "id", "site"))
})

test_that("simulation refuses the bootstrap instead of returning all ones", {
  d <- data.frame(id = 1:40, y = 1:40)
  des <- design_bootstrap(n_replicates = 20)
  expect_error(inclusion_prob(d, des, simulate = TRUE, R = 10),
               "not a probability sample")
  expect_error(joint_prob(d, des, rows = 1:5, simulate = TRUE, R = 10),
               "not a probability sample")
})

# ---- guards found by the audit round --------------------------------------

test_that("a design that samples with replacement refuses to weight", {
  df <- data.frame(id = 1:20, g = rep(c("a", "b"), each = 10), y = 1:20)
  for (des in list(design_simple(n = 8, replace = TRUE),
                   design_stratified("g", n = 8, replace = TRUE),
                   design_multistage("g", n_clusters = 1, n = 4,
                                     replace = TRUE))) {
    expect_error(draw(df, des, seed = 1, weights = TRUE), "with-replacement")
  }
})

test_that("duplicated column names are refused before they are silently renamed", {
  dd <- data.frame(a = 1:4, b = 5:8)
  names(dd) <- c("x", "x")
  expect_identical(names(draw(dd, design_simple(n = 2), seed = 1)),
                   c("x", "x"))
  expect_error(draw(dd, design_simple(n = 2), seed = 1, weights = TRUE),
               "duplicated column name")
})

test_that("bootstrap keeps its promises about class and column names", {
  skip_if_not_installed("tibble")
  expect_s3_class(draw(tibble::tibble(id = 1:9),
                       design_bootstrap(n_replicates = 2, n = 3), seed = 1),
                  "tbl_df")
  expect_error(draw(data.frame(id = 1:6, .replicate = 1:6),
                    design_bootstrap(n_replicates = 2, n = 3), seed = 1),
               "already has a column called `.replicate`")
  out <- draw(data.frame(id = 1:6), design_bootstrap(n_replicates = 2, n = 3),
              seed = 1)
  expect_identical(names(out), c(".replicate", "id"))
})

test_that("a certainty design refuses a `rest` that reshapes the frame", {
  expect_error(design_certainty("v", 50, design_bootstrap(n_replicates = 4)),
               "cannot be a bootstrap design")
})

test_that("counts are validated at construction, not left to become NA", {
  expect_error(design_simple(n = 3e9), "at most")
  expect_error(design_bootstrap(n_replicates = 3e9), "at most")
  expect_error(design_systematic(interval = 3e9), "at most")
})

test_that("a negative seed is a seed", {
  a <- draw(data.frame(id = 1:20), design_simple(n = 5), seed = -1)
  b <- draw(data.frame(id = 1:20), design_simple(n = 5), seed = -1)
  expect_identical(a, b)
  expect_error(draw(data.frame(id = 1:20), design_simple(n = 5), seed = 1.5),
               "whole number")
})

test_that("proportional allocation survives a frame big enough to overflow", {
  # `n * sizes` was evaluated in 32-bit integer arithmetic, so 3,000 rows out
  # of 2.4 million produced NAs and an unrelated error.
  big <- data.frame(g = rep(c("a", "b", "c"),
                            times = c(800000L, 900000L, 700000L)))
  out <- draw(big, design_stratified("g", n = 3000), seed = 1)
  expect_equal(nrow(out), 3000)
  expect_equal(sum(table(out$g)), 3000)
})

test_that("strata keys containing a dot stay distinct", {
  # interaction() joins levels with ".", so ("a.b", "c") and ("a", "b.c")
  # collapsed into one stratum -- wrong allocation, wrong probabilities, and
  # min_per_stratum silently unmet.
  d <- data.frame(id = 1:16,
                  g1 = c(rep("a.b", 6), rep("a", 2), rep("z", 8)),
                  g2 = c(rep("c", 6), rep("b.c", 2), rep("q", 8)))
  des <- design_stratified(c("g1", "g2"), n = 6, min_per_stratum = 1)
  expect_length(unique(inclusion_prob(d, des)), 3)
  for (s in 1:10) {
    got <- draw(d, des, seed = s)
    expect_true(any(got$g1 == "a" & got$g2 == "b.c"))
  }
})

test_that("designs validate what they read at construction where they can", {
  skip_if_not_installed("sf")
  expect_error(design_spatial(c("x", "y"), NULL, n = 5), "must be an sf")
  expect_error(design_spatial(c("x", "y"),
                              sf::st_sfc(sf::st_point(c(0, 0)), crs = 4326),
                              n = 5, crs = "nonsense"), "coordinate reference")
  expect_error(design_spatial(c("x", "y"), sf::st_sfc(sf::st_point(c(0, 0))),
                              n = 5), "no CRS")
  expect_error(design_bootstrap(10, method = "simple", block_length = 5),
               "only applies to method")
})

test_that("non-vector columns are caught wherever a design reads one", {
  d <- data.frame(id = 1:6)
  d$m <- matrix(1:12, nrow = 6)
  expect_error(draw(d, design_weighted("m", n = 2)), "flatten it first")
  expect_error(draw(d, design_certainty("m", 2, design_simple(n = 1))),
               "flatten it first")
})
