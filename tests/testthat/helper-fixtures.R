# Deterministic fixtures. Building a fixture from unseeded rnorm()/runif() means
# any future failure arrives with no way to reproduce it.

# Seeds the random number stream until the calling function or test ends, and
# then puts back the stream the caller had, or none if there was none. A bare
# set.seed() in a fixture or a test resets the stream every later test draws
# from, so a test that forgot its own seed would pass only after that one.
# Restores run last-in first-out, so several calls in one test unwind to the
# stream from before the first.
local_seed <- function(seed, envir = parent.frame()) {
  globals <- globalenv()
  had <- exists(".Random.seed", envir = globals, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globals, inherits = FALSE)
  restore <- function() {
    if (had) {
      assign(".Random.seed", old, envir = globals)
    } else if (exists(".Random.seed", envir = globals, inherits = FALSE)) {
      rm(".Random.seed", envir = globals)
    }
  }
  do.call(on.exit, list(as.call(list(restore)), add = TRUE, after = FALSE),
          envir = envir)
  set.seed(seed)
  invisible(seed)
}

make_df <- function(n = 100L) {
  local_seed(20250218)
  data.frame(
    id    = seq_len(n),
    value = stats::rnorm(n),
    site  = rep(letters[1:4], each = n / 4L),
    time  = format(
      rep(seq(as.POSIXct("2020-01-01 00:00:00", tz = "UTC"),
              by = "hour", length.out = n / 4L), 4L),
      "%Y-%m-%d %H:%M:%S"
    ),
    # Kept well inside a single hemisphere-width box: no polygon edge in these
    # tests needs to span 180 degrees of longitude, which is what collapses a
    # region under s2.
    lon   = stats::runif(n, -75, 75),
    lat   = stats::runif(n, -55, 55),
    stringsAsFactors = FALSE
  )
}

make_ts <- function() {
  data.frame(
    id = 1:48,
    ts = seq(as.POSIXct("2020-01-01", tz = "UTC"), by = "hour", length.out = 48)
  )
}

# A large but well-behaved region. Its edges are cut in degree space into
# pieces of at most 1 degree, so under s2 it is the box its corners describe:
# a single edge from -80 to 80 at latitude 60 would follow a great circle
# that bows about 24 degrees towards the pole.
wide_poly <- function() {
  skip_if_not_installed("sf")
  flat_box(c(-80, 80), c(-60, 60))
}

# A longitude/latitude box whose edges follow its straight lines in degrees,
# by the recipe in ?design_spatial.
flat_box <- function(lon, lat) {
  skip_if_not_installed("sf")
  box <- sf::st_sfc(sf::st_polygon(list(cbind(lon[c(1, 2, 2, 1, 1)],
                                              lat[c(1, 1, 2, 2, 1)]))))
  sf::st_set_crs(sf::st_segmentize(box, 1), 4326)
}

# A lon/lat rectangle whose edges span nearly the whole globe. Under s2 this
# collapses to a strip across the antimeridian; pinned so that stays detected.
antimeridian_poly <- function() {
  skip_if_not_installed("sf")
  sf::st_sfc(
    sf::st_polygon(list(cbind(
      c(-179, 179, 179, -179, -179),
      c(-89, -89, 89, 89, -89)
    ))),
    crs = 4326
  )
}

three_boxes <- function() {
  skip_if_not_installed("sf")
  lapply(list(c(0, 10), c(20, 30), c(40, 50)), function(r) flat_box(r, r))
}

expect_same_schema <- function(out, input) {
  expect_identical(names(out), names(input))
  expect_identical(class(out), class(input))
}

# Evaluates `code` under a collation that sorts "a" before "B", as most
# non-C locales do, and restores the session's collation afterwards. Skips
# where no such locale is installed. R reads the collation for text from the
# LC_ALL and LC_COLLATE environment variables before the locale, and testthat
# sets the second to "C", so both are set here too.
with_other_collation <- function(code) {
  old <- Sys.getlocale("LC_COLLATE")
  old_env <- Sys.getenv(c("LC_ALL", "LC_COLLATE"), unset = NA)
  on.exit({
    for (nm in names(old_env)) {
      if (is.na(old_env[[nm]])) {
        Sys.unsetenv(nm)
      } else {
        do.call(Sys.setenv, as.list(old_env[nm]))
      }
    }
    Sys.setlocale("LC_COLLATE", old)
  }, add = TRUE)
  Sys.unsetenv("LC_ALL")
  for (loc in c("en_US.UTF-8", "C.UTF-8", "C.utf8")) {
    Sys.setenv(LC_COLLATE = loc)
    if (nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", loc))) &&
        identical(sort(c("B", "a")), c("a", "B"))) {
      return(code)
    }
  }
  skip("No installed locale sorts \"a\" before \"B\".")
}

# A point frame for the tests that pass an sf object as `data`.
sf_points <- function() {
  skip_if_not_installed("sf")
  d <- data.frame(
    id = 1:60,
    site = rep(c("a", "b", "c"), times = c(30, 20, 10)),
    cl = rep(sprintf("c%02d", 1:12), each = 5),
    x = (1:60 %% 7) / 7,
    yy = (1:60 %% 11) / 11,
    y = 10 + (1:60 %% 13)
  )
  d$lon <- -60 + 2 * d$id
  d$lat <- -30 + d$id
  sf::st_as_sf(d, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
}

# ---- the design grid ---------------------------------------------------------

# A frame every design in design_rows() draws from: 240 rows in strata of
# 120, 80 and 40, 24 clusters of 10, a size `w` with a long right tail (about
# ten rows at 100 or more), hourly times over ten days, points in a lon/lat
# box and on the unit square, and a `y` that rises with `w` and `x` and
# shares a part within each cluster.
grid_frame <- function() {
  local_seed(20251112)
  n <- 240
  d <- data.frame(
    id = seq_len(n),
    site = rep(c("a", "b", "c"), times = c(120, 80, 40)),
    cl = rep(sprintf("c%02d", 1:24), each = 10),
    w = round(stats::rgamma(n, 2, 0.05), 1) + 1,
    ts = seq(as.POSIXct("2024-01-01", tz = "UTC"), by = "hour",
             length.out = n),
    lon = stats::runif(n, -60, 60),
    lat = stats::runif(n, -40, 40),
    x = stats::runif(n),
    yy = stats::runif(n),
    stringsAsFactors = FALSE
  )
  d$y <- 2 * d$w + 50 * d$x + rep(stats::rnorm(24, 0, 15), each = 10) +
    stats::rnorm(n, 0, 10)
  d$flag <- d$y > stats::median(d$y)
  d$day <- as.Date(d$ts)
  d
}

# One entry per design and common variant, each with the frame it draws from
# and what test-grid.R checks it against:
#   n          the rows a draw returns, and `fixed` whether every draw does
#   expected   the sum of the inclusion probabilities, `n` unless given
#   psu        the frame column whose values are the sampling units, or NULL
#              for rows; `strata` the strata columns, or a function of the
#              frame giving each row's stratum, or NULL for one stratum
#   method     the variance method that variance = "auto" reports
#   jackknife  "same" where the jackknife is the analytic estimator, "own"
#              where it is an estimate of its own, else the note it declines
#              with
#   survey     "same" where survey gives the same standard error, "close"
#              where it approximates the same estimator another way, "own"
#              where it has another estimator
#   refuse     the paths the design refuses on purpose, with the start of
#              the message: "weights", "inclusion" (also sampling_weight()),
#              "joint", "probability" (that plot view) and "plot" (every view)
#   warns      the class of a warning the design gives on every draw
#   conditional  whether a weighted sample is weighted on the clusters
#              selected, which leaves its expected size unknown
design_rows <- function() {
  f <- grid_frame()
  row <- function(design, n, fixed = TRUE, expected = n, frame = f,
                  psu = NULL, strata = NULL, method = "analytic",
                  jackknife = "same", survey = "same", refuse = list(),
                  warns = NULL, conditional = FALSE) {
    list(design = design, frame = frame, n = n, fixed = fixed,
         expected = expected, psu = psu, strata = strata, method = method,
         jackknife = jackknife, survey = survey, refuse = refuse,
         warns = warns, conditional = conditional)
  }
  no_joint <- "has no closed-form joint inclusion probability"
  fixed_start <- "with a fixed `start` is not a probability sample"
  resampling <- "it is not a probability sample of a finite population"
  successive <- "has no closed-form inclusion probability"
  systematic_jk <- "A systematic sample has one primary sampling unit"
  poisson_jk <- "Poisson sampling selects rows independently"
  hours <- function(d) {
    start <- as.POSIXct("2024-01-01", tz = "UTC")
    floor(as.numeric(difftime(d$ts, start, units = "hours")) / 24)
  }
  days <- function(d) floor(as.numeric(d$day - as.Date("2024-01-01")) / 2)
  takes <- sum(f$w >= 100)

  rows <- list(
    simple = row(design_simple(n = 30), 30),
    simple_replace = row(
      design_simple(n = 30, replace = TRUE), 30,
      expected = 240 * (1 - (1 - 1 / 240)^30),
      refuse = list(weights = "not available for a with-replacement design",
                    joint = no_joint)
    ),
    stratified = row(design_stratified("site", n = 30), 30, strata = "site"),
    stratified_equal = row(
      design_stratified("site", n = 30, allocation = "equal"), 30,
      strata = "site"
    ),
    stratified_neyman = row(
      design_stratified("site", n = 30, allocation = "neyman",
                        allocation_by = "w"), 30, strata = "site"
    ),
    stratified_two = row(
      design_stratified(c("site", "flag"), n = 30, min_per_stratum = 2), 30,
      strata = c("site", "flag")
    ),
    systematic = row(design_systematic(interval = 8), 30,
                     method = "successive difference",
                     jackknife = systematic_jk, survey = "own"),
    systematic_order = row(design_systematic(interval = 8, order_by = "w"), 30,
                           method = "successive difference",
                           jackknife = systematic_jk, survey = "own"),
    systematic_start = row(
      design_systematic(interval = 8, start = 3), 30,
      refuse = list(weights = fixed_start, inclusion = fixed_start,
                    joint = fixed_start, probability = fixed_start)
    ),
    cluster = row(design_cluster("cl", n_clusters = 5), 50, psu = "cl"),
    cluster_balanced = row(
      design_cluster("cl", n_clusters = 5, balanced = TRUE), 50, psu = "cl",
      jackknife = "own", refuse = list(joint = no_joint), conditional = TRUE
    ),
    multistage = row(design_multistage("cl", n_clusters = 6, n = 24), 24,
                     psu = "cl", jackknife = "own"),
    multistage_prop = row(
      design_multistage("cl", n_clusters = 6, n = 24,
                        allocation = "proportional", min_per_cluster = 1),
      24, psu = "cl", jackknife = "own", conditional = TRUE,
      refuse = list(inclusion = "has no closed-form inclusion probability",
                    joint = no_joint,
                    probability = "has no closed-form inclusion probability")
    ),
    weighted_successive = row(
      design_weighted("w", n = 30), 30,
      refuse = list(weights = successive, inclusion = successive,
                    joint = no_joint, probability = successive)
    ),
    weighted_systematic = row(
      design_weighted("w", n = 30, method = "systematic"), 30,
      method = "deville", jackknife = "own", survey = "close",
      refuse = list(joint = no_joint)
    ),
    weighted_poisson = row(
      design_weighted("w", n = 30, method = "poisson"), NA, fixed = FALSE,
      expected = 30, jackknife = poisson_jk
    ),
    certainty_simple = row(
      design_certainty("w", 100, design_simple(n = 20)), takes + 20
    ),
    certainty_stratified = row(
      design_certainty("w", 100, design_stratified("site", n = 20)),
      takes + 20, strata = "site"
    ),
    certainty_systematic = row(
      design_certainty("w", 100,
                       design_weighted("w", n = 20, method = "systematic")),
      takes + 20, method = "deville", jackknife = "own", survey = "close",
      refuse = list(joint = no_joint)
    ),
    certainty_poisson = row(
      design_certainty("w", 100,
                       design_weighted("w", n = 20, method = "poisson")),
      NA, fixed = FALSE, expected = takes + 20, jackknife = poisson_jk
    ),
    reservoir = row(design_reservoir(n = 30), 30),
    reservoir_capped = row(design_reservoir(n = 30, max_items = 150), 30,
                           warns = "drawn_stream_truncated"),
    reservoir_stream = row(
      design_reservoir(n = 30), 30, frame = as.list(f$id),
      refuse = list(weights = "needs a data frame, not a list",
                    inclusion = "must be a data frame, not a list",
                    joint = "must be a data frame, not a list",
                    plot = "must be a data frame, not a list")
    ),
    bootstrap = row(
      design_bootstrap(n_replicates = 20, n = 30), 600,
      refuse = list(weights = resampling, inclusion = resampling,
                    joint = "A bootstrap is not a probability sample",
                    probability = resampling)
    ),
    bootstrap_block = row(
      design_bootstrap(n_replicates = 20, n = 30, method = "block",
                       block_length = 5), 600,
      refuse = list(weights = resampling, inclusion = resampling,
                    joint = "A bootstrap is not a probability sample",
                    probability = resampling)
    ),
    temporal_hours = row(
      design_temporal("ts", from = "2024-01-01", to = "2024-01-11",
                      interval = 24, per_interval = 3, unit = "hours"),
      30, strata = hours
    ),
    temporal_days = row(
      design_temporal("day", from = "2024-01-01", to = "2024-01-11",
                      interval = 2, per_interval = 3, unit = "days"),
      15, strata = days
    ),
    spread = row(design_spread(c("x", "yy"), n = 30, scale = FALSE), 30,
                 method = "local mean", jackknife = "own", survey = "own",
                 refuse = list(joint = no_joint)),
    spread_size = row(
      design_spread(c("x", "yy"), n = 30, size = "w", scale = FALSE), 30,
      method = "local mean", jackknife = "own", survey = "own",
      refuse = list(joint = no_joint)
    )
  )
  if (requireNamespace("sf", quietly = TRUE)) {
    rows$spatial <- row(
      design_spatial(c("lon", "lat"), flat_box(c(-50, 50), c(-30, 30)),
                     n = 30), 30
    )
  }
  rows
}
