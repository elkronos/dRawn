# Deterministic fixtures. Building a fixture from unseeded rnorm()/runif() means
# any future failure arrives with no way to reproduce it.

make_df <- function(n = 100L) {
  set.seed(20250218)
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
