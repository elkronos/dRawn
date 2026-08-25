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

# A large but well-behaved region: every edge spans 160 degrees of longitude or
# less, so s2 draws it the way the coordinates read.
wide_poly <- function() {
  skip_if_not_installed("sf")
  sf::st_sfc(
    sf::st_polygon(list(cbind(
      c(-80, 80, 80, -80, -80),
      c(-60, -60, 60, 60, -60)
    ))),
    crs = 4326
  )
}

# The literal polygon from the retired UAT, kept so its behaviour stays pinned.
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
  lapply(list(c(0, 10), c(20, 30), c(40, 50)), function(r) {
    sf::st_sfc(
      sf::st_polygon(list(cbind(c(r[1], r[2], r[2], r[1], r[1]),
                                c(r[1], r[1], r[2], r[2], r[1])))),
      crs = 4326
    )
  })
}

expect_same_schema <- function(out, input) {
  expect_identical(names(out), names(input))
  expect_identical(class(out), class(input))
}
