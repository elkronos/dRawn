test_that("samples rows inside the region", {
  skip_if_not_installed("sf")
  df <- make_df()
  res <- draw(df, design_spatial(c("lon", "lat"), wide_poly(), n = 5), seed = 123)
  expect_equal(nrow(res), 5)
  expect_same_schema(res, df)
})

test_that("coords must name exactly two columns", {
  expect_error(design_spatial("lon", region = NULL, n = 5), "exactly two columns")
  expect_error(design_spatial(c("a", "b", "c"), region = NULL, n = 5),
               "exactly two columns")
})

test_that("a region containing no rows is reported clearly", {
  skip_if_not_installed("sf")
  df <- make_df()
  empty <- sf::st_sfc(
    sf::st_polygon(list(cbind(c(20, 21, 21, 20, 20), c(20, 20, 21, 21, 20)))),
    crs = 4326
  )
  expect_error(draw(df, design_spatial(c("lon", "lat"), empty, n = 1), seed = 1),
               "No rows fall inside")
})

test_that("a list of geometries is unioned, whatever its length", {
  skip_if_not_installed("sf")
  set.seed(7)
  df <- data.frame(
    id  = 1:60,
    lon = c(stats::runif(20, 1, 9), stats::runif(20, 21, 29),
            stats::runif(20, 41, 49))
  )
  df$lat <- df$lon

  res <- draw(df, design_spatial(c("lon", "lat"), three_boxes(), n = 30), seed = 1)
  expect_equal(nrow(res), 30)
  expect_equal(length(unique(cut(res$lat, c(0, 10, 30, 50)))), 3)
})

test_that("a region in another CRS is reprojected, and a CRS-less one errors", {
  skip_if_not_installed("sf")
  df <- data.frame(id = 1:20, lon = seq(1, 9, length.out = 20),
                   lat = seq(1, 9, length.out = 20))
  box <- sf::st_sfc(
    sf::st_polygon(list(cbind(c(0, 10, 10, 0, 0), c(0, 0, 10, 10, 0)))),
    crs = 4326
  )
  expect_equal(
    nrow(draw(df, design_spatial(c("lon", "lat"), sf::st_transform(box, 3857),
                                 n = 5), seed = 1)),
    5
  )
  expect_error(
    draw(df, design_spatial(c("lon", "lat"), sf::st_set_crs(box, NA), n = 5),
         seed = 1),
    "has no CRS"
  )
})

test_that("missing coordinates follow the na_rm contract", {
  skip_if_not_installed("sf")
  df <- make_df()
  df$lat[1] <- NA
  expect_error(draw(df, design_spatial(c("lon", "lat"), wide_poly(), n = 5),
                    seed = 1),
               "missing coordinates")
  expect_equal(
    nrow(draw(df, design_spatial(c("lon", "lat"), wide_poly(), n = 5,
                                 na_rm = TRUE), seed = 1)),
    5
  )
})

test_that("over-drawing the region is refused", {
  skip_if_not_installed("sf")
  df <- make_df()
  expect_error(draw(df, design_spatial(c("lon", "lat"), wide_poly(), n = 500),
                    seed = 1),
               "exceeds the")
})

test_that("a sample is refused once sf_use_s2() differs from the draw", {
  skip_if_not_installed("sf")
  # Under spherical geometry the region's top edge, from -80 to 80 degrees of
  # longitude at latitude 60, is a great circle that bows north past the rows
  # at latitude 65 and 70. Under planar geometry those rows are outside.
  grid <- expand.grid(lon = seq(-70, 70, by = 10),
                      lat = c(-50, -20, 0, 20, 50, 65, 70))
  grid$y <- seq_len(nrow(grid)) %% 7 + grid$lat / 10
  des <- design_spatial(c("lon", "lat"), wide_poly(), n = 20)
  with_s2 <- function(value, code) {
    old <- suppressMessages(sf::sf_use_s2(value))
    on.exit(suppressMessages(sf::sf_use_s2(old)))
    suppressMessages(code)
  }

  s <- with_s2(FALSE, draw(grid, des, seed = 1, weights = TRUE))
  expect_equal(unique(s$.prob), 20 / 75)
  planar <- with_s2(FALSE, ht_total(s, "y"))
  expect_error(with_s2(TRUE, ht_total(s, "y")),
               "gives 105 of its 105 frame rows a different inclusion probability",
               fixed = TRUE)
  expect_error(with_s2(TRUE, ht_mean(s, "y")),
               "sf::sf_use_s2() decides which rows of a spatial design fall",
               fixed = TRUE)
  expect_error(with_s2(TRUE, sample_summary(s)), "it is TRUE now", fixed = TRUE)
  # Back under the setting it was drawn with, nothing has changed.
  expect_equal(with_s2(FALSE, ht_total(s, "y"))$se, planar$se)

  skip_if_not_installed("survey")
  expect_error(with_s2(TRUE, as_svydesign(s)), "sf_use_s2", fixed = TRUE)
})
