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
