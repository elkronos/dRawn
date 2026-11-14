test_that("samples rows inside the region", {
  skip_if_not_installed("sf")
  df <- make_df()
  res <- draw(df, design_spatial(c("lon", "lat"), wide_poly(), n = 5), seed = 123)
  expect_equal(nrow(res), 5)
  expect_same_schema(res, df)
})

test_that("coords must name two columns or one POINT geometry column", {
  expect_error(design_spatial(c("a", "b", "c"), region = NULL, n = 5),
               "exactly two columns")
  # One name is a geometry column or nothing, which only the frame can say.
  skip_if_not_installed("sf")
  df <- make_df()
  expect_error(draw(df, design_spatial("lon", wide_poly(), n = 5)),
               "two coordinate columns, c(x, y) -- longitude then latitude -- or one POINT geometry column, not `lon`.",
               fixed = TRUE)
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
  local_seed(7)
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
  bowing <- sf::st_sfc(sf::st_polygon(list(cbind(c(-80, 80, 80, -80, -80),
                                                 c(-60, -60, 60, 60, -60)))),
                       crs = 4326)
  grid <- expand.grid(lon = seq(-70, 70, by = 10),
                      lat = c(-50, -20, 0, 20, 50, 65, 70))
  grid$y <- seq_len(nrow(grid)) %% 7 + grid$lat / 10
  des <- design_spatial(c("lon", "lat"), bowing, n = 20)
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

  # A sample saved by drawn 0.2.0 has no `.row` and no `drawn_pi`, but its
  # `.prob` shows the same change.
  old <- s
  old$.row <- NULL
  attr(old, "drawn_pi") <- NULL
  expect_error(with_s2(TRUE, ht_total(old, "y")),
               "gives 20 of its 20 rows a different inclusion probability",
               fixed = TRUE)
  expect_error(with_s2(TRUE, sample_summary(old)),
               "in a session with the same version of sf", fixed = TRUE)
  expect_equal(with_s2(FALSE, ht_total(old, "y"))$se, planar$se)

  skip_if_not_installed("survey")
  expect_error(with_s2(TRUE, as_svydesign(s)), "sf_use_s2", fixed = TRUE)
  expect_error(with_s2(TRUE, as_svydesign(old)), "sf_use_s2", fixed = TRUE)
})

# ---- region checks, edges and coordinates -------------------------------------

lonlat_box <- function(lon, lat, crs = 4326) {
  sf::st_sfc(sf::st_polygon(list(cbind(lon[c(1, 2, 2, 1, 1)],
                                       lat[c(1, 1, 2, 2, 1)]))), crs = crs)
}

test_that("a list region in mixed CRSs is reprojected member by member", {
  skip_if_not_installed("sf")
  # c() refuses geometries in different CRSs, which came back as "`region`
  # must be an sf or sfc geometry, or a list of them, not list".
  df <- data.frame(id = 1:40, lon = seq(1, 39, length.out = 40), lat = 5)
  a <- lonlat_box(c(0, 10), c(0, 10))
  b <- lonlat_box(c(20, 30), c(0, 10))
  mixed <- design_spatial(c("lon", "lat"), list(a, sf::st_transform(b, 3857)),
                          n = 4)
  same <- design_spatial(c("lon", "lat"), list(a, b), n = 4)
  expect_equal(inclusion_prob(df, mixed), inclusion_prob(df, same))
  expect_equal(sum(inclusion_prob(df, mixed) > 0), 20L)
  expect_error(design_spatial(c("lon", "lat"),
                              list(a, sf::st_set_crs(b, NA)), n = 4),
               "`region[[2]]` has no CRS", fixed = TRUE)
  expect_error(design_spatial(c("lon", "lat"), list(a, "b"), n = 4),
               "`region[[2]]` must be an sf or sfc geometry, not character.",
               fixed = TRUE)
  expect_error(design_spatial(c("lon", "lat"), list(), n = 4),
               "`region` is an empty list.", fixed = TRUE)
})

test_that("a region that encloses no area is refused", {
  skip_if_not_installed("sf")
  line <- sf::st_sfc(sf::st_linestring(rbind(c(0, 0), c(10, 10))), crs = 4326)
  expect_error(design_spatial(c("lon", "lat"), line, n = 1),
               "`region` must be POLYGON or MULTIPOLYGON geometry, not LINESTRING.",
               fixed = TRUE)
  expect_error(design_spatial(c("lon", "lat"),
                              sf::st_sfc(sf::st_point(c(0, 0)), crs = 4326),
                              n = 1),
               "not POINT. A point or a line encloses no area", fixed = TRUE)
  coll <- sf::st_sfc(sf::st_geometrycollection(list(
    sf::st_polygon(list(cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)))))),
    crs = 4326)
  expect_error(design_spatial(c("lon", "lat"), coll, n = 1),
               "sf::st_collection_extract", fixed = TRUE)
  expect_error(design_spatial(c("lon", "lat"),
                              list(lonlat_box(c(0, 1), c(0, 1)), line), n = 1),
               "`region[[2]]` must be POLYGON or MULTIPOLYGON", fixed = TRUE)
})

test_that("coordinates outside longitude/latitude range are refused", {
  skip_if_not_installed("sf")
  box <- lonlat_box(c(-102, -99), c(39, 42))
  # Latitude and longitude the wrong way round
  swapped <- data.frame(lat = c(40, 41), lon = c(-100, -101))
  des <- design_spatial(c("lat", "lon"), box, n = 1)
  expect_error(draw(swapped, des, seed = 1),
               "`lon` holds -101, outside -90 to 90, the range of a latitude",
               fixed = TRUE)
  expect_error(inclusion_prob(swapped, des), "(lat, lon) order, swap them",
               fixed = TRUE)
  expect_error(draw(data.frame(lon = c(10, 20), lat = c(95, 5)),
                    design_spatial(c("lon", "lat"), box, n = 1)),
               "`lat` holds 95", fixed = TRUE)
  # Projected coordinates under the default crs
  expect_error(inclusion_prob(data.frame(lon = 512000, lat = 40),
                              design_spatial(c("lon", "lat"), box, n = 1)),
               "`lon` holds 512,000, outside -180 to 360", fixed = TRUE)
  # Longitudes from 0 to 360 are longitudes, and s2 wraps them
  east <- data.frame(lon = c(259, 200), lat = c(40, 40))
  expect_equal(inclusion_prob(east, design_spatial(c("lon", "lat"), box,
                                                   n = 1)),
               c(1, 0))
  # A projected crs is not range-checked
  proj <- design_spatial(c("x", "y"), lonlat_box(c(0, 1e6), c(0, 1e6), 3857),
                         n = 1, crs = 3857)
  expect_equal(inclusion_prob(data.frame(x = 5e5, y = 5e5), proj), 1)
})

test_that("an edge that bows far from its straight line warns on every path", {
  skip_if_not_installed("sf")
  # The top edge, at latitude 30 from -50 to 50, follows a great circle that
  # reaches latitude 41.9 at longitude 0, so rows up to there are inside.
  box <- lonlat_box(c(-50, 50), c(-30, 30))
  pts <- data.frame(id = 1:4, lon = 0, lat = c(29, 35, 41, 43))
  des <- design_spatial(c("lon", "lat"), box, n = 2)
  expect_warning(p <- inclusion_prob(pts, des),
                 class = "drawn_geodesic_edges")
  expect_equal(p, c(2 / 3, 2 / 3, 2 / 3, 0))
  expect_warning(inclusion_prob(pts, des), "passes 11.9 degrees",
                 fixed = TRUE)
  expect_warning(joint_prob(pts, des), class = "drawn_geodesic_edges")
  expect_warning(draw(pts, des, seed = 1), class = "drawn_geodesic_edges")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_warning(plot(des, pts, type = "probability"),
                 class = "drawn_geodesic_edges")
  expect_warning(plot(des, pts, type = "map", seed = 1),
                 class = "drawn_geodesic_edges")
  # Once per call, not once per simulated draw
  caught <- 0L
  withCallingHandlers(
    inclusion_prob(pts, des, simulate = TRUE, R = 20, seed = 1),
    drawn_geodesic_edges = function(w) {
      caught <<- caught + 1L
      invokeRestart("muffleWarning")
    })
  expect_identical(caught, 1L)

  # The documented recipe: vertices every 0.1 degrees in degree space
  flat <- sf::st_set_crs(sf::st_segmentize(sf::st_set_crs(box, NA), 0.1), 4326)
  straight <- design_spatial(c("lon", "lat"), flat, n = 1)
  expect_silent(p <- inclusion_prob(pts, straight))
  expect_equal(p, c(1, 0, 0, 0))

  # Planar geometry draws the edges straight
  old <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  expect_silent(p <- suppressMessages(
    inclusion_prob(pts, design_spatial(c("lon", "lat"), box, n = 1))))
  expect_equal(p, c(1, 0, 0, 0))
})

test_that("the antimeridian warning gives both readings of the edge", {
  skip_if_not_installed("sf")
  # A box meant to cross the antimeridian is what s2 draws, so the warning
  # says which reading it took instead of calling the region too small.
  dateline <- lonlat_box(c(179, -179), c(0, 2))
  d <- data.frame(lon = c(179.5, -179.5, 0), lat = 1)
  des <- design_spatial(c("lon", "lat"), dateline, n = 2)
  w <- tryCatch(inclusion_prob(d, des), warning = identity)
  expect_s3_class(w, "drawn_antimeridian")
  expect_match(conditionMessage(w), "joins longitudes 179 and -179",
               fixed = TRUE)
  expect_match(conditionMessage(w), "a span of 2 degrees", fixed = TRUE)
  expect_match(conditionMessage(w), "span the 358 degrees of longitude",
               fixed = TRUE)
  expect_equal(suppressWarnings(inclusion_prob(d, des)), c(1, 1, 0))
  expect_warning(joint_prob(d, des), class = "drawn_antimeridian")
})

test_that("a region larger than a hemisphere works split into a list", {
  skip_if_not_installed("sf")
  # Unioned, the two halves came back as the rest of the globe under s2.
  halves <- list(flat_box(c(-179, 0), c(-60, 60)),
                 flat_box(c(0, 179), c(-60, 60)))
  d <- data.frame(lon = c(0, 100, -100, 179.5, 0),
                  lat = c(0, 10, -10, 0, 70))
  des <- design_spatial(c("lon", "lat"), halves, n = 1)
  expect_silent(p <- inclusion_prob(d, des))
  expect_equal(p > 0, c(TRUE, TRUE, TRUE, FALSE, FALSE))
})

test_that("a POINT geometry column can stand in for the coordinates", {
  skip_if_not_installed("sf")
  pts <- sf_points()
  plain <- sf::st_drop_geometry(pts)
  # No row lies on an edge, where planar and spherical tests can differ.
  box <- flat_box(c(-49, 41), c(-24.5, 20.5))
  by_cols <- design_spatial(c("lon", "lat"), box, n = 8)
  by_geom <- design_spatial("geometry", box, n = 8)
  expect_equal(inclusion_prob(pts, by_geom), inclusion_prob(plain, by_cols))
  s <- draw(pts, by_geom, seed = 3, weights = TRUE)
  expect_identical(s$id, draw(plain, by_cols, seed = 3)$id)
  expect_equal(ht_total(s, "y")$se,
               ht_total(draw(plain, by_cols, seed = 3, weights = TRUE),
                        "y")$se)
  # The column's own CRS is used, and the region is reprojected to it
  expect_equal(inclusion_prob(sf::st_transform(pts, 3857), by_geom),
               inclusion_prob(pts, by_geom))
  # An empty point is a missing coordinate
  gap <- pts
  sf::st_geometry(gap)[10] <- sf::st_point()
  expect_error(inclusion_prob(gap, by_geom),
               "1 row(s) have missing coordinates in `geometry`", fixed = TRUE)
  p <- inclusion_prob(gap, design_spatial("geometry", box, n = 8,
                                          na_rm = TRUE))
  expect_equal(p[10], 0)
  expect_equal(sum(p > 0), 44L)
  lines <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(
    sf::st_linestring(rbind(c(0, 0), c(1, 1))),
    sf::st_linestring(rbind(c(0, 1), c(1, 2))), crs = 4326))
  expect_error(draw(lines, design_spatial("geometry", box, n = 1)),
               "holding LINESTRING geometry", fixed = TRUE)
})
