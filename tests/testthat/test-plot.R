test_that("both plot types render without error", {
  pop <- data.frame(id = 1:200, site = rep(c("a", "b"), times = c(150, 50)),
                    cl = rep(paste0("c", 1:20), each = 10), w = 1:200)
  f <- tempfile(fileext = ".png")
  grDevices::png(f); on.exit({grDevices::dev.off(); unlink(f)}, add = TRUE)

  expect_invisible(plot(design_simple(n = 20), pop, seed = 1))
  expect_silent(plot(design_systematic(interval = 5), pop, seed = 1))
  expect_silent(plot(design_cluster("cl", n_clusters = 3), pop, seed = 1))
  expect_silent(plot(design_stratified("site", n = 20), pop, type = "probability"))
  expect_silent(plot(design_weighted("w", n = 20, method = "systematic"), pop,
                     type = "probability"))
})

test_that("plot validates its input", {
  pop <- data.frame(id = 1:50)
  f <- tempfile(fileext = ".png")
  grDevices::png(f); on.exit({grDevices::dev.off(); unlink(f)}, add = TRUE)

  expect_error(plot(design_simple(n = 5)), "second argument")
  expect_error(plot(design_simple(n = 5), 1:10), "must be a data frame")
  expect_error(plot(design_simple(n = 5), pop, ncol = 0), "at least 1")
  # A design with no closed-form probability cannot be plotted that way
  expect_error(plot(design_weighted("id", n = 5), pop, type = "probability"),
               "no closed-form inclusion probability")
})

test_that("plot does not disturb the caller's graphics layout", {
  pop <- data.frame(id = 1:100)
  f <- tempfile(fileext = ".png")
  grDevices::png(f); on.exit({grDevices::dev.off(); unlink(f)}, add = TRUE)
  op <- graphics::par(mfrow = c(2, 2))
  # Restore before the device closes, or restoring opens a fresh Rplots.pdf
  on.exit(graphics::par(op), add = TRUE, after = FALSE)
  before <- graphics::par("mfrow")
  plot(design_simple(n = 10), pop, seed = 1)
  expect_equal(graphics::par("mfrow"), before)
})

test_that("large frames are subsampled for display", {
  pop <- data.frame(id = 1:20000)
  f <- tempfile(fileext = ".png")
  grDevices::png(f); on.exit({grDevices::dev.off(); unlink(f)}, add = TRUE)
  expect_silent(plot(design_simple(n = 500), pop, seed = 1, max_dots = 1000))
})

test_that("the default title describes the design", {
  expect_match(drawn:::design_label(design_stratified("site", n = 60)),
               "^design_stratified\\(strata = \"site\", n = 60")
  expect_match(drawn:::design_label(design_cluster("cl", n_clusters = 4)),
               "n_clusters = 4")
})

test_that("plot leaves graphical parameters as it found them", {
  pop <- data.frame(id = 1:200, site = rep(c("a", "b"), each = 100),
                    w = stats::runif(200, 1, 5))
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f), add = TRUE)
  grDevices::png(f)
  on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)

  before <- graphics::par(c("mar", "bg"))
  plot(design_simple(n = 20), pop)
  expect_identical(graphics::par(c("mar", "bg")), before)
  plot(design_stratified("site", n = 20), pop, type = "probability")
  expect_identical(graphics::par(c("mar", "bg")), before)
})

test_that("plot reports the real reason it cannot draw probabilities", {
  pop <- data.frame(id = 1:50, site = rep(c("a", "b"), each = 25))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(plot(design_stratified("nope", n = 10), pop,
                    type = "probability"),
               "missing column")
  expect_error(plot(design_systematic(interval = 7, start = 3), pop,
                    type = "probability"),
               "fixed `start`")
  pop$w <- 1:50
  expect_error(plot(design_weighted("w", n = 10), pop,
                    type = "probability"),
               "no closed-form")
  # Weights draw() refuses are refused first, for the same reason.
  expect_error(plot(design_weighted("site", n = 10), pop,
                    type = "probability"),
               "`site` must be numeric, not character", fixed = TRUE)
})

test_that("a spatial map draws the region's outline as the design reads it", {
  skip_if_not_installed("sf")
  skip_if_not(sf::sf_use_s2(), "spherical geometry is switched off")
  box <- sf::st_sfc(sf::st_polygon(list(cbind(c(-50, 50, 50, -50, -50),
                                             c(-30, -30, 30, 30, -30)))),
                    crs = 4326)
  d <- data.frame(id = 1:200, lon = (1:200 * 37) %% 120 - 60,
                  lat = (1:200 * 53) %% 90 - 45)
  des <- design_spatial(c("lon", "lat"), box, n = 20)

  # Under s2 the top edge follows its great circle up to latitude 41.9.
  outline <- drawn:::region_outline(des, d)
  expect_gt(max(outline$y, na.rm = TRUE), 41.9)
  expect_lt(max(outline$y, na.rm = TRUE), 42)
  planar <- function() {
    old <- suppressMessages(sf::sf_use_s2(FALSE))
    on.exit(suppressMessages(sf::sf_use_s2(old)))
    drawn:::region_outline(des, d)
  }
  flat <- planar()
  expect_equal(max(flat$y), 30)
  expect_length(flat$x, 5)

  drawn_calls <- function(code) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    grDevices::dev.control("enable")
    code
    vapply(grDevices::recordPlot()[[1]], function(e) {
      f <- e[[2]][[1]]
      if (is.list(f)) f$name else as.character(f)[1]
    }, character(1))
  }
  calls <- drawn_calls(expect_warning(plot(des, d, type = "map", seed = 1),
                                      class = "drawn_geodesic_edges"))
  expect_true("C_polygon" %in% calls)
  calls <- drawn_calls(plot(design_simple(n = 20), d, type = "map",
                            coords = c("lon", "lat"), seed = 1))
  expect_false("C_polygon" %in% calls)
})

test_that("a design with reserves is plotted by its main sample", {
  pop <- data.frame(id = 1:100, x = (1:100) %% 10, y = (1:100) %/% 10)
  f <- tempfile(fileext = ".png")
  grDevices::png(f); on.exit({grDevices::dev.off(); unlink(f)}, add = TRUE)
  shown <- NULL
  local_mocked_bindings(draw_titles = function(main, sub, pal) {
    shown <<- c(main, sub)
  })

  des <- design_simple(n = 10, reserve = 5)
  expect_identical(plot(des, pop, seed = 1), des)
  expect_equal(shown, c("design_simple(n = 10, reserve = 5)",
                        "10 of 100 rows selected"))
  spread <- design_spread(c("x", "y"), n = 10, reserve = 5)
  expect_identical(plot(spread, pop, type = "map", seed = 1), spread)
  expect_match(shown[2], "^10 of 100")
})
