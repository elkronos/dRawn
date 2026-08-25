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
  on.exit(graphics::par(op), add = TRUE)
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
               "design_stratified\\(n = 60")
  expect_match(drawn:::design_label(design_cluster("cl", n_clusters = 4)),
               "n_clusters = 4")
})
