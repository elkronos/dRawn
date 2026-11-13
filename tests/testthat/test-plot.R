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

# The calls a plot records, with their arguments, on a null device.
recorded_calls <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  grDevices::dev.control("enable")
  code
  calls <- grDevices::recordPlot()[[1]]
  names(calls) <- vapply(calls, function(e) {
    f <- e[[2]][[1]]
    if (is.list(f)) f$name else as.character(f)[1]
  }, character(1))
  lapply(calls, function(e) e[[2]][-1])
}

test_that("the probability view counts and draws every row of a large frame", {
  n <- 10000
  pop <- data.frame(id = seq_len(n),
                    site = rep(c("a", "b", "c"), length.out = n),
                    w = rep(1:50, length.out = n))
  pop$site[c(2, 5000)] <- NA
  pop$w[c(3, 9999)] <- 5000
  des <- design_certainty("w", 1000,
                          design_stratified("site", n = 500, na_rm = TRUE))
  p <- inclusion_prob(pop, des)

  view <- drawn:::probability_view(p, max_dots = 4000)
  expect_equal(view$size, 3)
  expect_equal(c(view$n_zero, view$n_certain), c(2L, 2L))
  expect_equal(c(view$low, view$top, view$mean), c(0, 1, mean(p)))
  # Rows 2 and 5000 fall in different runs, so each gets its own mark, at
  # its frame row.
  expect_equal(view$zero, c(2, 5000))
  expect_equal(range(view$y), c(0, 1))
  expect_equal(range(view$x), c(0.5, n + 0.5))

  shown <- NULL
  local_mocked_bindings(draw_titles = function(main, sub, pal) shown <<- sub)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(des, pop, type = "probability")
  expect_equal(shown, paste("2 rows unreachable", "2 rows taken with certainty",
                            "range 0 to 1", sep = "   |   "))
})

test_that("a run of rows keeps the order of its lowest and highest values", {
  p <- seq(1, 0.01, length.out = 1000)
  view <- drawn:::probability_view(p, max_dots = 100)
  expect_equal(view$size, 10)
  # A falling frame draws a falling step, with no saw teeth.
  expect_true(all(diff(view$y) <= 0))
  expect_true(all(diff(view$x) >= 0))
})

test_that("the probability view gives every row a plateau of its own", {
  view <- drawn:::probability_view(c(1, 1, 1, 3) / 3, max_dots = 4000)
  expect_equal(view$x, c(0.5, 1.5, 1.5, 2.5, 2.5, 3.5, 3.5, 4.5))
  # The last row, taken whole, is as wide as the others.
  expect_equal(view$y[7:8], c(1, 1))

  expect_equal(drawn:::row_ticks(1), 1)
  expect_equal(drawn:::row_ticks(2), 1:2)
  expect_equal(drawn:::row_ticks(3), 1:3)
  expect_equal(drawn:::row_ticks(4), 1:4)
  expect_equal(drawn:::row_ticks(10000), seq(2000, 10000, by = 2000))

  # A one-row frame draws its step and its fill.
  calls <- recorded_calls(plot(design_simple(n = 1), data.frame(id = 1),
                               type = "probability"))
  expect_true("C_polygon" %in% names(calls))
  expect_equal(sum(names(calls) == "C_plotXY"), 2)
  calls <- recorded_calls(plot(design_simple(n = 1), data.frame(id = 1:2),
                               type = "probability"))
  axes <- calls[names(calls) == "C_axis"]
  side_1 <- axes[vapply(axes, function(a) a[[1]] == 1, logical(1))][[1]]
  expect_equal(side_1[[2]], 1:2)
  expect_equal(side_1[[3]], c("1", "2"))
  # The mean is labelled right-aligned above the end of its line, inside the
  # panel.
  label <- calls[["C_text"]]
  expect_equal(label[[1]]$x, 2.5)
  expect_equal(label[[3]], c(1, -0.4))
  expect_null(label[[4]])
})

test_that("the selection view counts every selected row of a large frame", {
  pop <- data.frame(id = 1:4000)
  shown <- NULL
  local_mocked_bindings(draw_titles = function(main, sub, pal) shown <<- sub)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  # An evenly spaced subset of 1,500 rows fell in step with every fourth row,
  # so some starts showed no selected row at all.
  filled <- integer(4)
  for (start in 1:4) {
    des <- design_systematic(interval = 4, start = start)
    plot(des, pop, max_dots = 1500)
    expect_equal(shown, "1,000 of 4,000 rows selected")
    counts <- tabulate(draw(pop, des)$id, 4000)
    view <- drawn:::selection_view(counts, 1500)
    expect_equal(view$size, 3)
    expect_match(view$key, "each dot is 3 rows")
    filled[start] <- sum(view$filled)
  }
  # Three of every four runs of three rows hold a selected row, whatever the
  # start.
  expect_true(all(abs(filled - 1334 * 3 / 4) <= 1))
})

test_that("a with-replacement selection counts draws and rings repeats", {
  pop <- data.frame(id = 1:240, x = (1:240) %% 16, y = (1:240) %/% 16)
  des <- design_simple(n = 60, replace = TRUE)
  ids <- draw(pop, des, seed = 3)$id
  counts <- tabulate(ids, 240)

  shown <- NULL
  local_mocked_bindings(draw_titles = function(main, sub, pal) shown <<- sub)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(des, pop, seed = 3)
  expect_equal(shown, paste0("60 draws, ", sum(counts > 0),
                             " distinct rows of 240"))
  view <- drawn:::selection_view(counts, 4000)
  expect_equal(which(view$ringed), which(counts > 1))
  expect_equal(view$key, "a ring marks a row drawn more than once")

  plot(des, pop, type = "map", coords = c("x", "y"), seed = 3)
  expect_equal(shown, paste0("60 draws, ", sum(counts > 0),
                             " distinct rows of 240   |   ringed rows were ",
                             "drawn more than once"))
})

test_that("a bootstrap is plotted by its first replicate", {
  pop <- data.frame(id = 1:240)
  des <- design_bootstrap(n_replicates = 50, n = 30)
  reps <- draw(pop, des, seed = 1)
  first <- reps$id[reps$.replicate == 1]

  shown <- NULL
  local_mocked_bindings(draw_titles = function(main, sub, pal) shown <<- sub)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(des, pop, seed = 1), des)
  expect_equal(shown, paste0("replicate 1 of 50   |   30 draws, ",
                             length(unique(first)), " distinct rows of 240"))
  # The matrix outputs give the same picture.
  plot(design_bootstrap(n_replicates = 50, n = 30, output = "weights"), pop,
       seed = 1)
  expect_equal(shown, paste0("replicate 1 of 50   |   30 draws, ",
                             length(unique(first)), " distinct rows of 240"))

  # The probability view has nothing to show and does not point to the
  # selection view either.
  err <- tryCatch(plot(des, pop, type = "probability"), error = identity)
  expect_match(conditionMessage(err), "has no inclusion probability")
  expect_no_match(conditionMessage(err), "type = \"selection\"", fixed = TRUE)
})
