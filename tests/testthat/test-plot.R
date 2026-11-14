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
  local_seed(54)
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

  # Both views of every output mark the rows of the first replicate. The
  # map drew the count matrix of output = "weights" as row numbers.
  located <- data.frame(id = 1:30, x = 1:30 %% 7, y = 1:30 %/% 7)
  boot <- draw(located, design_bootstrap(n_replicates = 20, n = 3), seed = 1)
  real <- plot_counts
  marked <- NULL
  local_mocked_bindings(plot_counts = function(x, y, seed) {
    out <- real(x, y, seed)
    marked <<- which(out > 0)
    out
  })
  for (output in c("rows", "indices", "weights")) {
    for (type in c("selection", "map")) {
      marked <- NULL
      plot(design_bootstrap(n_replicates = 20, n = 3, output = output),
           located, type = type, coords = c("x", "y"), seed = 1)
      expect_equal(marked, sort(unique(boot$id[boot$.replicate == 1])),
                   info = paste(output, type))
    }
  }

  # The probability view has nothing to show and does not point to the
  # selection view either.
  err <- tryCatch(plot(des, pop, type = "probability"), error = identity)
  expect_match(conditionMessage(err), "has no inclusion probability")
  expect_no_match(conditionMessage(err), "type = \"selection\"", fixed = TRUE)
})

test_that("palette, main and sub are checked before anything is drawn", {
  pop <- data.frame(id = 1:50)
  d <- design_simple(n = 5)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
  on.exit(graphics::par(op), add = TRUE, after = FALSE)
  before <- graphics::par(c("mar", "bg"))

  expect_error(plot(d, pop, palette = list(surface = "offwhite")),
               "`palette$surface` must be one colour", fixed = TRUE)
  expect_identical(graphics::par(c("mar", "bg")), before)
  expect_error(plot(d, pop, palette = list(acent = "red")),
               "`palette` has an entry `acent`, which the plot does not use",
               fixed = TRUE)
  expect_error(plot(d, pop, palette = list(accent = NULL)),
               "`palette$accent` must be one colour, a name from colors() or a code such as \"#1B6B52\", not NULL.",
               fixed = TRUE)
  expect_error(plot(d, pop, palette = list(accent = c("red", "blue"))),
               "`palette$accent` must be one colour", fixed = TRUE)
  expect_error(plot(d, pop, palette = "red"), "named list of colours")
  expect_error(plot(d, pop, main = c("a", "b")),
               "`main` must be a single string, an expression or NULL")
  expect_error(plot(d, pop, sub = 1), "`sub` must be a single string or NULL")
  expect_error(plot(d, pop, palette = c(accent = "red", "blue")),
               "`palette` entry 2 has no name. Name each entry after the colour it replaces",
               fixed = TRUE)
  # A function has no value to show. is.na() of it warned.
  expect_no_warning(
    expect_error(plot(d, pop, main = mean),
                 "`main` must be a single string, an expression or NULL, not a function.",
                 fixed = TRUE)
  )
  expect_no_warning(
    expect_error(plot(d, pop, palette = list(accent = mean)),
                 "`palette$accent` must be one colour, a name from colors() or a code such as \"#1B6B52\", not a function.",
                 fixed = TRUE)
  )
  expect_identical(graphics::par(c("mar", "bg")), before)

  # A name is plotmath, as a call is.
  expect_no_condition(plot(d, pop, seed = 1, main = quote(alpha)))

  # A named character vector is a palette too, and an expression a title.
  shown <- NULL
  local_mocked_bindings(draw_titles = function(main, sub, pal) {
    shown <<- list(main = main, accent = pal$accent)
  })
  plot(d, pop, seed = 1, palette = c(accent = "red"),
       main = expression(alpha))
  expect_identical(shown, list(main = expression(alpha), accent = "red"))
  plot(d, pop, seed = 1, main = quote(alpha))
  expect_identical(shown$main, quote(alpha))
})

test_that("an error while drawing leaves the margins as they were", {
  pop <- data.frame(id = 1:50, x = 1:50, y = 50:1)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
  on.exit(graphics::par(op), add = TRUE, after = FALSE)
  before <- graphics::par(c("mar", "bg"))
  local_mocked_bindings(draw_titles = function(main, sub, pal) stop("boom"))
  for (type in c("selection", "probability", "map")) {
    expect_error(plot(design_simple(n = 5), pop, type = type,
                      coords = c("x", "y"), seed = 1), "boom")
    expect_identical(graphics::par(c("mar", "bg")), before, info = type)
  }
  # A background par() refuses after it has set the margins. The margins
  # were left set for the next plot.
  local_mocked_bindings(check_palette = function(palette) {
    utils::modifyList(drawn:::drawn_pal(), list(surface = "offwhite"))
  })
  for (type in c("selection", "probability", "map")) {
    expect_error(plot(design_simple(n = 5), pop, type = type,
                      coords = c("x", "y"), seed = 1), "offwhite")
    expect_identical(graphics::par(c("mar", "bg")), before, info = type)
  }
})

test_that("`sub` is added to the line under the title in every view", {
  pop <- data.frame(id = 1:100, x = (1:100) %% 10, y = (1:100) %/% 10)
  d <- design_simple(n = 10)
  shown <- NULL
  local_mocked_bindings(draw_titles = function(main, sub, pal) shown <<- sub)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  plot(d, pop, seed = 1, sub = "pilot frame, May")
  expect_equal(shown, "10 of 100 rows selected   |   pilot frame, May")
  plot(d, pop, type = "probability", sub = "pilot frame, May")
  expect_equal(shown, "range 0.1 to 0.1   |   pilot frame, May")
  # The map stopped with 'graphical parameter "bg" has the wrong length'.
  plot(d, pop, type = "map", coords = c("x", "y"), seed = 1,
       sub = "pilot frame, May")
  expect_equal(shown, "10 of 100 rows selected   |   pilot frame, May")
})

test_that("graphical parameters in `...` replace the plot's own", {
  pop <- data.frame(id = 1:100, x = (1:100 %% 10) * 1000, y = (1:100) %/% 10)
  window <- function(calls) calls[["C_plot_window"]]

  # Each of these stopped with 'matched by multiple actual arguments'.
  calls <- recorded_calls(plot(design_simple(n = 10), pop,
                               type = "probability", ylim = c(0, 1),
                               xlim = c(1, 50)))
  expect_equal(window(calls)[[1]], c(1, 50))
  expect_equal(window(calls)[[2]], c(0, 1))
  # The rules and the axis follow the limits given.
  axes <- calls[names(calls) == "C_axis"]
  side_2 <- axes[vapply(axes, function(a) a[[1]] == 2, logical(1))][[1]]
  expect_equal(side_2[[2]], seq(0, 1, by = 0.2))
  expect_equal(calls[["C_segments"]][[1]], 1)
  expect_equal(calls[["C_segments"]][[3]], 50)

  calls <- recorded_calls(plot(design_simple(n = 10), pop, seed = 1,
                               xlim = c(0, 30)))
  expect_equal(window(calls)[[1]], c(0, 30))

  calls <- recorded_calls(plot(design_simple(n = 10), pop, type = "map",
                               coords = c("x", "y"), seed = 1, asp = 1,
                               xlim = c(0, 5000)))
  expect_equal(window(calls)[[1]], c(0, 5000))
  expect_equal(window(calls)[[4]], 1)

  # Coordinates a design measures distance in get one scale on both axes;
  # other columns, and scaled ones, fill the panel.
  asp_of <- function(design, ...) {
    window(recorded_calls(plot(design, pop, type = "map", seed = 1, ...)))[[4]]
  }
  expect_equal(asp_of(design_spread(c("x", "y"), n = 10, scale = FALSE)), 1)
  expect_true(is.na(asp_of(design_spread(c("x", "y"), n = 10))))
  expect_true(is.na(asp_of(design_simple(n = 10), coords = c("x", "y"))))
  expect_true(is.na(asp_of(design_spread(c("x", "y"), n = 10, scale = FALSE),
                           asp = NA)))
})

test_that("a longitude/latitude map is scaled by its middle latitude", {
  skip_if_not_installed("sf")
  pop <- data.frame(id = 1:100, lon = (1:100 %% 10) / 10,
                    lat = 59.5 + (1:100 %/% 10) / 10)
  box <- sf::st_sfc(sf::st_polygon(list(cbind(c(0, 1, 1, 0, 0),
                                             c(59, 59, 61, 61, 59)))),
                    crs = 4326)
  des <- design_spatial(c("lon", "lat"), box, n = 10)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  asp <- NULL
  local_mocked_bindings(panel_frame = function(defaults, dots) {
    asp <<- defaults$asp
    do.call(graphics::plot, defaults)
  })
  suppressWarnings(plot(des, pop, type = "map", seed = 1))
  expect_equal(asp, 1 / cos(mean(range(pop$lat)) * pi / 180))
})

test_that("`...` takes graphical parameters and nothing else", {
  pop <- data.frame(id = 1:50)
  d <- design_simple(n = 5)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(plot(d, pop, simulate = TRUE),
               "plot() for a design does not take `simulate`.", fixed = TRUE)
  expect_error(plot(d, pop, simulate = TRUE, R = 10),
               "does not take `simulate` or `R`.", fixed = TRUE)
  expect_error(plot(d, pop, panel.first = stop("evaluated")),
               "does not take `panel.first`", fixed = TRUE)
  expect_error(plot(d, pop, "selection", 1, NULL, 4000, NULL, NULL, NULL,
                    NULL, 5),
               "must be named graphical parameters")
  expect_silent(plot(d, pop, seed = 1, col = "red", cex = 2, las = 1))
})

test_that("the map names columns it leaves out and rows it cannot place", {
  pop <- data.frame(id = 1:240, x = (1:240 * 7) %% 31, yy = (1:240 * 11) %% 37,
                    z = (1:240 * 13) %% 41)
  shown <- NULL
  local_mocked_bindings(draw_titles = function(main, sub, pal) shown <<- sub)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  plot(design_spread(c("x", "yy", "z"), n = 30, scale = FALSE), pop,
       type = "map", seed = 5)
  expect_equal(shown,
               "30 of 240 rows selected   |   placed by x and yy, not z")

  d <- design_simple(n = 30)
  ids <- draw(pop, d, seed = 5)$id
  gaps <- pop
  gaps$x[ids[1:3]] <- NA
  gaps$yy[ids[4]] <- Inf
  plot(d, gaps, type = "map", coords = c("x", "yy"), seed = 5)
  expect_equal(shown, paste0("30 of 240 rows selected   |   4 selected rows ",
                             "not shown, missing x or yy"))

  # The background holds the 210 rows not selected, and 50 of them are drawn.
  plot(d, pop, type = "map", coords = c("x", "yy"), seed = 5, max_dots = 50)
  expect_equal(shown, paste0("30 of 240 rows selected   |   background ",
                             "thinned to 50 of 210"))

  gaps$x <- NA_real_
  expect_error(plot(d, gaps, type = "map", coords = c("x", "yy"), seed = 5),
               "`x` and `yy` are not both finite in any row of `y`",
               fixed = TRUE)
})

test_that("the default title gives the arguments that differ by design", {
  label <- drawn:::design_label
  expect_equal(label(design_stratified(c("site", "flag"), n = 30)),
               paste0("design_stratified(strata = c(\"site\", \"flag\"), ",
                      "n = 30, allocation = \"proportional\")"))
  expect_equal(label(design_simple(n = 10, replace = TRUE)),
               "design_simple(n = 10, replace = TRUE)")
  expect_equal(label(design_systematic(interval = 8, start = 3,
                                       order_by = "w")),
               "design_systematic(interval = 8, start = 3, order_by = \"w\")")
  expect_equal(label(design_reservoir(n = 30, max_items = 150)),
               "design_reservoir(n = 30, max_items = 150)")
  expect_equal(label(design_reservoir(n = 30)), "design_reservoir(n = 30)")
  expect_equal(label(design_certainty("w", 100, design_simple(n = 20))),
               paste0("design_certainty(above = \"w\", threshold = 100, ",
                      "rest = design_simple(n = 20))"))
  expect_equal(label(design_bootstrap(20, method = "block",
                                      block_length = 5)),
               paste0("design_bootstrap(n_replicates = 20, ",
                      "method = \"block\", block_length = 5)"))
  expect_equal(label(design_temporal("ts", "2024-01-01", "2024-01-02", 24,
                                     5)),
               paste0("design_temporal(time = \"ts\", interval = 24, ",
                      "unit = \"hours\", per_interval = 5)"))
  expect_equal(label(design_stratified("site", n = 30, allocation = "neyman",
                                       allocation_by = "w",
                                       min_per_stratum = 2)),
               paste0("design_stratified(strata = \"site\", n = 30, ",
                      "allocation = \"neyman\", allocation_by = \"w\", ",
                      "min_per_stratum = 2)"))
  expect_equal(label(design_cluster("cl", n_clusters = 4)),
               "design_cluster(clusters = \"cl\", n_clusters = 4)")
  expect_equal(label(design_spread(c("x", "y"), n = 10, scale = FALSE)),
               paste0("design_spread(across = c(\"x\", \"y\"), n = 10, ",
                      "scale = FALSE)"))

  # Without sf regions or functions, a label is the call that builds the
  # design.
  designs <- list(
    design_simple(n = 10, replace = TRUE),
    design_stratified(c("site", "flag"), n = 30, min_per_stratum = 1,
                      na_rm = TRUE),
    design_systematic(interval = 7, n_starts = 2),
    design_multistage("cl", n_clusters = 2, n = 6, min_per_cluster = 1),
    design_weighted("w", n = 5, method = "poisson"),
    design_certainty("w", 10, design_stratified("s", 5,
                                                allocation = c(a = 2L, b = 3L))),
    design_bootstrap(20, method = "circular", output = "indices"),
    design_reservoir(n = 3, max_items = 50, batches = TRUE)
  )
  for (des in designs) {
    expect_equal(eval(parse(text = label(des))), des, info = label(des))
  }
})

test_that("plot names `y` when it is not a data frame", {
  expect_error(plot(design_simple(n = 2), as.list(data.frame(id = 1:5))),
               "`y` must be a data frame, not a list.", fixed = TRUE)
  expect_error(plot(design_simple(n = 2), data.frame(id = integer())),
               "`y` has no rows.", fixed = TRUE)
})

test_that("a title too long for the figure is shrunk to fit", {
  pop <- data.frame(id = 1:50)
  title_of <- function(main) {
    calls <- recorded_calls(plot(design_simple(n = 5), pop, seed = 1,
                                 main = main))
    calls[names(calls) == "C_mtext"][[1]]
  }
  short <- title_of("short")
  expect_equal(short[[1]], "short")
  expect_equal(short[[8]], 0.9)
  main <- strrep("design_stratified(strata = \"band\") ", 3)
  long <- title_of(main)
  expect_equal(long[[1]], main)
  expect_lt(long[[8]], 0.9)
  expect_gt(long[[8]], 0.6)
  # Past two thirds of its size, the text is cut and says so.
  longest <- title_of(strrep("x", 1000))
  expect_equal(longest[[8]], 0.6)
  expect_match(longest[[1]], "^x+\\.\\.\\.$")
  expect_lt(nchar(longest[[1]]), 1000)
})

test_that("a frame with the column plot() tags rows with is refused", {
  pop <- data.frame(id = 1:20, x = 1:20, .drawn_plot_id = 1)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (type in c("selection", "map")) {
    expect_error(plot(design_simple(n = 5), pop, type = type, seed = 1,
                      coords = c("id", "x")),
                 "`y` already has a column called `.drawn_plot_id`.",
                 fixed = TRUE)
  }
  names(pop)[3] <- "tag"
  expect_invisible(plot(design_simple(n = 5), pop, seed = 1))
})
