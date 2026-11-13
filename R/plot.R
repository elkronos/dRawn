# Palette. One accent plus neutrals, because both views encode two states or one
# magnitude -- not a set of categories. Validated on a light surface: the accent
# and the recessive fill separate by dE 36 under protanopia and 39 under normal
# vision (targets >= 8 and >= 15), and accent, ink and secondary all clear WCAG
# AA text contrast against the surface. The recessive fill sits deliberately
# below 3:1 -- it carries no information on its own, and the count annotation
# supplies the relief.
drawn_pal <- function() {
  list(
    surface   = "#FCFCFB",
    ink       = "#14201C",
    secondary = "#5B6763",
    muted     = "#98A29E",
    recessive = "#D3DAD6",
    accent    = "#1B6B52",
    fill      = "#DCEAE3",
    rule      = "#E6EAE8"
  )
}

#' Picture a sampling design
#'
#' Draws the design against a population so you can see what it selects. Two
#' views, and each answers a question that a table answers slowly:
#'
#' \describe{
#'   \item{`"selection"`}{Every row of the frame as a dot, in frame order, with
#'     the selected ones filled in. Designs look distinct: simple random
#'     sampling scatters, systematic makes a lattice, cluster sampling takes
#'     solid contiguous runs, and probability-proportional-to-size thickens
#'     wherever the weight is large. If your frame is sorted by something
#'     meaningful, an unintended pattern shows up immediately. A design that
#'     samples with replacement rings each row it drew more than once, and the
#'     count under the title gives the draws and the distinct rows. A
#'     bootstrap is shown by its first replicate.}
#'   \item{`"probability"`}{Each row's inclusion probability across the frame,
#'     as a step, with frame rows along the axis. Flat means every row had the
#'     same chance; plateaus mean strata; a rise means size-proportional
#'     selection. Rows the design can never reach sit at zero and are marked,
#'     which is usually the thing worth finding out. The count of those rows,
#'     the count taken with certainty, the range and the mean cover every row
#'     of the frame.}
#'   \item{`"map"`}{Every row placed by two columns — coordinates, or any
#'     two numeric variables — with the selected ones filled in. This is the
#'     view for [design_spread()] and [design_spatial()], where what matters is
#'     how evenly the sample covers the space: compare a simple random sample's
#'     clumps and gaps against a spread one. For [design_spatial()] the
#'     outline of `region` is drawn as well, with each edge of a
#'     longitude/latitude region following its great circle under spherical
#'     geometry, which is how the design decides which rows are inside.}
#' }
#'
#' A design with `reserve` is shown by its main sample, without the reserves.
#'
#' Base graphics, so there is no plotting dependency to install.
#'
#' @param x A design object.
#' @param y The population data frame to draw against.
#' @param type `"selection"`, `"probability"` or `"map"`.
#' @param seed Optional seed for the `"selection"` draw, so the picture is
#'   reproducible.
#' @param ncol Dots per row in the `"selection"` grid. Defaults to whatever
#'   fills the panel at its current aspect ratio.
#' @param max_dots The most dots or steps to draw, which keeps them visible.
#'   A larger frame is shown in runs of consecutive rows: in `"selection"` a
#'   dot stands for a run and is filled when any row in it is selected, and in
#'   `"probability"` a step spans the lowest and highest probability in its
#'   run. `"map"` draws every selected row and an evenly spaced subset of the
#'   others. The counts under the title cover the whole frame.
#' @param main Title. Defaults to a description of the design.
#' @param palette Named list overriding any of `surface`, `ink`, `secondary`,
#'   `muted`, `recessive`, `accent`, `fill`, `rule`.
#' @param coords For `type = "map"`, the two columns to place rows by, `c(x,
#'   y)`, or one POINT geometry column. Defaults to the design's own `across`
#'   ([design_spread()]) or `coords` ([design_spatial()]); required for any
#'   other design. Date and date-time columns get date axes.
#' @param ... Passed to the underlying plot call.
#'
#' @return `x`, invisibly. Called for the plot.
#'
#' @examples
#' pop <- data.frame(
#'   id = 1:400,
#'   site = rep(c("a", "b", "c", "d"), times = c(200, 100, 60, 40)),
#'   cl = rep(paste0("c", 1:40), each = 10)
#' )
#'
#' op <- par(mfrow = c(2, 1))
#'
#' # Scattered, versus a visible lattice
#' plot(design_simple(n = 60), pop, seed = 1)
#' plot(design_systematic(interval = 7), pop, seed = 1)
#'
#' par(op)
#'
#' # Where the probabilities sit
#' plot(design_stratified("site", n = 60), pop, type = "probability")
#'
#' @seealso [draw()], [inclusion_prob()]
#' @export
plot.drawn_design <- function(x, y, type = c("selection", "probability", "map"),
                              seed = NULL, ncol = NULL, max_dots = 4000,
                              main = NULL, palette = NULL, coords = NULL, ...) {
  type <- check_choice(type, c("selection", "probability", "map"), "type")
  if (missing(y)) {
    stop("Pass the population data frame as the second argument: ",
         "plot(design, data).", call. = FALSE)
  }
  validate_data(y)
  max_dots <- check_count(max_dots, "max_dots", allow_zero = FALSE)
  if (!is.null(ncol)) ncol <- check_count(ncol, "ncol", allow_zero = FALSE)

  pal <- utils::modifyList(drawn_pal(), palette %||% list())
  main <- main %||% design_label(x)

  if (type == "probability") {
    return(plot_probability(x, y, max_dots, main, pal, ...))
  }
  # A design with `reserve` is shown by its main sample.
  if (type == "map") {
    n_pop <- nrow(y)
    shown <- if (n_pop > max_dots) {
      round(seq(1, n_pop, length.out = max_dots))
    } else {
      seq_len(n_pop)
    }
    sub <- if (length(shown) < n_pop) {
      paste0(fmt_n(length(shown)), " of ", fmt_n(n_pop), " rows shown")
    }
    plot_map(main_design(x), y, coords, shown, seed, main, sub, pal, ...)
  } else {
    plot_selection(main_design(x), y, max_dots, seed, ncol, main, pal, ...)
  }
  invisible(x)
}

#' Runs of consecutive frame rows, one per dot or step
#'
#' Every run holds `size` rows except the last, which holds what is left, so
#' there are at most `max_dots` runs. Runs are whole rows, so a periodic
#' selection cannot fall between the rows shown.
#'
#' @noRd
row_runs <- function(n_rows, max_dots) {
  size <- max(1, ceiling(n_rows / max_dots))
  k <- ceiling(n_rows / size)
  from <- (seq_len(k) - 1) * size + 1
  list(size = size, k = k, from = from, to = pmin(from + size - 1, n_rows))
}

#' The run each of `rows` falls in
#' @noRd
run_of <- function(rows, runs) (rows - 1) %/% runs$size + 1

#' How many times a design draws each frame row, for the selection and map
#' views
#'
#' A bootstrap is drawn as one replicate, as rows. The union of a thousand
#' replicates holds nearly every row, and the matrix outputs hold no row ids.
#'
#' @noRd
plot_counts <- function(x, y, seed) {
  key <- ".drawn_plot_id"
  if (key %in% names(y)) {
    stop("`y` already has a column called `", key, "`. Rename it.",
         call. = FALSE)
  }
  tagged <- y
  tagged[[key]] <- seq_len(nrow(y))
  if (inherits(x, "drawn_design_bootstrap")) {
    x$output <- "rows"
    x$n_replicates <- 1L
  }
  picked <- tryCatch({
    s <- draw(tagged, x, seed = seed)
    if (is.data.frame(s)) s[[key]] else unlist(s)
  }, error = function(e) {
    stop("Could not draw this design against `y`: ", conditionMessage(e),
         call. = FALSE)
  })
  tabulate(picked[!is.na(picked)], nbins = nrow(y))
}

#' The count under the title of the selection and map views
#'
#' Distinct rows for a design without replacement; draws and distinct rows
#' for one with replacement, where a row can be drawn more than once. A
#' bootstrap names the replicate shown.
#'
#' @noRd
selection_label <- function(design, counts) {
  boot <- inherits(design, "drawn_design_bootstrap")
  distinct <- sum(counts > 0)
  of <- fmt_n(length(counts))
  if (boot || samples_with_replacement(design) || any(counts > 1)) {
    draws <- sum(counts)
    out <- paste0(fmt_n(draws), " draw", if (draws != 1) "s", ", ",
                  fmt_n(distinct), " distinct row", if (distinct != 1) "s",
                  " of ", of)
  } else {
    out <- paste0(fmt_n(distinct), " of ", of, " rows selected")
  }
  if (boot) {
    out <- paste0("replicate 1 of ", fmt_n(design$n_replicates), "   |   ",
                  out)
  }
  out
}

#' One dot per run of rows, for the selection view
#'
#' A dot is filled when any row of its run was drawn, and ringed when any was
#' drawn more than once. `key` explains a dot that stands for more than one
#' row, and the ring, for the line under the panel.
#'
#' @noRd
selection_view <- function(counts, max_dots) {
  runs <- row_runs(length(counts), max_dots)
  filled <- tabulate(run_of(which(counts > 0), runs), runs$k) > 0
  ringed <- tabulate(run_of(which(counts > 1), runs), runs$k) > 0
  key <- NULL
  if (runs$size > 1) {
    key <- paste0("each dot is ", fmt_n(runs$size), " rows, filled when any ",
                  "of them is selected",
                  if (any(ringed)) {
                    " and ringed when one was drawn more than once"
                  })
  } else if (any(ringed)) {
    key <- "a ring marks a row drawn more than once"
  }
  list(filled = filled, ringed = ringed, size = runs$size, key = key)
}

#' What the probability view draws
#'
#' Each row is a plateau one unit wide, centred on its frame row, so the
#' first and last rows show as wide as the others. Above `max_dots` rows, a
#' run of rows is one step that passes through its lowest and highest
#' probability, in the order they come, so no unreachable or take-all row
#' drops out of the line. Zero marks go on every unreachable row, or on the
#' first such row of each run. The counts, range and mean come from every
#' row.
#'
#' @noRd
probability_view <- function(p, max_dots) {
  n <- length(p)
  runs <- row_runs(n, max_dots)
  left <- runs$from - 0.5
  right <- runs$to + 0.5
  if (runs$size == 1) {
    sx <- as.vector(rbind(left, right))
    sy <- rep(p, each = 2)
    zero <- which(p == 0)
  } else {
    ends <- vapply(seq_len(runs$k), function(r) {
      v <- p[runs$from[r]:runs$to[r]]
      lo <- which.min(v)
      hi <- which.max(v)
      if (!length(lo)) return(c(NA_real_, NA_real_, NA_real_))
      first_zero <- if (v[lo] == 0) runs$from[r] - 1 + lo else NA_real_
      c(if (lo <= hi) c(v[lo], v[hi]) else c(v[hi], v[lo]), first_zero)
    }, numeric(3))
    mid <- (left + right) / 2
    sx <- as.vector(rbind(left, mid, mid, right))
    sy <- as.vector(rbind(ends[1, ], ends[1, ], ends[2, ], ends[2, ]))
    zero <- ends[3, !is.na(ends[3, ])]
  }
  list(x = sx, y = sy, zero = zero, n = n, size = runs$size,
       n_zero = sum(p == 0, na.rm = TRUE),
       n_certain = sum(is_certain(p), na.rm = TRUE),
       low = min(p, na.rm = TRUE), top = max(p, na.rm = TRUE),
       mean = mean(p, na.rm = TRUE))
}

#' The line under the title of the probability view
#' @noRd
probability_label <- function(view) {
  rows <- function(k) paste0(fmt_n(k), " row", if (k != 1) "s")
  paste(c(
    if (view$n_zero) paste0(rows(view$n_zero), " unreachable"),
    if (view$n_certain) paste0(rows(view$n_certain), " taken with certainty"),
    paste0("range ", signif(view$low, 3), " to ", signif(view$top, 3))
  ), collapse = "   |   ")
}

#' Whole-row tick positions for the frame-row axis
#' @noRd
row_ticks <- function(n) {
  xt <- pretty(c(1, n), n = 5)
  xt[xt >= 1 & xt <= n & xt == round(xt)]
}

#' @noRd
fmt_n <- function(n) formatC(n, format = "d", big.mark = ",")

#' A one-line description of the design, for the default title
#' @noRd
design_label <- function(x) {
  p <- unclass(x)
  # One start and no reserves are the defaults, and say nothing about the
  # design.
  if (identical(p$n_starts, 1L)) p$n_starts <- NULL
  if (identical(p$reserve, 0L)) p$reserve <- NULL
  p <- p[!vapply(p, function(v) is.null(v) || inherits(v, c("sf", "sfc")),
                 logical(1))]
  keep <- intersect(c("strata", "clusters", "weights", "across", "size",
                      "time", "n", "reserve",
                      "n_clusters", "interval", "n_starts", "order_by",
                      "per_interval", "n_replicates",
                      "allocation", "method", "balanced", "within"), names(p))
  bits <- vapply(keep, function(k) paste0(k, " = ", fmt_param(p[[k]])),
                 character(1))
  paste0("design_", design_type(x), "(", paste(bits, collapse = ", "), ")")
}

#' Title above the panel, subtitle beneath it, both left-aligned
#' @noRd
draw_titles <- function(main, sub, pal) {
  graphics::mtext(main, side = 3, line = 1.55, adj = 0, cex = 0.9,
                  font = 2, col = pal$ink, xpd = NA)
  if (!is.null(sub) && nzchar(sub)) {
    graphics::mtext(sub, side = 3, line = 0.5, adj = 0, cex = 0.72,
                    col = pal$secondary, xpd = NA)
  }
}

#' @noRd
plot_selection <- function(x, y, max_dots, seed, ncol, main, pal, ...) {
  counts <- plot_counts(x, y, seed)
  view <- selection_view(counts, max_dots)
  sel <- view$filled

  k <- length(sel)
  # A second line under the grid needs room.
  bottom <- if (is.null(view$key)) 2.4 else 3.2
  op <- graphics::par(mar = c(bottom, 1.4, 3.6, 1.4), bg = pal$surface)
  on.exit(graphics::par(op), add = TRUE)

  # Shape the grid to the panel so it fills the space instead of floating in it.
  if (is.null(ncol)) {
    pin <- graphics::par("pin")
    aspect <- if (all(is.finite(pin)) && pin[2] > 0) pin[1] / pin[2] else 2
    ncol <- max(8L, min(k, as.integer(round(sqrt(k * aspect)))))
  }
  nr <- ceiling(k / ncol)
  gx <- ((seq_len(k) - 1L) %% ncol) + 1L
  gy <- nr - ((seq_len(k) - 1L) %/% ncol)

  graphics::plot(gx, gy, type = "n", axes = FALSE, ann = FALSE,
                 xlim = c(0.4, ncol + 0.6), ylim = c(0.4, nr + 0.6), ...)

  graphics::points(gx[!sel], gy[!sel], pch = 19, cex = 0.5, col = pal$recessive)
  # A surface-coloured ring goes on first, so adjacent selected dots stay
  # visually separate instead of merging into a blob.
  graphics::points(gx[sel], gy[sel], pch = 21, cex = 0.9, lwd = 1.5,
                   col = pal$surface, bg = pal$accent)
  # A second ring, in the accent, marks a row drawn more than once. It
  # shrinks with the spacing of the grid, so a dense grid is not covered in
  # large overlapping circles. A circle at cex 1 is about 0.45 * ps / 72
  # inches across.
  ring <- view$ringed
  if (any(ring)) {
    gap <- graphics::par("pin")[1] / diff(graphics::par("usr")[1:2])
    across <- 0.45 * graphics::par("ps") / 72 * graphics::par("cex")
    ring_cex <- min(1.5, max(1.1, 0.9 * gap / across))
    graphics::points(gx[ring], gy[ring], pch = 1, cex = ring_cex,
                     lwd = if (ring_cex < 1.5) 0.7 else 1, col = pal$accent)
  }

  draw_titles(main, selection_label(x, counts), pal)
  graphics::mtext("frame order, left to right, top to bottom", side = 1,
                  line = 0.8, adj = 0, cex = 0.68, col = pal$muted, xpd = NA)
  if (!is.null(view$key)) {
    graphics::mtext(view$key, side = 1, line = 1.6, adj = 0, cex = 0.68,
                    col = pal$muted, xpd = NA)
  }
  invisible(x)
}

#' @noRd
plot_probability <- function(x, y, max_dots, main, pal, ...) {
  # Keep the real reason. Assuming "no closed form" sent anyone with a missing
  # column, or a fixed `start`, off in entirely the wrong direction. A
  # bootstrap is not sent to the selection view, which shows one replicate
  # of many.
  p <- tryCatch(inclusion_prob(y, x), error = function(e) {
    stop("Cannot plot inclusion probabilities for this design: ",
         conditionMessage(e),
         if (!inherits(x, "drawn_design_bootstrap")) {
           "\nUse type = \"selection\" to see what it draws instead."
         },
         call. = FALSE)
  })
  view <- probability_view(p, max_dots)
  n <- view$n
  top <- view$top
  ylim <- c(0, if (!is.finite(top) || top <= 0) 1 else top * 1.14)
  xlim <- c(0.5, n + 0.5)

  op <- graphics::par(mar = c(3.6, 4.4, 3.6, 1.4), bg = pal$surface)
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(NA, xlim = xlim, ylim = ylim, axes = FALSE,
                 ann = FALSE, xaxs = "i", yaxs = "i", ...)

  # Rules sit behind the data, so the eye reads level before shape.
  ticks <- pretty(c(0, ylim[2]), n = 4)
  ticks <- ticks[ticks >= 0 & ticks <= ylim[2]]
  graphics::abline(h = ticks, col = pal$rule, lwd = 1)

  sx <- view$x
  sy <- view$y
  graphics::polygon(c(sx[1], sx, sx[length(sx)]), c(0, sy, 0),
                    col = pal$fill, border = NA)
  graphics::lines(sx, sy, col = pal$accent, lwd = 1.9)

  if (length(view$zero)) {
    graphics::points(view$zero, rep(0, length(view$zero)), pch = 19,
                     cex = 0.32, col = pal$muted)
  }

  m <- view$mean
  graphics::segments(xlim[1], m, xlim[2], m, col = pal$secondary, lty = 2,
                     lwd = 1.1)
  # Right-aligned above the end of the line, inside the panel.
  graphics::text(xlim[2], m, labels = paste0("mean ", signif(m, 3)),
                 adj = c(1, -0.4), cex = 0.68, col = pal$secondary, xpd = NA)

  graphics::axis(2, at = ticks, labels = format(ticks), las = 1,
                 cex.axis = 0.7, col = NA, col.ticks = NA,
                 col.axis = pal$secondary, line = -0.7)
  xt <- row_ticks(n)
  graphics::axis(1, at = xt, labels = fmt_n(xt), cex.axis = 0.7,
                 col = pal$rule, col.axis = pal$secondary, lwd = 1,
                 tck = -0.016)
  graphics::mtext(if (view$size > 1) {
                    paste0("row (frame order; a step per ", fmt_n(view$size),
                           " rows, through their lowest and highest)")
                  } else "row (frame order)",
                  side = 1, line = 2.1, cex = 0.72, col = pal$secondary)
  graphics::mtext("inclusion probability", side = 2, line = 3.0, cex = 0.72,
                  col = pal$secondary)

  draw_titles(main, probability_label(view), pal)
  invisible(x)
}

#' @noRd
plot_map <- function(x, y, coords, shown, seed, main, sub, pal, ...) {
  coords <- coords %||% x$across %||% x$coords
  if (!is.null(coords)) {
    validate_data(y, required_columns = coords)
    pos <- coordinate_columns(y, coords, "coords")
  }
  if (is.null(coords) || length(pos) < 2L) {
    stop("type = \"map\" needs two columns to place rows by. Pass ",
         "`coords = c(\"x\", \"y\")`, or the name of a POINT geometry ",
         "column.", call. = FALSE)
  }
  # A geometry column gives two coordinates, so a spread design's `across`
  # may hold more than two; the first two place the rows.
  pos <- pos[1:2]
  axis_ok <- function(v) is.numeric(v) || inherits(v, c("Date", "POSIXt"))
  if (!axis_ok(pos[[1]]) || !axis_ok(pos[[2]])) {
    stop("`", names(pos)[1], "` and `", names(pos)[2], "` must both be ",
         "numeric, dates or date-times.", call. = FALSE)
  }

  counts <- plot_counts(x, y, seed)
  sel <- counts > 0

  # Thinning a map would hide exactly the gaps it exists to show, so every
  # selected row is drawn and only the unselected background is thinned.
  bg <- setdiff(shown, which(sel))
  px <- pos[[1]]
  py <- pos[[2]]

  op <- graphics::par(mar = c(3.6, 3.8, 3.6, 1.4), bg = pal$surface)
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(px, py, type = "n", axes = FALSE, ann = FALSE, asp = NA, ...)
  graphics::box(col = pal$rule)
  # A spatial design's region, drawn under the rows so the rows inside it
  # and outside it can be told apart.
  if (inherits(x, "drawn_design_spatial") && identical(coords, x$coords)) {
    outline <- region_outline(x, y)
    graphics::polygon(outline$x, outline$y, border = pal$secondary,
                      col = NA, lwd = 1.2)
  }
  # Axis() labels a Date or date-time axis with dates.
  graphics::Axis(px, side = 1, cex.axis = 0.7, col = NA,
                 col.ticks = pal$rule, col.axis = pal$secondary, tck = -0.016)
  graphics::Axis(py, side = 2, cex.axis = 0.7, col = NA,
                 col.ticks = pal$rule, col.axis = pal$secondary, las = 1,
                 tck = -0.016)
  graphics::points(px[bg], py[bg], pch = 19, cex = 0.45, col = pal$recessive)
  graphics::points(px[sel], py[sel], pch = 21, cex = 1.05, lwd = 1.5,
                   col = pal$surface, bg = pal$accent)
  ring <- counts > 1
  graphics::points(px[ring], py[ring], pch = 1, cex = 1.7, lwd = 1,
                   col = pal$accent)
  graphics::mtext(names(pos)[1], side = 1, line = 2.1, cex = 0.72,
                  col = pal$secondary)
  graphics::mtext(names(pos)[2], side = 2, line = 2.6, cex = 0.72,
                  col = pal$secondary)
  repeats <- if (any(ring)) "ringed rows were drawn more than once"
  draw_titles(main, paste(c(selection_label(x, counts), repeats, sub),
                          collapse = "   |   "), pal)
  invisible(x)
}
