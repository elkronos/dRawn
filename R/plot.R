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
#'     meaningful, an unintended pattern shows up immediately.}
#'   \item{`"probability"`}{Each row's inclusion probability across the frame,
#'     as a step. Flat means every row had the same chance; plateaus mean
#'     strata; a rise means size-proportional selection. Rows the design can
#'     never reach sit at zero and are marked, which is usually the thing worth
#'     finding out.}
#' }
#'
#' Base graphics, so there is no plotting dependency to install.
#'
#' @param x A design object.
#' @param y The population data frame to draw against.
#' @param type `"selection"` or `"probability"`.
#' @param seed Optional seed for the `"selection"` draw, so the picture is
#'   reproducible.
#' @param ncol Dots per row in the `"selection"` grid. Defaults to whatever
#'   fills the panel at its current aspect ratio.
#' @param max_dots Frames larger than this are shown as an evenly spaced
#'   subset, noted under the title. Keeps individual dots visible.
#' @param main Title. Defaults to a description of the design.
#' @param palette Named list overriding any of `surface`, `ink`, `secondary`,
#'   `muted`, `recessive`, `accent`, `fill`, `rule`.
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
plot.drawn_design <- function(x, y, type = c("selection", "probability"),
                              seed = NULL, ncol = NULL, max_dots = 4000,
                              main = NULL, palette = NULL, ...) {
  type <- match.arg(type)
  if (missing(y)) {
    stop("Pass the population data frame as the second argument: ",
         "plot(design, data).", call. = FALSE)
  }
  validate_data(y)
  max_dots <- check_count(max_dots, "max_dots", allow_zero = FALSE)
  if (!is.null(ncol)) ncol <- check_count(ncol, "ncol", allow_zero = FALSE)

  pal <- utils::modifyList(drawn_pal(), palette %||% list())
  n_pop <- nrow(y)
  shown <- if (n_pop > max_dots) {
    round(seq(1, n_pop, length.out = max_dots))
  } else {
    seq_len(n_pop)
  }
  sub <- if (length(shown) < n_pop) {
    paste0(fmt_n(length(shown)), " of ", fmt_n(n_pop), " rows shown")
  } else {
    NULL
  }
  main <- main %||% design_label(x)

  if (type == "probability") {
    return(plot_probability(x, y, shown, main, sub, pal, ...))
  }
  plot_selection(x, y, shown, seed, ncol, main, sub, pal, ...)
}

#' @noRd
fmt_n <- function(n) formatC(n, format = "d", big.mark = ",")

#' A one-line description of the design, for the default title
#' @noRd
design_label <- function(x) {
  p <- unclass(x)
  p <- p[!vapply(p, function(v) is.null(v) || inherits(v, c("sf", "sfc")),
                 logical(1))]
  keep <- intersect(c("strata", "clusters", "weights", "time", "n",
                      "n_clusters", "interval", "per_interval", "n_replicates",
                      "allocation", "method", "balanced"), names(p))
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
plot_selection <- function(x, y, shown, seed, ncol, main, sub, pal, ...) {
  key <- ".drawn_plot_id"
  if (key %in% names(y)) {
    stop("`y` already has a column called `", key, "`. Rename it.",
         call. = FALSE)
  }
  tagged <- y
  tagged[[key]] <- seq_len(nrow(y))

  picked <- tryCatch({
    s <- draw(tagged, x, seed = seed)
    if (is.data.frame(s)) s[[key]] else unlist(s)
  }, error = function(e) {
    stop("Could not draw this design against `y`: ", conditionMessage(e),
         call. = FALSE)
  })

  sel <- rep(FALSE, nrow(y))
  sel[picked[!is.na(picked)]] <- TRUE
  sel <- sel[shown]

  k <- length(shown)
  op <- graphics::par(mar = c(2.4, 1.4, 3.6, 1.4), bg = pal$surface)
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

  draw_titles(main, paste0(fmt_n(sum(sel)), " of ", fmt_n(k), " rows selected",
                           if (!is.null(sub)) paste0("   |   ", sub) else ""),
              pal)
  graphics::mtext("frame order, left to right, top to bottom", side = 1,
                  line = 0.8, adj = 0, cex = 0.68, col = pal$muted, xpd = NA)
  invisible(x)
}

#' @noRd
plot_probability <- function(x, y, shown, main, sub, pal, ...) {
  # Keep the real reason. Assuming "no closed form" sent anyone with a missing
  # column, or a fixed `start`, off in entirely the wrong direction.
  p <- tryCatch(inclusion_prob(y, x), error = function(e) {
    stop("Cannot plot inclusion probabilities for this design: ",
         conditionMessage(e),
         "\nUse type = \"selection\" to see what it draws instead.",
         call. = FALSE)
  })
  pv <- p[shown]
  k <- length(pv)
  top <- max(pv, na.rm = TRUE)
  ylim <- c(0, if (!is.finite(top) || top <= 0) 1 else top * 1.14)

  op <- graphics::par(mar = c(3.6, 4.4, 3.6, 1.4), bg = pal$surface)
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(NA, xlim = c(1, max(k, 2)), ylim = ylim, axes = FALSE,
                 ann = FALSE, xaxs = "i", yaxs = "i", ...)

  # Rules sit behind the data, so the eye reads level before shape.
  ticks <- pretty(c(0, ylim[2]), n = 4)
  ticks <- ticks[ticks >= 0 & ticks <= ylim[2]]
  graphics::abline(h = ticks, col = pal$rule, lwd = 1)

  if (k >= 2L) {
    sx <- rep(seq_len(k), each = 2)[-1]
    sy <- rep(pv, each = 2)[-(2 * k)]
    graphics::polygon(c(sx[1], sx, sx[length(sx)]), c(0, sy, 0),
                      col = pal$fill, border = NA)
    graphics::lines(sx, sy, col = pal$accent, lwd = 1.9)
  }

  zero <- which(pv == 0)
  if (length(zero)) {
    graphics::points(zero, rep(0, length(zero)), pch = 19, cex = 0.32,
                     col = pal$muted)
  }

  m <- mean(pv, na.rm = TRUE)
  graphics::segments(1, m, k, m, col = pal$secondary, lty = 2, lwd = 1.1)
  graphics::text(k, m, labels = paste0("mean ", signif(m, 3), " "), pos = 3,
                 offset = 0.2, adj = 1, cex = 0.68, col = pal$secondary,
                 xpd = NA)

  graphics::axis(2, at = ticks, labels = format(ticks), las = 1,
                 cex.axis = 0.7, col = NA, col.ticks = NA,
                 col.axis = pal$secondary, line = -0.7)
  xt <- pretty(c(1, k), n = 5)
  xt <- xt[xt >= 1 & xt <= k]
  graphics::axis(1, at = xt, labels = fmt_n(xt), cex.axis = 0.7,
                 col = pal$rule, col.axis = pal$secondary, lwd = 1,
                 tck = -0.016)
  graphics::mtext(if (length(shown) < nrow(y)) "row (frame order, thinned)"
                  else "row (frame order)",
                  side = 1, line = 2.1, cex = 0.72, col = pal$secondary)
  graphics::mtext("inclusion probability", side = 2, line = 3.0, cex = 0.72,
                  col = pal$secondary)

  draw_titles(main, paste0(
    if (length(zero)) {
      paste0(fmt_n(length(zero)), " row", if (length(zero) > 1) "s" else "",
             " unreachable   |   ")
    } else "",
    "range ", signif(min(pv, na.rm = TRUE), 3), " to ", signif(top, 3),
    if (!is.null(sub)) paste0("   |   ", sub) else ""), pal)
  invisible(x)
}
