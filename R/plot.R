#' Picture a sampling design
#'
#' Draws the design against a population so you can see what it selects. Two
#' views, and both answer questions that a table of numbers answers slowly:
#'
#' \describe{
#'   \item{`"selection"`}{Every row of the frame as a dot, laid out left to
#'     right and top to bottom in frame order, with the selected ones filled in.
#'     Designs look distinct: simple random sampling scatters, systematic makes
#'     a lattice, cluster sampling takes solid contiguous runs, and
#'     probability-proportional-to-size thickens wherever the weight is large.
#'     If your frame is sorted by something meaningful, an unintended pattern
#'     shows up here immediately.}
#'   \item{`"probability"`}{Each row's inclusion probability against its
#'     position. Flat means every row had the same chance; steps mean strata;
#'     a slope means size-proportional selection. Rows the design can never
#'     reach sit at zero, which is usually the thing you wanted to find out.}
#' }
#'
#' Base graphics, so there is no plotting dependency to install.
#'
#' @param x A design object.
#' @param y The population data frame to draw against.
#' @param type `"selection"` or `"probability"`.
#' @param seed Optional seed for the `"selection"` draw, so the picture is
#'   reproducible.
#' @param ncol Dots per row in the `"selection"` grid.
#' @param max_dots Frames larger than this are shown as an evenly spaced
#'   subset, with a note under the plot. Keeps individual dots visible.
#' @param col,col_bg Colours for selected and unselected rows.
#' @param main Title. Defaults to a description of the design.
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
#' op <- par(mfrow = c(2, 1), mar = c(2, 1, 2, 1))
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
                              seed = NULL, ncol = 40, max_dots = 4000,
                              col = "#17594A", col_bg = "#C9D2CE",
                              main = NULL, ...) {
  type <- match.arg(type)
  if (missing(y)) {
    stop("Pass the population data frame as the second argument: ",
         "plot(design, data).", call. = FALSE)
  }
  validate_data(y)
  max_dots <- check_count(max_dots, "max_dots", allow_zero = FALSE)
  ncol <- check_count(ncol, "ncol", allow_zero = FALSE)

  n_pop <- nrow(y)
  shown <- if (n_pop > max_dots) {
    round(seq(1, n_pop, length.out = max_dots))
  } else {
    seq_len(n_pop)
  }
  sub <- if (length(shown) < n_pop) {
    paste0(format(length(shown), big.mark = ","), " of ",
           format(n_pop, big.mark = ","), " rows shown")
  } else {
    NULL
  }
  main <- main %||% design_label(x)

  if (type == "probability") {
    return(plot_probability(x, y, shown, main, sub, col, ...))
  }
  plot_selection(x, y, shown, seed, ncol, main, sub, col, col_bg, ...)
}

#' @noRd
design_label <- function(x) {
  p <- unclass(x)
  p <- p[!vapply(p, function(v) is.null(v) || inherits(v, c("sf", "sfc")),
                 logical(1))]
  keep <- intersect(c("n", "n_clusters", "interval", "per_interval",
                      "n_replicates", "strata", "clusters", "weights",
                      "allocation", "method", "balanced"), names(p))
  bits <- vapply(keep, function(k) paste0(k, " = ", fmt_param(p[[k]])),
                 character(1))
  paste0("design_", design_type(x), "(", paste(bits, collapse = ", "), ")")
}

#' @noRd
plot_selection <- function(x, y, shown, seed, ncol, main, sub, col, col_bg, ...) {
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
  nr <- ceiling(k / ncol)
  gx <- ((seq_len(k) - 1L) %% ncol) + 1L
  gy <- nr - ((seq_len(k) - 1L) %/% ncol)

  graphics::plot(gx, gy, type = "n", axes = FALSE, xlab = "", ylab = "",
                 main = main, xlim = c(0.5, ncol + 0.5),
                 ylim = c(0.5, nr + 0.5), ...)
  graphics::points(gx[!sel], gy[!sel], pch = 19, cex = 0.6, col = col_bg)
  graphics::points(gx[sel], gy[sel], pch = 19, cex = 0.75, col = col)
  graphics::mtext(
    paste0(sum(sel), " of ", format(length(shown), big.mark = ","),
           " selected", if (!is.null(sub)) paste0("  |  ", sub) else ""),
    side = 1, line = 0.4, cex = 0.72, col = "#5B6763"
  )
  invisible(x)
}

#' @noRd
plot_probability <- function(x, y, shown, main, sub, col, ...) {
  p <- tryCatch(inclusion_prob(y, x), error = function(e) {
    stop("This design has no closed-form inclusion probability, so there is ",
         "nothing to plot. Use type = \"selection\", or see ?inclusion_prob.",
         call. = FALSE)
  })
  pv <- p[shown]

  graphics::plot(seq_along(pv), pv, type = "h", lwd = 1, col = col,
                 xlab = "row (frame order)", ylab = "inclusion probability",
                 main = main, ylim = c(0, max(1e-9, max(pv))), bty = "n", ...)
  graphics::abline(h = mean(pv), lty = 2, col = "#5B6763")
  graphics::mtext(
    paste0("mean ", signif(mean(pv), 3),
           if (any(pv == 0)) paste0("  |  ", sum(pv == 0), " row(s) unreachable")
           else "",
           if (!is.null(sub)) paste0("  |  ", sub) else ""),
    side = 1, line = 3.8, cex = 0.7, col = "#5B6763"
  )
  invisible(x)
}
