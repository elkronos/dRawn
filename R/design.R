#' Sampling designs
#'
#' A design describes *how* to sample, independently of *what* you are sampling.
#' Build one with a `design_*()` constructor, then hand it to [draw()] with the
#' data. The same design can be reused across data sets, stored, printed and
#' passed around.
#'
#' @section The designs:
#' \describe{
#'   \item{[design_simple()]}{Rows uniformly at random.}
#'   \item{[design_stratified()]}{A share of each stratum.}
#'   \item{[design_systematic()]}{Every *k*-th row.}
#'   \item{[design_cluster()]}{Whole clusters.}
#'   \item{[design_multistage()]}{Clusters, then rows within them.}
#'   \item{[design_weighted()]}{Rows with probability governed by a weight.}
#'   \item{[design_reservoir()]}{A fixed-size sample from a stream.}
#'   \item{[design_bootstrap()]}{Resampled replicates.}
#'   \item{[design_temporal()]}{A share of each time interval.}
#'   \item{[design_spatial()]}{Rows inside a region.}
#' }
#'
#' @section What every design guarantees:
#' Arguments mean the same thing everywhere they appear. `n` is always the total
#' number of rows drawn, never a per-group figure; `allocation` always says how a
#' total is split across groups; `na_rm` always decides whether missing keys are
#' dropped or raise an error; and [draw()] always restores the caller's random
#' number stream before returning.
#'
#' @name designs
#' @seealso [draw()]
NULL

#' Build a design object
#' @noRd
new_design <- function(type, params) {
  structure(
    params,
    class = c(paste0("drawn_design_", type), "drawn_design")
  )
}

#' Is this a sampling design?
#'
#' @param x An object.
#' @return `TRUE` if `x` was built by one of the `design_*()` constructors.
#' @examples
#' is_design(design_simple(n = 10))
#' is_design(mtcars)
#' @export
is_design <- function(x) inherits(x, "drawn_design")

#' @noRd
design_type <- function(x) {
  sub("^drawn_design_", "", class(x)[1])
}

#' @export
print.drawn_design <- function(x, ...) {
  cat("<sampling design: ", design_type(x), ">\n", sep = "")
  params <- unclass(x)
  params <- params[!vapply(params, is.null, logical(1))]
  if (length(params) == 0L) {
    return(invisible(x))
  }
  width <- max(nchar(names(params)))
  for (nm in names(params)) {
    cat("  ", formatC(nm, width = -width), "  ", fmt_param(params[[nm]]), "\n",
        sep = "")
  }
  invisible(x)
}

#' @noRd
fmt_param <- function(v) {
  if (inherits(v, c("sf", "sfc"))) {
    return(paste0("<", class(v)[1], ">"))
  }
  if (is.function(v)) return("<function>")
  if (inherits(v, "POSIXt") || inherits(v, "Date")) return(format(v))
  if (is.character(v)) {
    return(paste0("\"", v, "\"", collapse = ", "))
  }
  if (length(v) > 6L) {
    return(paste0("<", class(v)[1], "[", length(v), "]>"))
  }
  paste0(format(v), collapse = ", ")
}

