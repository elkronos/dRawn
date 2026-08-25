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
#'   \item{[design_certainty()]}{Everything above a threshold, plus a sample of
#'     the rest.}
#'   \item{[design_reservoir()]}{A fixed-size sample from a stream.}
#'   \item{[design_bootstrap()]}{Resampled replicates.}
#'   \item{[design_temporal()]}{A share of each time interval.}
#'   \item{[design_spatial()]}{Rows inside a region.}
#' }
#'
#' @section What every design guarantees:
#' Arguments mean the same thing everywhere they appear:
#'
#' * `n` is always the **total** number of rows drawn, never a per-group figure.
#'   [design_temporal()] uses `per_interval` precisely because that one *is* per
#'   group, and [design_cluster()] has no `n` at all because the row count
#'   follows from which clusters were selected.
#' * `allocation` always says how a total is split across groups, either
#'   `"proportional"` or `"equal"`. Splitting is exact: the largest-remainder
#'   method is used, so a request for 60 rows returns 60 rather than drifting
#'   with per-group rounding.
#' * `na_rm` always decides whether rows with a missing key are dropped or raise
#'   an error. It is never silently assumed either way.
#' * `replace` always means sampling with replacement within whatever group the
#'   design works on, and always rules out `draw(weights = TRUE)`: an inclusion
#'   probability describes distinct units, and a sample holding duplicates
#'   cannot be weighted by one.
#' * [draw()] always restores the caller's random number stream before
#'   returning, and always gives back a data frame with the input's class and
#'   column order.
#' * Rows come back in frame order for every design that selects a set of rows.
#'   The two exceptions are the ones where draw order is meaningful:
#'   [design_simple()] and [design_weighted()] return rows in the order they
#'   were drawn, and [design_bootstrap()] returns replicates in order with a
#'   leading `.replicate` column.
#'
#' @section Estimating from a sample:
#' Because a design is a value, it can report its own inclusion probabilities
#' before any sampling happens. [inclusion_prob()] gives first-order
#' probabilities and [sampling_weight()] their reciprocals — the number of
#' population rows each sampled row stands for. [joint_prob()] gives
#' second-order probabilities, and [ht_total()] and
#' [ht_mean()] combine them into a population total or mean with a standard
#' error. [deff()] reports what the design cost in precision against simple
#' random sampling, and [sample_summary()] reports what was actually drawn
#' against what was in the frame.
#'
#' Going the other way, [plan_size()] solves for the sample size a given margin
#' of error requires — the step *before* choosing a design. For analysis this
#' package does not do, [as_svydesign()] hands a sample to the `survey` package.
#'
#' Not every design has a closed form for these — see [inclusion_prob()] for
#' which, why, and what to do instead. Those that do not say so rather than
#' returning an approximation.
#'
#' @section Choosing one:
#' \describe{
#'   \item{Every unit equally likely}{[design_simple()], or
#'     [design_systematic()] when the frame has a useful order, or
#'     [design_reservoir()] when it does not fit in memory.}
#'   \item{Guaranteed coverage of subgroups}{[design_stratified()], with
#'     `min_per_stratum` if rare groups must appear.}
#'   \item{Fieldwork cost matters more than efficiency}{[design_cluster()] or
#'     [design_multistage()] — visiting five sites is cheaper than visiting
#'     fifty, at the cost of precision.}
#'   \item{Large units matter more}{[design_weighted()], with
#'     `method = "systematic"` or `"poisson"` if you intend to estimate.}
#'   \item{A few units dominate the total}{[design_certainty()] — take those
#'     with certainty and sample the tail, which removes them from the variance
#'     entirely.}
#'   \item{Coverage across time or space}{[design_temporal()],
#'     [design_spatial()].}
#'   \item{Uncertainty of a statistic, not a population total}{
#'     [design_bootstrap()].}
#' }
#'
#' @name designs
#' @seealso [draw()], [inclusion_prob()], [ht_total()], [plan_size()],
#'   [sample_summary()]
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

