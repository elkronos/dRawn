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
#'   \item{[design_stratified()]}{A share of each stratum, at random or
#'     walked in order.}
#'   \item{[design_systematic()]}{Every *k*-th row, or `n` rows at an even
#'     step, from one random start or several.}
#'   \item{[design_cluster()]}{Whole clusters.}
#'   \item{[design_multistage()]}{Clusters, at random or in proportion to
#'     size, then rows within them.}
#'   \item{[design_weighted()]}{Rows with probability governed by a weight,
#'     within strata if asked.}
#'   \item{[design_certainty()]}{Everything above a threshold, plus a sample of
#'     the rest.}
#'   \item{[design_reservoir()]}{A fixed-size sample from a stream.}
#'   \item{[design_bootstrap()]}{Resampled replicates.}
#'   \item{[design_temporal()]}{A share of each time interval.}
#'   \item{[design_spatial()]}{Rows inside a region.}
#'   \item{[design_spread()]}{A sample spread evenly across a map or across
#'     auxiliary variables, within strata if asked.}
#' }
#'
#' @section What every design guarantees:
#' Arguments mean the same thing everywhere they appear:
#'
#' * `n` is always the **total** number of rows drawn, never a per-group figure.
#'   [design_temporal()] uses `per_interval` precisely because that one *is* per
#'   group, and [design_cluster()] has no `n` at all because the row count
#'   follows from which clusters were selected.
#' * `allocation` always says how a total is split across groups:
#'   `"proportional"`, `"equal"` or, in [design_stratified()] only,
#'   `"neyman"`. [design_stratified()], [design_spread()] and
#'   [design_weighted()] also take the count for each stratum in place of a
#'   rule, and `n` is then their total.
#'   Splitting is exact: the largest-remainder method is used, so a request
#'   for 60 rows returns 60 rather than drifting with per-group rounding.
#' * `na_rm` always decides whether rows with a missing key are dropped or raise
#'   an error. It is never silently assumed either way.
#' * Groups always come in key order, whether strata, `by` domains or the
#'   rows of [sample_summary()]: a factor's levels, or the sorted values of any
#'   other column, with text sorted byte by byte as in the C locale. A seed
#'   therefore draws the same rows in every locale, and an `order_by` column
#'   of text is walked in the same order everywhere.
#' * `replace` always means sampling with replacement within whatever group the
#'   design works on, and always rules out `draw(weights = TRUE)`: an inclusion
#'   probability describes distinct units, and a sample holding duplicates
#'   cannot be weighted by one.
#' * [draw()] always restores the caller's random number stream before
#'   returning, and always gives back a data frame with the input's class and
#'   column order. [design_bootstrap()] with `output = "indices"` or
#'   `"weights"` gives a matrix instead, and [design_reservoir()] a list
#'   when `data` is not a data frame, unless `batches = TRUE`.
#' * Rows come back in frame order for every design that selects a set of rows,
#'   including [design_systematic()] with `order_by`.
#'   The two exceptions are the ones where draw order is meaningful:
#'   [design_simple()] and [design_weighted()] return rows in the order they
#'   were drawn, and [design_bootstrap()] returns replicates in order with a
#'   leading `.replicate` column. A design with `reserve` returns its main
#'   sample and then its reserves, in the order they are used, with a
#'   leading `.reserve` column.
#'
#' @section Estimating from a sample:
#' Because a design is a value, it can report its own inclusion probabilities
#' before any sampling happens. [inclusion_prob()] gives first-order
#' probabilities and [sampling_weight()] their reciprocals — the number of
#' population rows each sampled row stands for, and [certainty_set()] the rows
#' a design that selects in proportion to size takes with certainty.
#' [allocation_table()] shows how a design with strata splits `n` over the
#' strata of a frame, and which strata it holds at the floor, takes whole or
#' leaves empty. [joint_prob()] gives
#' second-order probabilities, and [ht_total()] and [ht_mean()] turn a sample
#' into a population total, mean or proportion — overall or by domain — with a
#' standard error from the variance estimator suited to the design and a
#' confidence interval on the design's degrees of freedom. [deff()] reports
#' what the design cost in precision against simple random sampling, and
#' [sample_summary()] reports what was actually drawn against what was in the
#' frame. [add_fieldwork()] joins fieldwork results to a sample by an id,
#' [use_reserves()] replaces the units of a list from [design_simple()] or
#' [design_spread()] with `reserve` that were found unavailable,
#' [as_drawn()] rebuilds a sample read back from a file or a database, and
#' [compact_sample()] drops the frame columns its design does not read before
#' it is saved or shared.
#'
#' Going the other way, [plan_size()] solves for the sample size a given margin
#' of error requires — the step *before* choosing a design. For analysis this
#' package does not do, [as_svydesign()] hands a sample to the `survey` package,
#' [replicate_weights()] builds replicate weights for it from the design, and
#' [frame_totals()] gives the frame's totals to calibrate it to.
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
#'     fifty, at the cost of precision. With `size`, [design_multistage()]
#'     selects the sites in proportion to size and takes the same number of
#'     rows from each.}
#'   \item{Large units matter more}{[design_weighted()], with
#'     `method = "systematic"` or `"poisson"` if you intend to estimate, and
#'     `strata` to sample in proportion to size within each stratum.}
#'   \item{A few units dominate the total}{[design_certainty()] — take those
#'     with certainty and sample the tail, which removes them from the variance
#'     entirely.}
#'   \item{Coverage across time or space}{[design_temporal()] for time;
#'     [design_spread()] for an even spread over a map, which is usually far
#'     more precise than a simple random sample when what you measure varies
#'     smoothly over space; [design_spatial()] to restrict sampling to a
#'     region.}
#'   \item{Balance on known covariates without choosing strata}{
#'     [design_spread()] across those covariates.}
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
  shown <- lapply(params, function(v) if (!is.null(v)) fmt_param(v))
  # A temporal design without `tz` takes its zone from the data.
  if (inherits(x, "drawn_design_temporal") && is.null(params$tz)) {
    shown$tz <- "from the time column, else UTC"
  }
  shown <- shown[!vapply(shown, is.null, logical(1))]
  if (length(shown) == 0L) {
    return(invisible(x))
  }
  width <- max(nchar(names(shown)))
  for (nm in names(shown)) {
    cat("  ", formatC(nm, width = -width), "  ", shown[[nm]], "\n", sep = "")
  }
  invisible(x)
}

#' @noRd
fmt_param <- function(v) {
  if (inherits(v, c("sf", "sfc"))) {
    return(paste0("<", class(v)[1], ">"))
  }
  if (is.function(v)) return("<function>")
  if (is_design(v)) return(design_label(v))
  # Counts given per stratum as `allocation`.
  if (is.data.frame(v) || (is.numeric(v) && !is.null(names(v)))) {
    return(fmt_counts(v))
  }
  if (inherits(v, "POSIXt") || inherits(v, "Date")) return(format(v))
  if (is.character(v)) {
    return(paste0("\"", v, "\"", collapse = ", "))
  }
  if (length(v) > 6L) {
    return(paste0("<", class(v)[1], "[", length(v), "]>"))
  }
  paste0(format(v), collapse = ", ")
}

