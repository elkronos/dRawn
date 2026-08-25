#' Hand a sample to the survey package
#'
#' Builds a [survey::svydesign()] object from a drawn sample, so the analysis
#' this package does not do — subpopulation estimates, regression, calibration,
#' quantiles with proper standard errors — can be done by the package that does.
#'
#' The two packages compute variance from different starting points. This one
#' uses the design's joint inclusion probabilities; `survey` reconstructs the
#' variance from the design's *shape*. So the job here is to express each design
#' in `survey`'s own terms rather than hand over a weight column and hope.
#'
#' @section What maps to what:
#' \tabular{lll}{
#'   **Design** \tab **Expressed as** \tab **Standard errors** \cr
#'   [design_simple()], [design_reservoir()], [design_spatial()] \tab `ids = ~1` with `fpc` the frame size \tab identical \cr
#'   [design_stratified()] \tab `strata` from the strata columns, `fpc` each stratum's size \tab identical \cr
#'   [design_temporal()] \tab `strata` from the sampling intervals, `fpc` each interval's size \tab identical \cr
#'   [design_cluster()] \tab `ids` the cluster column, `fpc` the number of clusters \tab identical \cr
#'   [design_weighted()], `"systematic"` \tab `ids = ~1` with `fpc` the frame size \tab identical \cr
#'   [design_weighted()], `"poisson"` \tab `survey::poisson_sampling()`, which models the random size \tab identical \cr
#'   [design_certainty()] \tab the certainty rows as their own stratum, taken whole \tab identical \cr
#'   [design_multistage()] \tab `ids` the cluster column \tab **differ by a few percent** \cr
#'   [design_systematic()] \tab `ids = ~1` with `fpc` the frame size \tab **`survey` returns one; this package declines** \cr
#' }
#'
#' "Identical" means to floating point, and is checked by this package's tests
#' against `survey::svytotal()`. The two exceptions are real and worth knowing:
#'
#' * **Multistage.** `survey` uses the ultimate-cluster approximation, which
#'   attributes all the variance to the first stage and ignores sampling within
#'   clusters. [ht_total()] uses the exact two-stage form. `survey`'s is the
#'   smaller of the two, by around 5–10% on a typical frame.
#' * **Systematic.** Most pairs of rows can never co-occur, so no
#'   design-unbiased variance exists and [ht_total()] returns `NA` with a note.
#'   `survey`, having only been told `ids = ~1`, computes the simple-random
#'   variance — which is the conservative substitute, not the design's own.
#'
#' Certainty rows arrive in a stratum where `n == N`, so `survey`'s own finite
#' population correction zeroes them out, matching this package's treatment.
#'
#' Two compositions have no single `survey` design and are refused rather than
#' approximated: [design_certainty()] over a cluster, multistage or Poisson
#' `rest`, where the certainty rows and the rest are different kinds of sampling
#' unit. [design_bootstrap()] is refused outright — it resamples the sample, so
#' there is no finite population for `svydesign()` to represent.
#'
#' @param sample A data frame returned by [draw()] with `weights = TRUE`.
#' @param ... Passed to [survey::svydesign()]. Anything named here overrides
#'   what the mapping above would have supplied, so `nest = TRUE` or a
#'   replacement `fpc` is yours to set.
#'
#' @return A `survey.design` object.
#'
#' @examplesIf requireNamespace("survey", quietly = TRUE)
#' set.seed(1)
#' pop <- data.frame(
#'   id = 1:400,
#'   site = rep(c("a", "b", "c", "d"), times = c(200, 100, 60, 40)),
#'   spend = round(stats::runif(400, 10, 500))
#' )
#' s <- draw(pop, design_stratified("site", n = 60), seed = 1, weights = TRUE)
#'
#' des <- as_svydesign(s)
#' survey::svytotal(~spend, des)
#'
#' # The same total, and the same standard error
#' ht_total(s, "spend")
#'
#' # Now the analysis this package does not do
#' survey::svyby(~spend, ~site, des, survey::svymean)
#'
#' @seealso [ht_total()], [ht_mean()]
#' @export
as_svydesign <- function(sample, ...) {
  require_suggested("survey", "handing a sample to the survey package")
  if (!is.data.frame(sample) || is.null(attr(sample, "drawn_design"))) {
    stop("`sample` must be a data frame from draw(..., weights = TRUE).",
         call. = FALSE)
  }
  design <- attr(sample, "drawn_design")
  pop <- attr(sample, "drawn_population")
  if (design_type(design) == "bootstrap") {
    stop("A bootstrap design has no survey equivalent: it resamples the ",
         "sample rather than drawing a probability sample of a finite ",
         "population.", call. = FALSE)
  }

  dat <- as.data.frame(sample)
  parts <- survey_parts(design, dat, pop, rows = attr(sample, "drawn_rows"))
  dat$.fpc <- parts$fpc
  if (!is.null(parts$stratum)) dat$.stratum <- parts$stratum
  args <- c(list(data = dat, weights = stats::as.formula("~.weight"),
                 ids = parts$ids, fpc = stats::as.formula("~.fpc")),
            parts$extra)
  if (!is.null(parts$stratum)) {
    args$strata <- stats::as.formula("~.stratum")
  }
  do.call(survey::svydesign, utils::modifyList(args, list(...)))
}

#' Strata, population counts and cluster ids for one design
#'
#' `survey` reconstructs a variance from the design's shape, not from the
#' inclusion probabilities, so each design has to be expressed in its own terms:
#' a stratified design as strata, a cluster design as ids, a temporal design as
#' strata over its intervals, and a certainty design as a stratum that was taken
#' whole -- where `n == N` makes `survey`'s own correction zero out, matching
#' this package's treatment of certainty rows.
#'
#' @noRd
refuse_nested <- function(prefix, what, why) {
  if (nzchar(prefix)) {
    stop("A certainty design over ", what, " has no single survey equivalent: ",
         why, ". Estimate it with ht_total() or ht_mean() instead.",
         call. = FALSE)
  }
  invisible(NULL)
}

#' @noRd
survey_parts <- function(design, dat, pop, rows, prefix = "") {
  one <- stats::as.formula("~1")
  by_group <- function(cols) {
    sizes <- table(group_key(pop, cols))
    lab <- as.character(group_key(dat, cols))
    list(stratum = paste0(prefix, lab),
         fpc = as.integer(sizes[lab]), ids = one)
  }

  switch(design_type(design),
    stratified = by_group(design$strata),
    temporal = {
      sizes <- table(temporal_bucket(design, pop))
      lab <- temporal_bucket(design, dat)
      list(stratum = paste0(prefix, lab), fpc = as.integer(sizes[lab]),
           ids = one)
    },
    weighted = {
      if (design$method != "poisson") {
        return(list(stratum = NULL, fpc = rep(nrow(pop), nrow(dat)), ids = one))
      }
      # Poisson sampling has a random size, so no finite population correction
      # over a fixed n describes it. `survey` models it directly.
      refuse_nested(prefix, "a Poisson `rest`",
                    "its sample size is random, which no stratum can express")
      list(stratum = NULL, fpc = dat$.prob, ids = one,
           extra = list(pps = survey::poisson_sampling(dat$.prob)))
    },
    cluster = ,
    multistage = {
      refuse_nested(prefix, "a cluster or multistage `rest`",
                    paste0("the certainty rows and the clusters are different ",
                           "kinds of sampling unit"))
      list(stratum = NULL,
           fpc = rep(length(unique(pop[[design$clusters]])), nrow(dat)),
           ids = stats::as.formula(paste0("~", design$clusters)))
    },
    spatial = {
      # The design's population is the rows inside the region, not the frame.
      # Handing `survey` the frame size overstated its standard error by 6-11%.
      inside <- sum(spatial_inside(design, pop))
      list(stratum = NULL, fpc = rep(inside, nrow(dat)), ids = one)
    },
    certainty = {
      sp <- certainty_split(design, pop)
      certain_pop <- sp$keep[sp$take]
      below <- sp$data[sp$rest, , drop = FALSE]
      # Membership of the certainty part, not `.prob == 1`: `rest` can put a
      # row at probability 1 of its own accord -- a dominant PPS unit, or an
      # `n` equal to the rows below the threshold -- and filing those into the
      # certainty stratum made `svydesign()` reject the whole design with
      # "FPC implies >100% sampling in some strata".
      is_certain <- rows %in% sp$keep[sp$take]

      out <- list(stratum = character(nrow(dat)), fpc = numeric(nrow(dat)),
                  ids = one)
      out$stratum[is_certain] <- paste0(prefix, "(certain)")
      out$fpc[is_certain] <- length(certain_pop)
      if (any(!is_certain)) {
        rest <- survey_parts(design$rest, dat[!is_certain, , drop = FALSE],
                             below, rows = match(rows[!is_certain],
                                                 sp$keep[sp$rest]),
                             prefix = paste0(prefix, "(rest) "))
        out$stratum[!is_certain] <- rest$stratum %||%
          paste0(prefix, "(rest)")
        out$fpc[!is_certain] <- rest$fpc
      }
      out
    },
    list(stratum = NULL, fpc = rep(nrow(pop), nrow(dat)), ids = one)
  )
}
