#' Hand a sample to the survey package
#'
#' Builds a [survey::svydesign()] object from a drawn sample, so the analysis
#' this package does not do — regression, calibration, quantiles with proper
#' standard errors, replicate weights — can be done by the package that does.
#'
#' The two packages compute variance from different starting points. This one
#' uses the design's inclusion probabilities; `survey` reconstructs the
#' variance from the design's *shape*. So the job here is to express each design
#' in `survey`'s own terms rather than hand over a weight column and hope.
#'
#' @section What maps to what:
#' \tabular{lll}{
#'   **Design** \tab **Expressed as** \tab **Standard errors** \cr
#'   [design_simple()], [design_spatial()] \tab `ids = ~1` with `fpc` the rows the design can reach \tab identical \cr
#'   [design_reservoir()] \tab `ids = ~1` with `fpc` the rows the stream reached \tab identical \cr
#'   [design_stratified()] \tab `strata` from the strata columns, `fpc` each stratum's size \tab identical \cr
#'   [design_temporal()] \tab `strata` from the sampling intervals, `fpc` each interval's size \tab identical \cr
#'   [design_cluster()] \tab `ids` a `.cluster` column numbering the clusters, `fpc` the number of clusters \tab identical \cr
#'   [design_multistage()] \tab two stages, `ids` the `.cluster` numbers then the rows, `fpc` the number of clusters and each cluster's size \tab identical \cr
#'   [design_weighted()], `"poisson"` \tab `survey::poisson_sampling()`, which models the random size; with `strata`, the strata too, which change only the degrees of freedom \tab identical \cr
#'   [design_certainty()] \tab the certainty rows as one unit in a stratum of their own, taken whole \tab as `rest` \cr
#'   [design_weighted()], `"systematic"` \tab `pps = "brewer"` with the inclusion probabilities; with `strata`, `strata` from the strata columns and Brewer's within them \tab close for totals, less so for means \cr
#'   [design_multistage()] with `size` \tab two stages, `ids` the `.cluster` numbers then the rows, `pps = "brewer"` with each cluster's stage-one probability and `m / N_i` at the second stage; each cluster taken with certainty a stratum of its own, whose rows are its units \tab close \cr
#'   [design_spread()] \tab as its first-order design: `ids = ~1`, or `pps = "brewer"` with `size`; with `strata`, `strata` from the strata columns and `fpc` each stratum's size, or Brewer's within them \tab **`survey`'s is larger** \cr
#'   [design_systematic()] with `n_starts` of 2 or more \tab `ids` a `.cluster` column numbering the replicates, `fpc` the `n_starts * interval` possible starts \tab identical \cr
#'   [design_systematic()] \tab `ids = ~1` with `fpc` the rows the walk can reach \tab **differ** \cr
#'   [design_stratified()] with `within` \tab `strata` from the strata columns, `fpc` each stratum's rows walked \tab **differ** \cr
#' }
#'
#' "Identical" means to floating point, and is checked by this package's tests
#' against `survey::svytotal()`, as are domain estimates against
#' `survey::svyby()` and degrees of freedom against `survey::degf()`. The
#' exceptions are real and worth knowing:
#'
#' * **Systematic PPS.** Neither package has its joint probabilities. This one
#'   uses Deville's approximation and `survey` uses Brewer's. In simulations
#'   on skewed frames the standard errors of totals agreed within about 1
#'   percent, and usually within 0.3 percent. Those of means differed more,
#'   typically by 1 to 4 percent and by up to 7 percent at a sampling fraction
#'   of one half. On a frame whose units below certainty had probabilities
#'   from 0.05 to 0.95, with about 17 of them in each sample, totals differed
#'   by 2 percent typically and by up to 8 percent. The degrees of freedom
#'   agree. With `strata`, each stratum is handed over the same way.
#'   [design_multistage()] with `size` meets the same two approximations at
#'   the first stage, while the variance within the selected clusters is the
#'   same in both packages. Where the clusters below certainty had
#'   stage-one probabilities of 0.12 to 0.41, the standard errors of totals
#'   agreed within 0.3 percent over 200 samples. With a cluster at 0.93 among
#'   them, `survey`'s came out 0.97 to 1.48 times this package's: Brewer's
#'   approximation centres on the plain mean of the clusters' totals over
#'   their probabilities, which that cluster pulls, where Deville's weights
#'   each by `1 - pi`.
#' * **Spread.** `survey` has no estimator for a spatially balanced sample, so
#'   it is told only the first-order design and computes a variance that
#'   ignores the spreading. That errs conservative, often substantially.
#' * **Systematic.** `survey`, told only `ids = ~1`, computes the simple-random
#'   variance; [ht_total()] uses the successive-difference approximation, which
#'   can see a trend along the sort order. Neither is design-unbiased, because
#'   no such estimator exists for one systematic sample. Within strata,
#'   `survey` computes the stratified simple-random variance and [ht_total()]
#'   adds up the strata's successive-difference approximations. With
#'   `n_starts` of 2 or more the replicates are clusters and the two agree.
#'   A replicated sample with a start that drew no rows, which happens only
#'   when `n_starts * interval` exceeds the rows walked, is refused, because
#'   `survey` would count one replicate fewer.
#'
#' Units taken with certainty are in every sample and add nothing to the
#' variance: the rows above a [design_certainty()] threshold, strata and
#' intervals taken whole, and the units [design_weighted()] and
#' [design_spread()] with `size` put at probability 1. Each group of them goes
#' to `survey` as one primary sampling unit in a stratum of its own, with an
#' `fpc` of 1, so `survey` gives them no variance and `survey::degf()` leaves
#' them out of the degrees of freedom, as [ht_total()] does. For the PPS
#' designs this also keeps them out of Brewer's approximation, which left
#' among the other units inflated `survey`'s standard error, by a factor of 23
#' on a frame with four dominant units.
#'
#' A sample in which every row was taken with certainty, a census, goes to
#' `survey` as `ids = ~1` with `fpc` the number of rows. Both packages then
#' report a variance of 0, and `survey::degf()` gives `n - 1` where
#' [ht_total()] gives `Inf`. A sample with no rows is refused, and so is a
#' certainty sample whose `rest` drew none of the rows below the threshold,
#' whose total `survey` would report as exact where [ht_total()] reports no
#' variance.
#'
#' Two compositions have no single `survey` design and are refused rather than
#' approximated: [design_certainty()] over a cluster, multistage or Poisson
#' `rest`, or a systematic `rest` with several starts, where the certainty rows
#' and the rest are different kinds of sampling unit. [design_bootstrap()] is
#' refused outright — it resamples the sample, so there is no finite
#' population for `svydesign()` to represent.
#'
#' A stratum or interval with a single sampled row makes `survey` stop with
#' "Stratum has only one PSU" unless `options(survey.lonely.psu)` says
#' otherwise; [ht_total()] declines in the same situation. Draw at least two
#' per stratum.
#'
#' @section Weights and design columns:
#' `survey` gets the sample's `.weight` column as its weights. Weights adjusted
#' after the draw, for nonresponse or by trimming, and written into `.weight`
#' are therefore estimated here, while [ht_total()] and [ht_mean()] keep to the
#' design probabilities in `.prob` and warn that `.weight` differs. `.prob`
#' itself must stay as [draw()] set it: it supplies the probabilities of the
#' PPS mappings above, and a sample whose `.prob` was changed is refused.
#'
#' Strata, intervals, cluster labels and population counts come from the copy
#' of the frame stored with the sample, not from the sample's own columns, and
#' are passed in the columns `.stratum`, `.cluster`, `.psu`, `.fpc`, `.fpc2`
#' and `.unit`. Renaming or recoding a design column in the sample after the draw
#' does not change the survey design.
#'
#' @param sample A data frame returned by [draw()] with `weights = TRUE`. See
#'   "Weights and design columns".
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
#' survey::svyquantile(~spend, des, quantiles = 0.5)
#'
#' @seealso [ht_total()], [ht_mean()]
#' @export
as_svydesign <- function(sample, ...) {
  require_suggested("survey", "handing a sample to the survey package")
  state <- drawn_state(sample)
  design <- state$design
  pop <- state$pop
  if (design_type(design) == "bootstrap") {
    stop("A bootstrap design has no survey equivalent: it resamples the ",
         "sample rather than drawing a probability sample of a finite ",
         "population.", call. = FALSE)
  }

  if (!nrow(sample)) {
    stop_not_drawn(paste0(
      "has no rows, so there is no design to hand to `survey`. For an empty ",
      "draw from a design of random size, ht_total() gives the total, 0, ",
      "without a variance."))
  }

  # `survey` keeps the data it is given by reference, so it gets its own copy.
  # An in-place reorder of the sample afterwards would otherwise reach it.
  dat <- snapshot_frame(as.data.frame(sample))
  parts <- survey_parts(design, dat, state$units, pop, rows = state$rows)
  if (all(is_certain(dat$.prob))) parts <- census_parts(dat)
  dat$.fpc <- parts$fpc
  # Clusters come from the frame, so a cluster column renamed or recoded in
  # the sample after the draw still gives `survey` the clusters that were
  # drawn.
  if (!is.null(parts$cluster)) dat$.cluster <- parts$cluster
  fpc <- stats::as.formula("~.fpc")
  if (!is.null(parts$fpc2)) {
    # Two stages: clusters, then rows within them, each with its own
    # population count, so `survey` computes the full two-stage variance.
    dat$.fpc2 <- parts$fpc2
    dat$.unit <- seq_len(nrow(dat))
    fpc <- stats::as.formula("~.fpc + .fpc2")
  }
  if (!is.null(parts$stratum)) dat$.stratum <- parts$stratum
  ids <- parts$ids
  if (any(parts$whole)) {
    # Units taken with certainty share one primary sampling unit per stratum,
    # with an `fpc` of one unit, so `survey` treats each group as taken whole.
    dat$.psu <- shared_psu(parts$stratum, parts$whole)
    ids <- stats::as.formula("~.psu")
  }
  args <- c(list(data = dat, weights = stats::as.formula("~.weight"),
                 ids = ids, fpc = fpc),
            parts$extra)
  if (!is.null(parts$stratum)) {
    args$strata <- stats::as.formula("~.stratum")
  }
  do.call(survey::svydesign, utils::modifyList(args, list(...)))
}

#' Refuse a certainty design whose `rest` has no place in one survey design
#' @noRd
refuse_nested <- function(prefix, what, why) {
  if (nzchar(prefix)) {
    stop("A certainty design over ", what, " has no single survey equivalent: ",
         why, ". Estimate it with ht_total() or ht_mean() instead.",
         call. = FALSE)
  }
  invisible(NULL)
}

#' Strata, population counts and cluster ids for one design
#'
#' `survey` reconstructs a variance from the design's shape, not from the
#' inclusion probabilities, so each design has to be expressed in its own terms:
#' a stratified design as strata, a cluster design as ids, a temporal design as
#' strata over its intervals, and a certainty design as a stratum that was taken
#' whole.
#'
#' Units taken with certainty -- a certainty design's rows above the threshold,
#' strata and intervals taken whole, and the units a size-based design puts at
#' probability 1 -- are flagged in `whole`. Each stratum of them becomes one
#' primary sampling unit with an `fpc` of 1, which reads the same as a count of
#' units and as a probability. `survey` gives such a stratum no variance, and
#' `survey::degf()` no degree of freedom, as `design_df()` does.
#'
#' `dat` is the sample handed to `survey`; `units` is the frame's copy of the
#' same rows, in the same order, and is where strata and intervals are read.
#' Clusters are read from the frame at `rows` and returned as integer codes in
#' `cluster`, for a column of their own in `dat`.
#'
#' @noRd
survey_parts <- function(design, dat, units, pop, rows, prefix = "") {
  one <- stats::as.formula("~1")
  # `code` holds each frame row's group as a position in `labels`; the
  # sampled units are frame rows `rows`.
  by_group <- function(code, labels) {
    sizes <- tabulate(code, length(labels))
    unit <- code[rows]
    whole <- taken_whole(unit, dat$.prob)
    fpc <- sizes[unit]
    fpc[whole] <- 1L
    list(stratum = paste0(prefix, labels[unit]), fpc = fpc, ids = one,
         whole = whole)
  }

  switch(design_type(design),
    stratified = {
      if (walks_strata(design)) {
        # Rows the walk drops for a missing `order_by` are in no stratum.
        al <- stratum_allocation(design, pop)
        return(by_group(al$code, unique_labels(
          structure(integer(0), levels = al$labels, class = "factor"))))
      }
      key <- group_key(pop, design$strata)
      by_group(as.integer(key), unique_labels(key))
    },
    temporal = {
      bucket <- temporal_bucket(design, pop)
      labels <- unique(bucket[!is.na(bucket)])
      by_group(match(bucket, labels), labels)
    },
    weighted = {
      if (!is.null(design$strata)) {
        return(weighted_strata_parts(design, dat, pop, prefix, by_group))
      }
      if (design$method != "poisson") {
        # Neither package has the joint probabilities of systematic PPS.
        # `survey` uses Brewer's approximation and this package Deville's.
        return(pps_parts(dat, prefix, "brewer"))
      }
      # Poisson sampling has a random size, so no finite population correction
      # over a fixed n describes it. `survey` models it directly.
      refuse_nested(prefix, "a Poisson `rest`",
                    "its sample size is random, which no stratum can express")
      pps_parts(dat, prefix, survey::poisson_sampling(dat$.prob))
    },
    cluster = ,
    multistage = {
      refuse_nested(prefix, "a cluster or multistage `rest`",
                    paste0("the certainty rows and the clusters are different ",
                           "kinds of sampling unit"))
      # Integer codes, not labels: `survey` makes a factor of the ids, which
      # compares them as text and merges two numbers that print alike.
      if (!is.null(design$size)) {
        return(pps_stage_parts(design, units, pop))
      }
      cl <- count_clusters(design, pop)
      n_clusters <- cl$total
      cluster <- cl$code[rows]
      if (design_type(design) == "cluster") {
        return(list(stratum = NULL, fpc = rep(n_clusters, nrow(dat)),
                    ids = stats::as.formula("~.cluster"), cluster = cluster))
      }
      sizes <- tabulate(cl$code, n_clusters)
      list(stratum = NULL, fpc = rep(n_clusters, nrow(dat)),
           fpc2 = sizes[cluster],
           ids = stats::as.formula("~.cluster + .unit"), cluster = cluster)
    },
    reservoir = {
      # Rows past `max_items` are never read, so the population is the rows
      # the stream reached, not the whole frame.
      list(stratum = NULL, fpc = rep(reservoir_reach(design, nrow(pop)),
                                     nrow(dat)), ids = one)
    },
    systematic = {
      if (replicated(design)) {
        refuse_nested(prefix, "a systematic `rest` with several starts",
                      paste0("the certainty rows and the replicates are ",
                             "different kinds of sampling unit"))
        return(replicate_parts(design, dat, pop, rows))
      }
      # Rows `na_rm` drops from the walk can never be drawn, so the
      # population is the rows with a sort key, not the whole frame.
      list(stratum = NULL, fpc = rep(reachable_rows(design, pop), nrow(dat)),
           ids = one)
    },
    spread = {
      # Neither package has the pivotal method's joint probabilities. Handed
      # to `survey` as the unequal- or equal-probability draw it is at first
      # order, which does not credit the spreading and so errs conservative.
      fr <- spread_frame(design, pop)
      if (!is.null(design$strata)) {
        return(spread_strata_parts(design, dat, nrow(pop), fr, prefix,
                                   by_group))
      }
      if (is.null(design$size)) {
        list(stratum = NULL, fpc = rep(length(fr$keep), nrow(dat)), ids = one)
      } else {
        pps_parts(dat, prefix, "brewer")
      }
    },
    spatial = {
      # The design's population is the rows inside the region, not the frame.
      # Handing `survey` the frame size overstated its standard error by 6-11%.
      inside <- sum(spatial_inside(design, pop))
      list(stratum = NULL, fpc = rep(inside, nrow(dat)), ids = one)
    },
    certainty = {
      sp <- certainty_split(design, pop)
      below <- sp$data[sp$rest, , drop = FALSE]
      # Membership of the certainty part, not `.prob == 1`: `rest` can put a
      # row at probability 1 of its own accord -- a dominant PPS unit, or an
      # `n` equal to the rows below the threshold -- and such rows belong to
      # `rest`'s strata, which flag them in `whole` themselves.
      above <- rows %in% sp$keep[sp$take]
      if (!any(!above) && length(sp$rest)) {
        stop_not_drawn(paste0(
          "holds none of the ", rows_below(design, length(sp$rest)),
          ", because `rest` drew none of them, and `survey` would report the ",
          "total as exact. ht_total() and ht_mean() report no variance for ",
          "such a sample; use a `rest` with a fixed sample size to estimate ",
          "one."))
      }

      out <- list(stratum = character(nrow(dat)), fpc = numeric(nrow(dat)),
                  ids = one, whole = above)
      out$stratum[above] <- paste0(prefix, "(certain)")
      out$fpc[above] <- 1
      if (any(!above)) {
        rest <- survey_parts(design$rest, dat[!above, , drop = FALSE],
                             units[!above, , drop = FALSE], below,
                             rows = match(rows[!above], sp$keep[sp$rest]),
                             prefix = paste0(prefix, "(rest) "))
        out$stratum[!above] <- rest$stratum %||% paste0(prefix, "(rest)")
        out$fpc[!above] <- rest$fpc
        out$whole[!above] <- rest$whole %||% FALSE
        out$extra <- rest$extra
      }
      out
    },
    list(stratum = NULL, fpc = rep(nrow(pop), nrow(dat)), ids = one)
  )
}

#' A multistage sample with clusters selected in proportion to size
#'
#' Two stages, with Brewer's approximation at the first. Clusters below
#' certainty share the stratum `(pps)`, with their stage-one probabilities
#' as `fpc` and `m / N_i` at the second stage, so `survey` adds the
#' variance within each cluster divided by its probability, as
#' `pps_stage_variance()` does. A cluster taken with certainty is a stratum
#' of its own whose units are its sampled rows, with `m / N_i` at the first
#' stage and nothing left to sample at the second; one taken whole as well
#' is a single unit with an `fpc` of 1. `survey::degf()` then counts what
#' `design_df()` counts.
#'
#' @noRd
pps_stage_parts <- function(design, units, pop) {
  cu <- cluster_units(design, units, pop)
  st <- cluster_stage(design, pop)
  p1 <- st$pi1[cu$code]
  share <- st$m / st$sizes[cu$code]
  certain <- is_certain(p1)
  list(stratum = ifelse(certain, paste0("(certain) ", cu$code), "(pps)"),
       fpc = ifelse(certain, pmin(share, 1), p1),
       fpc2 = ifelse(certain, 1, share),
       ids = stats::as.formula("~.cluster + .unit"), cluster = cu$psu,
       extra = list(pps = "brewer"))
}

#' A systematic sample with several starts, as a sample of clusters
#'
#' The `m` starts are a simple random sample of the `m * interval` possible
#' ones, and each replicate is the cluster of rows its start takes. Handed
#' to `survey` as clusters numbered in `.cluster`, with that count of
#' possible starts as `fpc`, `survey`'s one-stage variance is the random
#' group variance of `replicate_variance()`. A replicate whose start fell
#' past the last row walked has no rows to pass, and `survey` would count
#' one replicate fewer, so such a sample is refused.
#'
#' @noRd
replicate_parts <- function(design, dat, pop, rows) {
  m <- systematic_starts(design)
  walk <- systematic_positions(design, pop)
  cluster <- systematic_replicates(design, walk, rows)
  got <- length(unique(cluster))
  # An interval of 1 takes every row, and a census needs no replicates.
  if (got < m && !all(is_certain(dat$.prob))) {
    stop("A systematic sample with an empty replicate has no `survey` ",
         "equivalent: ", m - got, " of its ", m, " starts fell past the last ",
         "of the ", sum(!is.na(walk)), " rows walked and drew none, and ",
         "`survey` would count only the ", got, " with rows. Estimate it ",
         "with ht_total() or ht_mean(), which count an empty replicate as a ",
         "total of 0, or draw with `n_starts * interval` of at most the rows ",
         "walked.", call. = FALSE)
  }
  list(stratum = NULL, fpc = rep(m * design$interval, nrow(dat)),
       ids = stats::as.formula("~.cluster"), cluster = cluster)
}

#' A spread design with strata, as stratified sampling or stratified PPS
#'
#' Each stratum is handed over as the draw it is at first order: without
#' `size`, a simple random sample of its rows, and with it, Brewer's
#' approximation over their probabilities, as `strata_pps_parts()` gives
#' it. `by_group` is the helper of `survey_parts()`, which reads the sampled
#' frame rows, and `n_pop` the number of rows in the frame.
#'
#' @noRd
spread_strata_parts <- function(design, dat, n_pop, fr, prefix, by_group) {
  # Each frame row's stratum; rows the design drops count in none.
  code <- rep(NA_integer_, n_pop)
  code[fr$keep] <- fr$code
  labels <- unique_labels(structure(integer(0), levels = fr$labels,
                                    class = "factor"))
  if (is.null(design$size)) return(by_group(code, labels))
  strata_pps_parts(code, labels, dat, prefix, by_group, "brewer")
}

#' A weighted design with strata, as stratified PPS
#'
#' Systematic PPS goes to `survey` as Brewer's approximation within each
#' stratum, and Poisson sampling as `survey::poisson_sampling()` with the
#' strata, whose variance does not depend on them but whose degrees of
#' freedom count them, as `design_df()` does.
#'
#' @noRd
weighted_strata_parts <- function(design, dat, pop, prefix, by_group) {
  pps <- if (design$method == "poisson") {
    refuse_nested(prefix, "a Poisson `rest`",
                  "its sample size is random, which no stratum can express")
    survey::poisson_sampling(dat$.prob)
  } else {
    "brewer"
  }
  st <- weighted_setup(design, pop)
  # Each frame row's stratum; rows the design drops count in none.
  code <- rep(NA_integer_, nrow(pop))
  code[st$keep] <- st$code
  labels <- unique_labels(structure(integer(0), levels = st$labels,
                                    class = "factor"))
  strata_pps_parts(code, labels, dat, prefix, by_group, pps)
}

#' Strata sampled in proportion to size, with a stratum for certain units
#'
#' `code` gives each frame row's stratum as a position in `labels`. Rows
#' taken with certainty go in a stratum of their own, `(certain)`, as
#' `pps_parts()` puts them, so that `survey`'s approximation runs over the
#' other rows of each stratum. Every other row has its inclusion probability
#' as `fpc`, for `pps`.
#'
#' @noRd
strata_pps_parts <- function(code, labels, dat, prefix, by_group, pps) {
  out <- by_group(code, labels)
  certain <- is_certain(dat$.prob)
  # A stratum of the frame could be called `(certain)` too.
  cert <- paste0(prefix, "(certain)")
  while (cert %in% paste0(prefix, labels)) cert <- paste0(cert, " ")
  out$stratum[certain] <- cert
  out$fpc <- ifelse(certain, 1, dat$.prob)
  out$whole <- certain
  out$extra <- list(pps = pps)
  out
}

#' Which rows lie in a stratum taken whole
#'
#' Rows of a stratum or interval are drawn with one probability, so a stratum
#' is taken whole when every sampled row in it is certain.
#'
#' @noRd
taken_whole <- function(stratum, prob) {
  certain <- is_certain(prob)
  stratum %in% stratum[certain] & !(stratum %in% stratum[!certain])
}

#' A size-based design, with its certain units in a stratum of their own
#'
#' `survey`'s Brewer approximation and its Poisson variance read `fpc` as
#' probabilities. A unit at probability 1 left among the others still enters
#' Brewer's centring on the mean of `y / pi` and its `n / (n - 1)` scaling,
#' which put `survey`'s standard error at 23 times this package's on one frame
#' with four dominant units. In a stratum of their own, taken whole, they add
#' nothing and Brewer runs over the units below certainty.
#'
#' @noRd
pps_parts <- function(dat, prefix, pps) {
  certain <- is_certain(dat$.prob)
  list(stratum = if (any(certain)) {
         paste0(prefix, ifelse(certain, "(certain)", "(pps)"))
       },
       fpc = ifelse(certain, 1, dat$.prob), ids = stats::as.formula("~1"),
       whole = certain, extra = list(pps = pps))
}

#' One primary sampling unit for each stratum of units taken with certainty
#'
#' Every other row is its own unit. The ids are unique across strata, so the
#' units nest in their strata without `nest = TRUE`.
#'
#' @noRd
shared_psu <- function(stratum, whole) {
  psu <- seq_along(whole)
  psu[whole] <- length(whole) + match(stratum[whole], unique(stratum[whole]))
  psu
}

#' A sample whose every row was taken with certainty
#'
#' One stratum, every row its own unit, and `fpc` the number of rows, so that
#' `survey`'s finite population correction is 0. `svydesign()` reads an `fpc`
#' made only of 1s as neither counts nor probabilities, and stopped with
#' "subscript out of bounds" on the probabilities of 1 the PPS mappings gave
#' a census.
#'
#' @noRd
census_parts <- function(dat) {
  n <- nrow(dat)
  list(stratum = NULL, fpc = rep(n, n), ids = stats::as.formula("~1"))
}
