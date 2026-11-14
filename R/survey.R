#' Hand a sample to the survey package
#'
#' Builds a [survey::svydesign()] object from a drawn sample, so the analysis
#' this package does not do — regression, calibration, analysis with
#' replicate weights — can be done by the package that does.
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
#'   [design_multistage()], and [design_cluster()] with `balanced = TRUE` \tab two stages, `ids` the `.cluster` numbers then the rows, `fpc` the number of clusters and each cluster's size; with one sampled row from a selected cluster of several, one stage, `ids` the `.cluster` numbers and no `fpc` \tab identical \cr
#'   [design_weighted()], `"poisson"` \tab `survey::poisson_sampling()`, which models the random size; with `strata`, the strata too, which change only the degrees of freedom \tab identical \cr
#'   [design_certainty()] \tab the certainty rows as one unit in a stratum of their own, taken whole, beside `rest` as it maps alone; over a Poisson `rest`, `survey::poisson_sampling()` from the whole sample's `.prob` \tab as `rest` \cr
#'   [design_weighted()], `"systematic"` \tab `pps = "brewer"` with the inclusion probabilities; with `strata`, `strata` from the strata columns and Brewer's within them \tab close for totals, less so for means \cr
#'   [design_multistage()] with `size` \tab two stages, `ids` the `.cluster` numbers then the rows, `pps = "brewer"` with each cluster's stage-one probability and `m / N_i` at the second stage; each cluster taken with certainty a stratum of its own, whose rows are its units; with one sampled row from a selected cluster of several, one stage, `ids` the `.cluster` numbers and no `fpc` \tab close; identical in one stage \cr
#'   [design_spread()] \tab as its first-order design: `ids = ~1`, or `pps = "brewer"` with `size`; with `strata`, `strata` from the strata columns and `fpc` each stratum's size, or Brewer's within them \tab **differ** \cr
#'   [design_systematic()] with `n_starts` of 2 or more \tab `ids` a `.cluster` column numbering the replicates, `fpc` the `n_starts * interval` possible starts \tab identical \cr
#'   [design_systematic()] \tab `ids = ~1` with `fpc` the rows the walk can reach \tab **differ** \cr
#'   [design_stratified()] with `within` \tab `strata` from the strata columns, `fpc` each stratum's rows walked \tab **differ** \cr
#' }
#'
#' "Identical" means to floating point under `survey`'s default options, and
#' is checked by this package's tests against `survey::svytotal()`, as are
#' domain estimates against `survey::svyby()` and degrees of freedom against
#' `survey::degf()`. The exceptions are real and worth knowing:
#'
#' * **Systematic PPS.** Neither package has its joint probabilities. This one
#'   uses Deville's approximation and `survey` uses Brewer's. In simulations
#'   on frames of 100 to 300 rows with gamma and lognormal sizes, the
#'   standard errors of totals usually agreed within 0.5 percent, and within
#'   3 percent in every sample. Those of means differed more, and more as the
#'   sampling fraction grew: typically by under 2 percent at a fraction of a
#'   tenth or a fifth, 2 to 4 percent at a third and 3 to 7 percent at a
#'   half, where single samples differed by up to 12 percent. On a frame
#'   whose units below certainty had probabilities from 0.05 to 0.95, with
#'   about 17 of them in each sample, totals differed by 2 percent typically
#'   and by up to 8 percent. The degrees of freedom agree. With `strata`,
#'   each stratum is handed over the same way.
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
#'   ignores the spreading. Relative to the true variance that errs
#'   conservative, often substantially. Relative to [ht_total()]'s local mean
#'   approximation it goes either way. In simulations drawing 30 of 300
#'   points spread over a square, `survey`'s standard error was the larger in
#'   every sample when `y` varied smoothly over the square, three times as
#'   large at the median. When `y` had no spatial pattern it was the smaller
#'   in 67 to 79 percent of samples, by 4 to 6 percent at the median, and
#'   drawing 150 of the 300 it was the smaller in every sample, by 30 percent
#'   at the median, because the local mean has no finite population
#'   correction.
#' * **Systematic.** `survey`, told only `ids = ~1`, computes the simple-random
#'   variance; [ht_total()] uses the successive-difference approximation, which
#'   can see a trend along the sort order. Neither is design-unbiased, and at
#'   a step of 2 rows or more no such estimator exists for one systematic
#'   sample. Within strata, `survey` computes the stratified simple-random
#'   variance and [ht_total()] adds up the strata's successive-difference
#'   approximations. With `n_starts` of 2 or more the replicates are clusters
#'   and the two agree. A replicated sample with a start that drew no rows,
#'   which happens only when `n_starts * interval` exceeds the rows walked, is
#'   refused, because `survey` would count one replicate fewer.
#'
#' Units taken with certainty are in every sample and add nothing to the
#' variance: the rows above a [design_certainty()] threshold, strata and
#' intervals taken whole, the units [design_weighted()] and
#' [design_spread()] with `size` put at probability 1, and the clusters
#' [design_multistage()] and [design_cluster()] with `balanced = TRUE` take
#' whole when they select every cluster. Each group of them goes to `survey`
#' as one primary sampling unit in a stratum of its own, with an `fpc` of 1,
#' so `survey` gives them no variance and `survey::degf()` leaves them out
#' of the degrees of freedom, as [ht_total()] does. For the PPS designs this
#' also keeps them out of Brewer's approximation, which left among the other
#' units inflated `survey`'s standard error, by a factor of 23 on a frame
#' with four dominant units. Beside a `rest` of clusters, or of the
#' replicates of several systematic starts, the certainty rows form one more
#' primary sampling unit. In two stages, such a unit is taken whole at the
#' second stage too.
#'
#' A sample in which every row was taken with certainty, a census, goes to
#' `survey` as `ids = ~1` with `fpc` the number of rows. Both packages then
#' report a variance of 0, and `survey::degf()` gives `n - 1` where
#' [ht_total()] gives `Inf`. A sample with no rows is refused, and so is a
#' certainty sample whose `rest` drew none of the rows below the threshold,
#' whose total `survey` would report as exact where [ht_total()] reports no
#' variance.
#'
#' [design_certainty()] over a Poisson `rest` goes to `survey` as
#' `survey::poisson_sampling()` built from the whole sample's `.prob`, so the
#' certainty rows enter at probability 1 and add no variance.
#' [design_bootstrap()] is refused outright — it resamples the sample, so
#' there is no finite population for `svydesign()` to represent.
#'
#' A stratum or interval with a single sampled row makes `survey` stop with
#' "Stratum has only one PSU" unless `options(survey.lonely.psu)` says
#' otherwise; [ht_total()] declines in the same situation. Draw at least two
#' per stratum. Under [design_multistage()] with `size`, a cluster taken with
#' certainty is a stratum of its sampled rows, so one that gives a single row
#' of several is such a stratum.
#' `options(survey.ultimate.cluster = TRUE)` makes `survey` leave out the
#' second stage of every two-stage design, so its standard errors then
#' differ from [ht_total()]'s. The agreement the table states holds under
#' the default, `FALSE`.
#'
#' A sample weighted conditional on the clusters selected (see "Design
#' weights" in [draw()]) goes over in two stages too. `survey` works out
#' each cluster's second-stage sampling fraction from the rows sampled in
#' it, which is the take given that selection, so its estimate and standard
#' error are those of [ht_total()].
#'
#' One sampled row from a selected cluster of several leaves the variation
#' within clusters unmeasured. [ht_total()] then falls back to the
#' delete-a-cluster jackknife, which is the ultimate-cluster variance with no
#' finite population correction, and `survey` gets that design: one stage,
#' `ids` the `.cluster` numbers and no `fpc`. So does a certainty design's
#' `rest` of balanced clusters, which [draw()] weights without regard to the
#' clusters selected, leaving no two-stage estimator. Under
#' [design_multistage()] with `size` the jackknife runs over the clusters
#' below certainty, and a cluster taken whole is held fixed in a stratum of
#' its own in both packages.
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
#' does not change the survey design. Those six names are reserved: a sample
#' with a column of one of them is refused rather than overwritten. Rename
#' the column in the sample, which leaves the design as drawn.
#'
#' @section Replicate weights and calibration:
#' `survey::as.svrepdesign()` turns the design returned here into one with
#' replicate weights. For samples from [design_simple()],
#' [design_stratified()] without `within`, [design_cluster()] without
#' `balanced`, [design_temporal()] and [design_spatial()], its jackknife
#' gives the standard errors of `ht_total(variance = "jackknife")`. It stops
#' with an error on the designs handed over with Brewer's approximation or
#' as Poisson sampling, those of [design_weighted()], [design_spread()] and
#' [design_multistage()] with `size`, and [design_certainty()] over them,
#' because `survey` builds no replicates for its own
#' probability-proportional-to-size and Poisson designs. That is a
#' limitation of `survey`. For a systematic sample, alone or within strata,
#' it builds a jackknife that ignores the order of the walk. For
#' [design_cluster()] with `balanced = TRUE` it keeps the finite population
#' correction of the first stage, which `ht_total(variance = "jackknife")`
#' leaves out, and warns that it dropped the corrections after the first
#' stage. Its variance is then the jackknife's times `1 - n / N`, for `n`
#' of the `N` clusters selected.
#'
#' [replicate_weights()] builds replicate weights from the sample itself,
#' for [survey::svrepdesign()]: the jackknife of
#' `ht_total(variance = "jackknife")`, for balanced clusters and the designs
#' handed over with Brewer's approximation too, and successive differences
#' along the walk for systematic samples. Poisson sampling, which the
#' jackknife understates, has none, and `survey` estimates it from the
#' design returned here.
#'
#' [frame_totals()] gives the population totals that `survey::calibrate()`
#' and `survey::postStratify()` take, from the frame stored with the
#' sample, over the rows the design could reach. Calibration is left to
#' `survey`, on the design returned here.
#'
#' @param sample A data frame returned by [draw()] with `weights = TRUE`. See
#'   "Weights and design columns".
#' @param ... Passed to [survey::svydesign()] by name, and an unnamed
#'   argument is refused. Anything named here overrides what the mapping
#'   above would have supplied, so `nest = TRUE` or a replacement `fpc` is
#'   yours to set.
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
#' coef(summary(survey::svyglm(spend ~ site, des)))
#'
#' @seealso [ht_total()], [ht_mean()], [replicate_weights()],
#'   [frame_totals()]
#' @export
as_svydesign <- function(sample, ...) {
  require_suggested("survey", "handing a sample to the survey package")
  overrides <- list(...)
  check_named_dots(overrides)
  # drawn_state() refuses a bootstrap design, which has no survey equivalent.
  state <- drawn_state(sample)
  design <- state$design
  pop <- state$pop

  if (!nrow(sample)) {
    stop_not_drawn(paste0(
      "has no rows, so there is no design to hand to `survey`. For an empty ",
      "draw from a design of random size, ht_total() gives the total, 0, ",
      "without a variance."))
  }
  check_survey_columns(sample)

  # `survey` keeps the data it is given by reference, so it gets its own copy.
  # An in-place reorder of the sample afterwards would otherwise reach it.
  dat <- snapshot_frame(as.data.frame(sample))
  parts <- survey_parts(design, dat, state$units, pop, rows = state$rows)
  if (all(is_certain(dat$.prob))) parts <- census_parts(dat)
  # Clusters come from the frame, so a cluster column renamed or recoded in
  # the sample after the draw still gives `survey` the clusters that were
  # drawn.
  if (!is.null(parts$cluster)) dat$.cluster <- parts$cluster
  fpc <- NULL
  if (!is.null(parts$fpc)) {
    dat$.fpc <- parts$fpc
    fpc <- stats::as.formula("~.fpc")
  }
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
    # The other units keep their own ids, which are rows, or the clusters or
    # replicates numbered in `.cluster`.
    dat$.psu <- shared_psu(parts$stratum, parts$whole, parts$cluster)
    ids <- if (is.null(parts$fpc2)) {
      stats::as.formula("~.psu")
    } else {
      stats::as.formula("~.psu + .unit")
    }
  }
  if (identical(parts$extra$pps, "poisson")) {
    # Built from the whole sample, so rows a certainty design took above its
    # threshold enter at probability 1 beside its Poisson `rest`.
    parts$extra$pps <- survey::poisson_sampling(dat$.prob)
  }
  args <- c(list(data = dat, weights = stats::as.formula("~.weight"),
                 ids = ids, fpc = fpc),
            parts$extra)
  if (!is.null(parts$stratum)) {
    args$strata <- stats::as.formula("~.stratum")
  }
  do.call(survey::svydesign, utils::modifyList(args, overrides))
}

#' Replicate weights for a sample
#'
#' Builds replicate weights from a drawn sample for [survey::svrepdesign()],
#' from the design the sample was drawn with. Each replicate is the
#' sample's `.weight` column times a factor, and the variance of a total
#' over the replicates is that of one of [ht_total()]'s estimators:
#'
#' * For most designs, the delete-a-unit jackknife of
#'   `ht_total(variance = "jackknife")`. Each replicate deletes one primary
#'   sampling unit, a cluster where the design has them and a row
#'   otherwise, and reweights the other units of its stratum; strata are
#'   the design's strata, the intervals of [design_temporal()] and the
#'   strata of [design_weighted()] and [design_spread()]. `rscales` holds
#'   `(m - 1) / m` for a stratum of `m` sampled units, times its finite
#'   population correction.
#'   Units taken with certainty keep their weight in every replicate, as
#'   the jackknife holds them fixed. Where the default of [ht_total()] is
#'   another estimator, such as Deville's approximation for
#'   [design_weighted()] with `method = "systematic"`, the local mean
#'   approximation for [design_spread()] or the two-stage variance for
#'   [design_multistage()] and balanced [design_cluster()], the replicates
#'   give the jackknife's variance, not that default.
#' * For [design_systematic()] with one start, and [design_stratified()]
#'   with `within`, successive-difference replicates in the order the
#'   design walked, whose variance is the successive-difference
#'   approximation [ht_total()] gives by default, with its term for a
#'   random sample size. Each pair of neighbours in the walk shares a row of
#'   a Hadamard matrix, which raises one of the two and lowers the other.
#'   The first and last rows of a walk are not paired with each other: the
#'   circular form of successive-difference replication pairs them, which
#'   adds the squared difference between the two ends of the sort order. On
#'   a frame of 240 rows sorted on a size measure related to `y`, at an
#'   interval of 8, that more than doubled the variance in each of the 8
#'   possible samples.
#' * For [design_systematic()] with `n_starts` of 2 or more, one replicate
#'   per start, each deleting the rows its start drew. The variance is the
#'   random group variance [ht_total()] gives by default.
#'
#' Rows above the threshold of [design_certainty()] keep their weight in
#' every replicate, and the replicates for the rows below it are those of
#' `rest`. A sample whose every row was taken with certainty gets one
#' replicate equal to its weights, and a variance of 0.
#'
#' Poisson sampling, alone or as the `rest` of [design_certainty()], is
#' refused: its sample size is random, and deleting rows from it
#' understates the variance. [as_svydesign()] hands it to `survey`, which
#' models it. A sample the jackknife declines, such as one with a single
#' unit in a stratum, is refused with the jackknife's reason.
#'
#' Every replicate is a multiple of `.weight`, so an adjustment written into
#' `.weight` after the draw, such as trimming, is carried into each one
#' unchanged. One that depends on which rows responded is more accurately
#' made again within each replicate.
#'
#' The matrix of weights takes `8 * n * R` bytes for `n` sampled rows and
#' `R` replicates, and building it about twice that. The jackknife has a
#' replicate for each sampled unit, and successive differences have a
#' power of 2 above the number of rows, so the size grows with the square
#' of the sample. Past about 1 GB, `replicate_weights()` refuses: from
#' about 11,200 rows for a jackknife over rows, and about 8,200 for
#' successive differences. Such a sample can still be estimated with
#' [ht_total()], or handed to `survey` with [as_svydesign()].
#'
#' @param sample A data frame returned by [draw()] with `weights = TRUE`.
#' @param type The form of the jackknife. `"JKn"`, the default, deletes
#'   units within strata and keeps each stratum's correction in `rscales`.
#'   `"JK1"` deletes across the whole sample, with the correction in
#'   `scale`, and is given only where the jackknife has one stratum; with
#'   more, it would count the differences between strata as sampling
#'   variance, and is refused. A systematic sample with one start, and
#'   strata walked in order, get successive-difference replicates of type
#'   `"other"` whichever is given.
#'
#' @return A numeric matrix with a row for each row of `sample`, in its
#'   order, and a column for each replicate, holding replicate weights.
#'   Its attributes `type`, `scale` and `rscales` are the arguments of the
#'   same names of [survey::svrepdesign()], which with
#'   `combined.weights = TRUE` computes the variance as
#'   `scale * sum(rscales * (t_r - t)^2)` over the replicate estimates
#'   `t_r`, with `t` their mean. For a total that mean is the sample's
#'   estimate.
#'
#' @examplesIf requireNamespace("survey", quietly = TRUE)
#' pop <- data.frame(
#'   id = 1:400,
#'   site = rep(c("a", "b", "c", "d"), times = c(200, 100, 60, 40)),
#'   spend = 10 + (1:400 * 37) %% 490
#' )
#' s <- draw(pop, design_stratified("site", n = 60), seed = 1, weights = TRUE)
#'
#' rw <- replicate_weights(s)
#' des <- survey::svrepdesign(
#'   data = s, repweights = rw, weights = ~.weight,
#'   type = attr(rw, "type"), scale = attr(rw, "scale"),
#'   rscales = attr(rw, "rscales"), combined.weights = TRUE
#' )
#' survey::svytotal(~spend, des)
#'
#' # The same standard error as the jackknife here
#' ht_total(s, "spend", variance = "jackknife")
#'
#' # A systematic sample on a sorted frame: the replicates follow the walk
#' sys <- draw(pop, design_systematic(interval = 8, order_by = "spend"),
#'             seed = 1, weights = TRUE)
#' rw <- replicate_weights(sys)
#' des <- survey::svrepdesign(
#'   data = sys, repweights = rw, weights = ~.weight,
#'   type = attr(rw, "type"), scale = attr(rw, "scale"),
#'   rscales = attr(rw, "rscales"), combined.weights = TRUE
#' )
#' survey::SE(survey::svytotal(~spend, des))
#' ht_total(sys, "spend")$se
#'
#' @seealso [as_svydesign()], [ht_total()]
#' @export
replicate_weights <- function(sample, type = c("JKn", "JK1")) {
  type <- check_choice(type, c("JKn", "JK1"), "type")
  # drawn_state() refuses a bootstrap design, which has no replicate weights.
  state <- drawn_state(sample)
  design <- state$design
  if (!nrow(sample)) {
    stop_not_drawn(paste0(
      "has no rows, so there are no weights to replicate. For an empty draw ",
      "from a design of random size, ht_total() gives the total, 0, without ",
      "a variance."))
  }
  if (is.null(state$pi_all)) check_prob(sample)
  w <- check_replicate_base(sample)
  reps <- tryCatch(
    replicate_factors(design, state$units, state$pop, state$rows,
                      sample$.prob, type),
    drawn_declined = function(e) list(note = conditionMessage(e)))
  if (is.null(reps$factors)) {
    stop("`sample` has no replicate weights. ", reps$note, call. = FALSE)
  }
  out <- w * reps$factors
  attr(out, "type") <- reps$type
  attr(out, "scale") <- reps$scale
  attr(out, "rscales") <- reps$rscales
  out
}

#' The `.weight` column replicate weights are built on
#'
#' A replicate multiplies each row's weight, so every row needs a finite one.
#'
#' @noRd
check_replicate_base <- function(sample) {
  fix <- paste0(" Replicate weights multiply `.weight`; set it to ",
                "1 / `.prob` for the design weights.")
  w <- sample[[".weight"]]
  if (is.null(w)) {
    stop_not_drawn(paste0("has no `.weight` column, so there are no ",
                          "weights to replicate.", fix))
  }
  if (!is.numeric(w)) {
    stop_not_drawn(paste0("has ", a_class(w), " `.weight` column, where ",
                          "draw() puts a number.", fix))
  }
  bad <- !is.finite(w)
  if (any(bad)) {
    stop_not_drawn(paste0("has `.weight` values that are missing or not ",
                          "finite in ", sum(bad), " row(s): ",
                          format_rows(which(bad)), ".", fix))
  }
  w
}

#' Population totals from the frame behind a sample
#'
#' Totals the columns of the frame a sample was drawn from, as stored with
#' the sample, in the form [survey::calibrate()] takes as `population`. A
#' design handed over by [as_svydesign()] can then be calibrated, or
#' post-stratified, to the frame without building the totals by hand.
#'
#' The totals run over the rows the design could have drawn, those with an
#' inclusion probability above 0 when the sample was drawn: the rows inside
#' the region of [design_spatial()] or the window of [design_temporal()],
#' those [design_reservoir()] read before `max_items`, and not the rows
#' `na_rm` dropped. Those are the rows the sample's weights stand for. For
#' a sample weighted conditional on the clusters selected (see "Design
#' weights" in [draw()]), the rows of the clusters left out count too,
#' since another selection could have drawn them. A missing value of a
#' variable in `formula` in any of those rows is refused, since its total
#' is unknown.
#'
#' A value of a column in `formula` that no sampled row has cannot be
#' calibrated to. `survey::calibrate()` stops on such a factor level, and
#' leaves out the total of such a character value with a warning that it
#' reordered the totals. Calibration itself is left to `survey`, which gives
#' the calibrated estimates their standard errors.
#'
#' @param sample A data frame returned by [draw()] with `weights = TRUE`.
#' @param formula A one-sided formula of frame columns, such as
#'   `~ w + site`, as [survey::calibrate()] takes it. A name in it that is
#'   not a column of the frame, such as `cutoff` in `~ I(w > cutoff)`, is
#'   read from the environment of `formula`, as [stats::model.frame()]
#'   reads it.
#' @param table `FALSE`, the default, for the column totals of the model
#'   matrix of `formula`, which [survey::calibrate()] takes. `TRUE` for the
#'   count of rows in each cell of the variables in `formula`, as a table,
#'   which [survey::postStratify()] and [survey::rake()] take.
#'
#' @return With `table = FALSE`, a named numeric vector: the totals of the
#'   columns of `model.matrix(formula)` over the rows the design could
#'   reach, `(Intercept)` being the number of those rows. With
#'   `table = TRUE`, a table of counts of those rows from [stats::xtabs()].
#'
#' @examplesIf requireNamespace("survey", quietly = TRUE)
#' pop <- data.frame(
#'   id = 1:400,
#'   site = rep(c("a", "b", "c", "d"), times = c(200, 100, 60, 40)),
#'   size = 20 + (1:400 * 53) %% 300
#' )
#' pop$spend <- 3 * pop$size + (1:400 * 37) %% 90
#' s <- draw(pop, design_simple(n = 60), seed = 1, weights = TRUE)
#' des <- as_svydesign(s)
#'
#' # Calibrate to the frame's row count and total size
#' totals <- frame_totals(s, ~size)
#' totals
#' cal <- survey::calibrate(des, ~size, population = totals)
#' survey::svytotal(~spend, cal)
#' ht_total(s, "spend")
#'
#' # Post-stratify on site
#' counts <- frame_totals(s, ~site, table = TRUE)
#' post <- survey::postStratify(des, ~site, counts)
#' survey::svytotal(~spend, post)
#'
#' @seealso [as_svydesign()]
#' @export
frame_totals <- function(sample, formula, table = FALSE) {
  check_flag(table, "table")
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    stop("`formula` must be a one-sided formula of frame columns, such as ",
         "`~ w + site`, not ",
         if (inherits(formula, "formula")) {
           paste0("`", deparse1(formula), "`, which has a left-hand side")
         } else {
           a_class(formula)
         }, ".", call. = FALSE)
  }
  vars <- setdiff(all.vars(formula), ".")
  if (table && !length(vars)) {
    stop("`formula` (`", deparse1(formula), "`) names no column, and ",
         "`table = TRUE` counts rows by the columns it names. Name one, as ",
         "in `~ site`.", call. = FALSE)
  }
  state <- drawn_state(sample)
  pop <- state$pop
  unknown <- setdiff(vars, names(pop))
  gone <- intersect(unknown, state$compact$dropped)
  if (length(gone)) {
    stop("`formula` uses ", quote_names(gone), ", which compact_sample() ",
         "dropped from the frame stored with `sample`. Name ",
         if (length(gone) > 1L) "them" else "it", " in `keep` when ",
         "compacting, or rebuild the sample from the full frame with ",
         "as_drawn().", call. = FALSE)
  }
  # A name that is not a column, such as `cutoff` in `I(w > cutoff)`, is
  # read from where the formula was written, as model.frame() reads it. A
  # function of that name is a misspelt column, not a value.
  env <- environment(formula)
  defined <- vapply(unknown, function(v) {
    !is.null(env) && exists(v, envir = env) &&
      !is.function(get(v, envir = env))
  }, logical(1))
  unknown <- unknown[!defined]
  if (length(unknown)) {
    many <- length(unknown) > 1L
    stop("`formula` uses ", quote_names(unknown), ", which ",
         if (many) "are not columns" else "is not a column",
         " of the frame `sample` was drawn from, or ",
         if (many) "values" else "a value", " defined where `formula` was ",
         "written. Use columns of that frame; the totals come from the copy ",
         "stored with the sample.", call. = FALSE)
  }
  # A sample saved by drawn 0.2.0 records no probabilities, so they are
  # worked out again from the design.
  pi_all <- state$pi_all %||% exact_inclusion(state$design, pop)
  reach <- which(reachable(pi_all))
  rows <- as.data.frame(pop)[reach, , drop = FALSE]
  mf <- stats::model.frame(formula, rows, na.action = stats::na.pass)
  gaps <- vapply(mf, function(x) anyNA(x), logical(1))
  if (any(gaps)) {
    n_bad <- sum(!stats::complete.cases(mf))
    stop(quote_names(names(mf)[gaps]), if (sum(gaps) > 1L) " are" else " is",
         " missing in ", n_bad, " of the ", length(reach), " frame rows the ",
         "design could reach, so ", if (sum(gaps) > 1L) "their" else "its",
         " population total is unknown. Leave ",
         if (sum(gaps) > 1L) "them" else "it", " out of `formula`, or ",
         "supply the population totals yourself.", call. = FALSE)
  }
  if (table) return(stats::xtabs(formula, data = rows))
  colSums(stats::model.matrix(formula, mf))
}

#' Backticked names for a message: "`a`", "`a` and `b`", "`a`, `b` and `c`"
#' @noRd
quote_names <- function(x) {
  shown <- paste0("`", x, "`")
  if (length(shown) < 2L) return(shown)
  paste(paste(shown[-length(shown)], collapse = ", "), "and",
        shown[length(shown)])
}

#' Columns `as_svydesign()` writes the design into
#'
#' Each holds what the mapping needs: strata, cluster numbers, primary
#' sampling units, the two stages' population counts and the row ids of the
#' second stage. They are reserved in the sample handed over.
#'
#' @noRd
survey_columns <- function() {
  c(".stratum", ".cluster", ".psu", ".fpc", ".fpc2", ".unit")
}

#' Refuse a sample with a column `as_svydesign()` would overwrite
#'
#' Writing over a user column of one of these names would hand `survey` the
#' design's values in its place. A design column of such a name is in the
#' sample too, unless dropped, and is read from the frame stored with the
#' sample, so renaming it in the sample leaves the design unchanged.
#'
#' @noRd
check_survey_columns <- function(sample) {
  clash <- intersect(survey_columns(), names(sample))
  if (!length(clash)) return(invisible(NULL))
  shown <- quote_names(clash)
  it <- if (length(clash) == 1L) "it" else "them"
  stop_not_drawn(paste0(
    "already has ", shown, "; as_svydesign() would overwrite ", it,
    ". Rename ", it, " in the sample; the design keeps its strata and ",
    "clusters, which come from the copy of the frame stored with the ",
    "sample."))
}

#' Refuse unnamed arguments to `as_svydesign()`
#'
#' `utils::modifyList()` drops unnamed elements, so `as_svydesign(s, ~cl)`
#' would hand `survey` the mapped design and ignore `~cl`.
#'
#' @noRd
check_named_dots <- function(dots) {
  nm <- names(dots) %||% character(length(dots))
  bad <- which(is.na(nm) | !nzchar(nm))
  if (!length(bad)) return(invisible(NULL))
  where <- if (length(bad) > 1L) {
    paste("Arguments", paste(bad[-length(bad)], collapse = ", "), "and",
          bad[length(bad)], "of `...` have no name")
  } else {
    paste("Argument", bad, "of `...` has no name")
  }
  stop(where, ", and as_svydesign() passes `...` to survey::svydesign() by ",
       "name. Name ", if (length(bad) > 1L) "them" else "it", ", as in ",
       "`nest = TRUE` or `ids = ~cl`.", call. = FALSE)
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
#' strata, intervals and clusters taken whole, and the units a size-based
#' design puts at probability 1 -- are flagged in `whole`. Each stratum of
#' them becomes one primary sampling unit with an `fpc` of 1, which reads the
#' same as a count of units and as a probability. `survey` gives such a
#' stratum no variance, and `survey::degf()` no degree of freedom, as
#' `design_df()` does.
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
      by_group(as.integer(bucket), levels(bucket))
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
      # over a fixed n describes it. `survey` models it directly, with
      # `survey::poisson_sampling()` over the whole sample's probabilities,
      # which as_svydesign() builds.
      pps_parts(dat, prefix, "poisson")
    },
    cluster = ,
    multistage = {
      # The ids are integer codes. `survey` makes a factor of the ids, which
      # compares them as text and would merge two numbers that print alike.
      if (!is.null(design$size)) {
        return(pps_stage_parts(design, units, pop, dat$.prob, prefix))
      }
      cl <- count_clusters(design, pop)
      n_clusters <- cl$total
      cluster <- cl$code[rows]
      ids <- stats::as.formula("~.cluster")
      if (design_type(design) == "cluster" && !design$balanced) {
        return(list(stratum = NULL, fpc = rep(n_clusters, nrow(dat)),
                    ids = ids, cluster = cluster))
      }
      sizes <- tabulate(cl$code, n_clusters)
      fixed <- design_type(design) == "multistage" &&
        multistage_fixed(design, sizes)
      two_stage <- fixed ||
        weighted_on_selection(design, cluster, sizes, dat$.prob)
      if (!two_stage || one_row_clusters(design, cluster, sizes, fixed)) {
        # The two-stage estimator does not apply. One row of a cluster
        # leaves the variation within it unmeasured, or a take that depends
        # on the clusters selected is weighted without regard to them, as in
        # a certainty design's `rest`. ht_total() falls back to the
        # delete-a-cluster jackknife, which is the ultimate-cluster variance
        # with no finite population correction, and `survey` computes it
        # from the clusters alone.
        return(whole_cluster_parts(
          list(stratum = NULL, fpc = NULL, ids = ids, cluster = cluster),
          dat$.prob, prefix))
      }
      # Two stages, also for balanced clusters, which subsample. `survey`
      # takes each cluster's sampling fraction from the rows sampled in it,
      # which is the take given the clusters selected, so its estimate is the
      # one weighted conditional on them.
      whole_cluster_parts(
        list(stratum = NULL, fpc = rep(n_clusters, nrow(dat)),
             fpc2 = sizes[cluster],
             ids = stats::as.formula("~.cluster + .unit"), cluster = cluster),
        dat$.prob, prefix)
    },
    simple = {
      # Rows use_reserves() found unavailable are outside the population.
      list(stratum = NULL,
           fpc = rep(length(available_rows(design, nrow(pop))), nrow(dat)),
           ids = one)
    },
    reservoir = {
      # Rows past `max_items` are never read, so the population is the rows
      # the stream reached, not the whole frame.
      list(stratum = NULL, fpc = rep(reservoir_reach(design, nrow(pop)),
                                     nrow(dat)), ids = one)
    },
    systematic = {
      if (replicated(design)) {
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
        # A `rest` of clusters with no finite population correction gets an
        # infinite population, which `survey` reads as none.
        out$fpc[!above] <- rest$fpc %||% Inf
        out$whole[!above] <- rest$whole %||% FALSE
        out$ids <- rest$ids
        out$extra <- rest$extra
        if (!is.null(rest$cluster)) {
          # The clusters, or the replicates of several starts, are `rest`'s
          # primary sampling units, and the certainty rows one more beside
          # them.
          out$cluster <- rep(NA_integer_, nrow(dat))
          out$cluster[!above] <- rest$cluster
        }
        if (!is.null(rest$fpc2)) {
          # The certainty rows are one unit taken whole at both stages: all
          # its rows, whether `survey` reads the second stage as counts or
          # as probabilities.
          counts <- any(c(rest$fpc, rest$fpc2) > 1)
          out$fpc2 <- numeric(nrow(dat))
          out$fpc2[above] <- if (counts) sum(above) else 1
          out$fpc2[!above] <- rest$fpc2
        }
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
#' One row from a selected cluster of several leaves the variation within
#' clusters unmeasured, and `ht_total()` falls back to the delete-a-cluster
#' jackknife over the clusters below certainty, with no finite population
#' correction. That is `survey`'s one-stage design of the clusters with no
#' `fpc`, which such a sample gets, with any cluster taken whole in a stratum
#' of its own. A cluster taken with certainty but not whole then has a single
#' unit, and both packages decline; it keeps the two stages, where `survey`
#' stops on it. `prob` holds the sampled rows' inclusion probabilities.
#'
#' @noRd
pps_stage_parts <- function(design, units, pop, prob, prefix = "") {
  cu <- cluster_units(design, units, pop)
  st <- cluster_stage(design, pop)
  p1 <- st$pi1[cu$code]
  N <- st$sizes[cu$code]
  share <- st$m / N
  certain <- is_certain(p1)
  taken <- tabulate(cu$code, length(st$pi1))[cu$code]
  if (any(taken < 2L & taken < N) && !any(certain & taken < N)) {
    return(whole_cluster_parts(
      list(stratum = NULL, fpc = NULL, ids = stats::as.formula("~.cluster"),
           cluster = cu$psu),
      prob, prefix))
  }
  list(stratum = paste0(prefix, ifelse(certain, paste0("(certain) ", cu$code),
                                       "(pps)")),
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
  pps <- if (design$method == "poisson") "poisson" else "brewer"
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

#' Clusters taken whole, in a stratum of their own
#'
#' A cluster is taken whole when every one of its sampled rows is certain,
#' which needs every cluster selected and all of the cluster's rows drawn.
#' Such a cluster is in every sample and adds no variance, and `design_df()`
#' leaves it out of the degrees of freedom, as the jackknife holds it fixed.
#' `survey` would count it as one more primary sampling unit at the first
#' stage, and, in the one-stage form with no `fpc`, give it a share of the
#' variance between clusters.
#'
#' The clusters taken whole go in the stratum `(certain)` with an `fpc` of 1,
#' where `as_svydesign()` makes them one primary sampling unit, and the others
#' in the stratum `(clusters)`. In two stages that stratum's `fpc` is the
#' number of clusters left, all of them selected, and the shared unit is
#' taken whole at the second stage too: its `fpc2` is its number of rows. In
#' one stage, with no finite population correction, the other clusters get an
#' infinite population, which `survey` reads as none. `parts` is the mapping
#' without them, and `prob` the sampled rows' inclusion probabilities.
#'
#' @noRd
whole_cluster_parts <- function(parts, prob, prefix = "") {
  whole <- taken_whole(parts$cluster, prob)
  if (!any(whole) || all(whole)) return(parts)
  parts$stratum <- paste0(prefix, ifelse(whole, "(certain)", "(clusters)"))
  parts$whole <- whole
  if (is.null(parts$fpc)) {
    parts$fpc <- ifelse(whole, 1, Inf)
    return(parts)
  }
  parts$fpc <- ifelse(whole, 1, length(unique(parts$cluster[!whole])))
  parts$fpc2[whole] <- sum(whole)
  parts
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
#' Every other row keeps its unit in `unit`: its cluster's number, or the row
#' itself when `unit` is `NULL`. The shared units are numbered past the
#' others, so the ids are unique across strata and the units nest in their
#' strata without `nest = TRUE`.
#'
#' @noRd
shared_psu <- function(stratum, whole, unit = NULL) {
  psu <- unit %||% seq_along(whole)
  last <- max(0L, psu[!whole])
  psu[whole] <- last + match(stratum[whole], unique(stratum[whole]))
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
