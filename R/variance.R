# ---- Horvitz-Thompson estimation -------------------------------------------

#' Estimate a population total from a sample
#'
#' Forms the Horvitz-Thompson total `sum(y / pi)` with a standard error, a
#' confidence interval and a design effect, using a variance estimator matched
#' to how the design actually randomised.
#'
#' @section Variance:
#' Each design gets the estimator the survey-sampling literature recommends for
#' it, and `method` in the result says which one was used:
#'
#' \describe{
#'   \item{`"analytic"`}{Exact and design-unbiased. Fixed-size designs with
#'     known joint probabilities use the Sen-Yates-Grundy estimator. For
#'     simple random, reservoir and spatial samples, strata and intervals it
#'     reduces to the textbook formula built from the variance of `y` within
#'     each stratum or interval, and for an equal-take [design_multistage()]
#'     to the two-stage formula built from the variance between the estimated
#'     cluster totals and within clusters (Särndal, Swensson and Wretman
#'     1992). Those formulas are what is computed, in time and memory that
#'     grow with the number of sampled rows rather than with its square.
#'     Poisson sampling, whose rows are independent, uses
#'     `sum((1 - pi) / pi^2 * y^2)`; cluster designs apply the estimator at
#'     the cluster level, because taking whole clusters makes the row count
#'     random; certainty designs hand the problem to `rest`, since certainty
#'     rows are in every sample and add nothing. A [design_systematic()]
#'     sample with `n_starts` of 2 or more uses the variance between its
#'     replicates' totals, the random group estimator (Wolter 2007), on
#'     `n_starts - 1` degrees of freedom.}
#'   \item{`"deville"`}{Systematic probability-proportional-to-size
#'     ([design_weighted()] with `method = "systematic"`). Its joint
#'     probabilities have no closed form, so this uses Deville's (1999)
#'     approximation, which needs only the first-order probabilities and is the
#'     best-performing of the standard approximations for high-entropy designs
#'     (Matei and Tillé 2005). Randomised systematic selection is close to high
#'     entropy (Hartley and Rao 1962). Unlike a plain jackknife, each row's
#'     contribution is scaled by its own `1 - pi`, so rows near certainty stop
#'     adding variance they do not have. With `strata` it is added up over
#'     the strata. [design_multistage()] with `size` applies it to the
#'     estimated totals of the clusters, and adds the variance within them.}
#'   \item{`"successive difference"`}{[design_systematic()]. No
#'     design-unbiased estimator exists — the design has one random start, and
#'     most pairs of rows can never be drawn together — so this uses the
#'     successive-difference approximation, which compares each sampled row with
#'     its neighbour in the order the design walked (Wolter 2007). It removes a
#'     trend along that order rather than counting it as noise, which is why it
#'     beats treating the sample as simple random: on a frame sorted by a
#'     variable related to `y`, the simple random formula overstated the
#'     variance two- to thirtyfold in the package's simulations, where this one
#'     stayed near the truth. Two limits no estimator from a single systematic
#'     sample can escape: it is understated if the frame cycles with a period
#'     matching the step, and when a smooth trend so dominates `y` that the
#'     random start is almost the only source of variation, it runs low (0.3 to
#'     0.7 of the truth when the trend spanned 25 noise standard deviations).
#'     When `interval` does not divide the number of rows walked, the sample
#'     size is random, one of two values, and the total moves with it in a way
#'     successive differences cannot see. A term for that,
#'     `mean(y / pi)^2 * (r / k) * (1 - r / k)` with `k` the interval and `r`
#'     the remainder, is added, and the note says so. Without it a constant
#'     `y` had a standard error of 0, and a variable with a coefficient of
#'     variation of 0.02 an interval that never covered the truth. A design
#'     given `n` draws exactly `n` rows at a step of `N / n` and needs no
#'     such term. [design_stratified()] with
#'     `within = design_systematic()` adds up this estimator over its
#'     strata.}
#'   \item{`"local mean"`}{[design_spread()]. A spatially balanced sample is
#'     well spread precisely because neighbours are rarely taken together, so
#'     this compares each sampled row with the local mean of its nearest
#'     sampled neighbours (Stevens and Olsen 2003; Grafström and Schelin 2014).
#'     With `strata` it runs within each stratum, and the strata's variances
#'     are added.}
#'   \item{`"jackknife"`}{The delete-one-unit jackknife, by stratum (JKn) where
#'     the design has strata, and deleting whole clusters where it has them.
#'     Used when requested, and as the fallback under `variance = "auto"`.}
#' }
#'
#' Rows taken with certainty are in every possible sample and add nothing to
#' the variance. That covers rows above a [design_certainty()] threshold, rows
#' whose size puts them past certainty in [design_weighted()] or
#' [design_spread()], strata taken whole and a census. Every estimator leaves
#' them out. The jackknife holds them fixed instead of deleting them, and its
#' finite population correction counts only the other units. A sample made of
#' nothing else reports a variance of 0, with a note that the total is exact.
#'
#' The analytic estimator declines — returning `NA` with the reason in `note`
#' — rather than understate. Only a decline sends `variance = "auto"` to the
#' jackknife. Any other error inside an estimator, such as running out of
#' memory, stops `ht_total()` with that error. Two cases matter in practice:
#'
#' * **A stratum or interval with a single sampled row.** There is no
#'   variation inside it to measure, and simply leaving it out, as the
#'   Sen-Yates-Grundy sum does, reports a standard error that is too small —
#'   zero, if every stratum is like that. Draw at least two per stratum
#'   (`min_per_stratum = 2`, `per_interval = 2`).
#' * **One row per cluster in a multistage design.** Two rows of one cluster
#'   are never drawn together, so the within-cluster variation is
#'   unmeasurable. With `f1 = n_clusters / N_clusters`, the share of clusters
#'   selected, the exact estimator understates by `f1` times the
#'   within-cluster component of the variance. `variance = "auto"` falls back
#'   to the delete-a-cluster jackknife, the ultimate-cluster approximation
#'   standard in survey practice. It has no first-stage finite population
#'   correction, so it overstates by `f1` times the between-cluster component,
#'   the variance the cluster totals would have if clusters were drawn with
#'   replacement. The excess is small when few clusters are selected. With 30
#'   of 40 clusters whose totals differed widely, the jackknife averaged 3.3
#'   times the true variance. Whenever `f1` is above 0.2, the note gives the
#'   factor by which the jackknife inflates the between-cluster part. When
#'   every cluster is selected, the clusters are strata and the variance
#'   between them is zero by design, so [design_stratified()] with the cluster
#'   column as `strata` and two rows in each gives an exact variance instead.
#'
#' A Sen-Yates-Grundy estimate that comes out negative is treated the same way:
#' a failure of the estimator on an unlucky sample, not a variance.
#'
#' @section Confidence intervals:
#' Intervals use Student's t with the design's degrees of freedom — the number
#' of primary sampling units minus the number of strata — rather than the
#' normal distribution (Korn and Graubard 1999). The difference is negligible
#' with hundreds of units and decisive with a handful: at four clusters a
#' normal interval covers about 87% where it claims 95%, and the t interval
#' restores most of that. The degrees of freedom used are reported as `df`;
#' pass `df = Inf` for the normal interval.
#'
#' Units taken with certainty count in neither number. They are in every
#' sample and carry no variance, and counting them made the interval too
#' narrow. With eight certain units among ten, systematic PPS reported 9
#' degrees of freedom where the variance rests on two units, and its interval
#' covered 79%. A stratum or interval with no other unit drops out of the
#' count of strata. With no unit below certainty the variance is 0 and `df`
#' is `Inf`. [as_svydesign()] hands certain units to `survey` in a form that
#' makes `survey::degf()` give the same count, except for a census, where it
#' gives `n - 1`.
#'
#' An estimate whose variance is 0 is exact and gets a zero-width interval
#' whatever `df` is. Otherwise a design that leaves 0 degrees of freedom, such
#' as Poisson sampling that drew one row, gives a standard error but no
#' interval, and `note` says why. `df` set by hand must be at least 1. A
#' proportion of 0 or 1 is the exception to the first rule; see
#' "Proportions" in [ht_mean()].
#'
#' `alternative = "less"` gives a one-sided upper bound, the figure an auditor
#' signs off against, and `"greater"` a lower bound. The other end of `ci` is
#' then `-Inf` or `Inf`, or 0 or 1 for a proportion, and the printout labels
#' the bound. A one-sided bound at `level` is the end of a two-sided interval
#' at `1 - 2 * (1 - level)`, so the 95% upper bound is the upper end of the
#' 90% interval, which `level = 0.90` also gives.
#'
#' @section Rows the design cannot reach:
#' Rows with inclusion probability 0, such as rows past `max_items`, outside
#' a region or time window, dropped by `na_rm` or in a stratum allocated no
#' rows, can never be drawn, so a total covers the other rows only.
#' [draw()] warns about some of them by name. The estimate keeps their
#' number in `unreachable` and prints it on a line of its own, since a sample
#' saved and read back, or drawn with the warning muffled, would otherwise
#' give no sign of them. The design effect is measured against a simple
#' random sample of the rows the design can reach; see [deff()].
#'
#' @section An empty sample:
#' A design of random size, such as Poisson sampling or a systematic sample
#' whose random start falls past the last row, can draw nothing. `ht_total()`
#' then returns a total of 0, which keeps the estimator unbiased over repeated
#' draws, with `df` 0, the variance, standard error and interval `NA`, and a
#' note that the sample is empty. [ht_mean()] refuses the Hajek mean of an
#' empty sample, which is 0 / 0; `estimator = "ht"` gives 0. An empty sample
#' from a design that gives every row a probability of 0, such as one with
#' `n = 0`, is refused, and so is any empty sample in [as_svydesign()].
#'
#' @section Domains:
#' `by` estimates for each level of one or more columns — sites, regions,
#' months — without subsetting the sample first. Subsetting a sample and
#' treating the piece as its own sample gets the variance wrong, because the
#' number of rows that happen to fall in a domain is itself random. Here each
#' domain's total is estimated from `y * (row in domain)` over the whole sample,
#' so the same design and the same variance machinery apply (Särndal, Swensson
#' and Wretman 1992, section 10.3).
#'
#' Domains come out in key order: a factor's levels, numbers and dates in
#' sequence, and text sorted byte by byte as in the C locale, the same in every
#' session. A domain with no sampled row is left out, as `survey::svyby()`
#' leaves it out.
#'
#' A stratum or interval with a single sampled row stops a domain's variance,
#' as it stops the whole sample's, unless it holds no unit of the domain. That
#' is read from the `by` columns of the frame stored with the sample, so it
#' needs them in the frame, with the same values as in the sample. A domain
#' made of strata taken whole is then exact, with a zero-width interval. A
#' lonely stratum whose one sampled row lies outside the domain still stops
#' it when other units of the stratum belong to the domain, because their
#' variation went unmeasured.
#'
#' A domain's Hajek mean ([ht_mean()]) needs sampled rows in at least two
#' primary sampling units, rows or clusters. Its linearised residuals sum to
#' zero inside the unit that holds the whole domain, so with one unit every
#' estimator gives a standard error of about 0, and `survey` does the same.
#' Such a domain gets no standard error, with a note. Totals and the
#' `estimator = "ht"` mean are not affected. A domain with a handful of
#' sampled rows has a standard error that is itself unreliable, and the
#' interval covers less often than `level` says. In one simulation of a
#' skewed variable, a 95% interval for the mean of a domain with two or three
#' sampled rows covered 65% of the time, and 85% with four or five.
#'
#' A sample grouped with `dplyr::group_by()` is estimated as a whole, with a
#' warning of class `drawn_grouped_df`; the groups play no part. Use `by`.
#'
#' The strata, clusters and other pieces the variance reads from the frame
#' are worked out once for all the domains. For simple random, reservoir,
#' spatial, stratified and temporal samples every domain's variance comes
#' from one pass over the sample, in time that grows with the number of
#' sampled rows rather than with rows times domains.
#'
#' @section Contrasts:
#' `contrast` weighs the domain estimates of `by` and gives their weighted
#' sum, such as the difference between two sites' totals or means, with a
#' standard error, an interval and `df` as for any other estimate. Each
#' domain's estimate is the total, over the whole sample, of a variable that
#' is 0 outside the domain: `y` for a total, the linearised residual for a
#' Hajek mean and `y / N_d` for a Horvitz-Thompson mean. The weighted sum of
#' the estimates is the total of the same weighted sum of those variables, so
#' one run of the design's variance estimator gives its variance, with the
#' covariances between domains included. For simple random, stratified,
#' cluster and multistage samples the standard error equals that of
#' `survey::svycontrast()` applied to `survey::svyby(..., covmat = TRUE)`.
#'
#' The interval is symmetric, for a proportion too, since a difference of
#' proportions can be negative. Domains left out of `contrast` have a weight
#' of 0. A name that is not a domain with a sampled row is refused, and so is
#' a contrast whose weights are all 0. When a Hajek mean that the contrast
#' weighs rests on one primary sampling unit, the contrast gets no standard
#' error, for the reason given under "Domains".
#'
#' @section Which columns are read:
#' `ht_total()` and [ht_mean()] divide by `.prob`, the inclusion probabilities
#' [draw()] attached, and read strata, clusters, time stamps and the other
#' design columns from the copy of the frame stored with the sample. From the
#' sample itself they read only `y` and `by`, so renaming, recoding or dropping
#' a design column in the sample changes nothing.
#'
#' `.weight` is not used. A weight adjusted after the draw, for nonresponse or
#' by trimming, and written into `.weight` is ignored with a warning of class
#' `drawn_weight_edited`; [as_svydesign()] hands `.weight` to `survey`, which is
#' where adjusted weights are estimated from. A sample whose `.prob` was changed
#' or removed is refused. So is one whose design gives the frame different
#' inclusion probabilities in the current session than it did at draw time,
#' which a change of [sf::sf_use_s2()] does to a [design_spatial()] sample: the
#' variance is recomputed from the design and would describe another one.
#'
#' @param sample A data frame returned by [draw()] with `weights = TRUE`. See
#'   "Which columns are read".
#' @param y The variable to total: a column name, or a numeric vector as long as
#'   `sample`.
#' @param variance How to compute it. `"auto"` uses the design's own estimator
#'   (see "Variance") and falls back to the jackknife when that declines;
#'   `"analytic"` insists on the design's own estimator, returning `NA` with
#'   the reason in `note` rather than falling back; `"jackknife"` always
#'   resamples; `"none"` skips it. The result reports which was used in
#'   `method`, and `"none"` when nothing could produce a figure.
#' @param level Confidence level for the interval.
#' @param by Optional column name(s) in `sample` defining domains. When given,
#'   the result is one row per domain. See "Domains".
#' @param df Degrees of freedom for the interval. `NULL` uses the design's
#'   (see "Confidence intervals"); `Inf` gives a normal interval; otherwise a
#'   single number of 1 or more.
#' @param alternative `"two.sided"` for a confidence interval, `"less"` for a
#'   one-sided upper bound or `"greater"` for a lower bound, each at `level`.
#'   See "Confidence intervals".
#' @param contrast Optional weights over the `by` domains, as a named numeric
#'   vector such as `c(a = 1, b = -1)` for domain `a` less domain `b`. A name
#'   is the domain's value of a single `by` column, or the values of several
#'   joined by `" . "`, as in `"a . east"`, the way the notes of a `by` table
#'   name domains. Needs `by`. See "Contrasts".
#'
#' @return A list with a `print()` method, holding:
#'   \describe{
#'     \item{`total`}{The Horvitz-Thompson total, `sum(y / pi)`.}
#'     \item{`variance`, `se`, `ci`, `level`, `alternative`, `df`}{Its
#'       estimated variance, standard error, confidence interval or one-sided
#'       bound, and the degrees of freedom the interval used. `NA` where the
#'       design supports no variance.}
#'     \item{`ci_scale`}{The scale the interval was formed on: `"linear"`
#'       for a total, and `"logit"` or `"beta"` for a proportion from
#'       [ht_mean()].}
#'     \item{`n`, `design`}{Rows used, and the design's type.}
#'     \item{`deff`}{The design effect — see [deff()].}
#'     \item{`method`}{Which estimator produced the variance. See "Variance".}
#'     \item{`note`}{Why a variance is missing, what approximation was used,
#'       or which fallback was taken. `NULL` when an exact estimator applied
#'       cleanly.}
#'     \item{`unreachable`}{Rows of the frame with inclusion probability 0,
#'       which the estimate does not cover; `NA` for a sample saved by drawn
#'       0.2.0. See "Rows the design cannot reach".}
#'   }
#'   With `by`, a data frame of class `drawn_by` instead: one row per domain,
#'   with the domain columns, `n`, `total`, `se`, `ci_lower`, `ci_upper` and
#'   `method`, and the level, alternative, scale, notes and unreachable count
#'   as attributes.
#'
#'   With `contrast`, a list of class `drawn_contrast` with a `print()`
#'   method, holding `contrast`, the weighted sum; `weights` and `estimates`,
#'   the weight and estimate of each domain it weighs; `what` (`"total"` or
#'   `"mean"`), `estimator` for a mean, and `by`; and `variance`, `se`, `ci`,
#'   `level`, `alternative`, `df`, `n`, `ci_scale`, `design`, `method`, `note`
#'   and `unreachable` as above. `deff` is `NA`.
#'
#' @references
#' Horvitz, D. G. and Thompson, D. J. (1952). A generalization of sampling
#' without replacement from a finite universe. *Journal of the American
#' Statistical Association*, 47, 663–685.
#'
#' Sen, A. R. (1953). On the estimate of the variance in sampling with varying
#' probabilities. *Journal of the Indian Society of Agricultural Statistics*,
#' 5, 119–127.
#'
#' Yates, F. and Grundy, P. M. (1953). Selection without replacement from
#' within strata with probability proportional to size. *Journal of the Royal
#' Statistical Society B*, 15, 253–261.
#'
#' Hartley, H. O. and Rao, J. N. K. (1962). Sampling with unequal
#' probabilities and without replacement. *Annals of Mathematical Statistics*,
#' 33, 350–374.
#'
#' Deville, J.-C. (1999). Variance estimation for complex statistics and
#' estimators: linearization and residual techniques. *Survey Methodology*,
#' 25, 193–203.
#'
#' Matei, A. and Tillé, Y. (2005). Evaluation of variance approximations and
#' estimators in maximum entropy sampling with unequal probability and fixed
#' sample size. *Journal of Official Statistics*, 21, 543–570.
#'
#' Wolter, K. M. (2007). *Introduction to Variance Estimation*, 2nd ed.
#' Springer.
#'
#' Korn, E. L. and Graubard, B. I. (1999). *Analysis of Health Surveys*.
#' Wiley.
#'
#' Särndal, C.-E., Swensson, B. and Wretman, J. (1992). *Model Assisted Survey
#' Sampling*. Springer.
#'
#' @examples
#' set.seed(1)
#' pop <- data.frame(
#'   id = 1:200,
#'   site = rep(c("a", "b"), times = c(150, 50)),
#'   spend = round(stats::runif(200, 10, 500))
#' )
#'
#' s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
#' ht_total(s, "spend")
#'
#' sum(pop$spend)   # the truth
#'
#' # One estimate per site, from the same sample
#' ht_total(s, "spend", by = "site")
#' tapply(pop$spend, pop$site, sum)
#'
#' # How much more site a spent than site b
#' ht_total(s, "spend", by = "site", contrast = c(a = 1, b = -1))
#'
#' @seealso [ht_mean()], [deff()], [joint_prob()], [inclusion_prob()]
#' @export
ht_total <- function(sample, y, variance = c("auto", "analytic", "jackknife",
                                            "none"), level = 0.95,
                     by = NULL, df = NULL,
                     alternative = c("two.sided", "less", "greater"),
                     contrast = NULL) {
  variance <- check_choice(variance,
                           c("auto", "analytic", "jackknife", "none"),
                           "variance")
  alternative <- check_choice(alternative, c("two.sided", "less", "greater"),
                              "alternative")
  check_contrast_by(contrast, by)
  parts <- ht_prepare(sample, y, level, df)
  if (!is.null(by)) {
    return(estimate_by(parts, by, variance, level, df, what = "total",
                       alternative = alternative, contrast = contrast))
  }
  warn_grouped(sample, "ht_total")
  yv <- parts$y
  pi_i <- parts$pi
  total <- sum(yv / pi_i)
  estimate_one(parts, total, z = yv, variance, level, df, what = "total",
               class = "drawn_ht", deff_y = yv, alternative = alternative)
}

#' Variance, interval and design effect for one estimate
#'
#' `z` is the variable whose Horvitz-Thompson total has the estimate's
#' variance: `y` itself for a total, the linearised residual for a ratio.
#' `var_out`, when given, is a variance already settled, such as the decline
#' from `ratio_decline()`, and no estimator is run. `divide_by` is what a mean
#' divides its total by, for the design effect.
#'
#' @noRd
estimate_one <- function(parts, est, z, variance, level, df, what, class,
                         deff_y = NULL, extra = list(), proportion = FALSE,
                         ci_method = "logit", alternative = "two.sided",
                         var_out = NULL, divide_by = parts$reach) {
  var_out <- var_out %||% if (nrow(parts$sample)) {
    ht_variance_dispatch(parts$design, parts$units, parts$pop, parts$rows, z,
                         parts$pi, variance)
  } else {
    empty_sample()
  }
  df <- df %||% design_df(parts$design, parts$units, parts$pop, parts$rows,
                          parts$pi)
  out <- finish_estimate(est, var_out, parts, level, df, what = what,
                         class = class, y = deff_y, proportion = proportion,
                         ci_method = ci_method, alternative = alternative,
                         divide_by = divide_by)
  for (nm in names(extra)) out[[nm]] <- extra[[nm]]
  out
}

#' Route to the requested variance estimator, falling back where allowed
#'
#' `units` holds the frame's copy of each sampled row, in the sample's order.
#' The jackknife reads the clustering and stratification columns off it, and
#' without them it would delete one row at a time across strata and overstate
#' the variance many times over. They come from the frame rather than the
#' sample, so a column renamed or recoded in the sample after the draw cannot
#' change which rows share a stratum or cluster.
#'
#' An estimator declines by returning `NA` with a note, or by raising a
#' condition of class `drawn_declined` (see `stop_declined()`), as the joint
#' probabilities of a design without a closed form do. Only that class is
#' caught. Any other error reaches the caller, so running out of memory or a
#' bug cannot pass for a decline and send the estimate to the jackknife.
#'
#' `domain` says which rows of the frame belong to the domain being
#' estimated: `NULL` for the whole sample, a logical vector over the frame, or
#' `NA` when the frame cannot say. A stratum with one sampled row and no unit
#' of the domain adds nothing to the domain's variance, so the estimators skip
#' it rather than declining.
#'
#' `ctx` is a `variance_context()` shared by the calls for one sample, or
#' `NULL` for an estimator to work out what it needs from the frame itself.
#'
#' @noRd
ht_variance_dispatch <- function(design, units, pop, rows, z, pi_i, variance,
                                 domain = NULL, ctx = NULL) {
  if (variance == "none") {
    return(list(variance = NA_real_, method = "none",
                note = "Variance not requested."))
  }
  if (variance == "jackknife") {
    jk <- jackknife_variance(design, units, pop, z, pi_i, domain, ctx)
    if (is.na(jk$variance)) jk$method <- "none"
    return(jk)
  }

  # A negative Sen-Yates-Grundy estimate, a declining estimator and an
  # overflow are failures rather than variances. Each gets a reason instead
  # of a number that sqrt() and deff() would then have to guess about.
  a <- tryCatch(ht_variance(design, pop, rows, z, pi_i, domain, ctx),
                drawn_declined = function(e) {
                  list(variance = NA_real_, note = conditionMessage(e))
                })
  a$method <- a$method %||% "analytic"
  reason <- if (is.na(a$variance) && !is.nan(a$variance)) {
    a$note %||% "The design's own estimator returned no figure."
  } else if (!is.finite(a$variance)) {
    overflow_note()
  } else if (a$variance < 0) {
    paste0("The Sen-Yates-Grundy estimator returned a negative variance (",
           signif(a$variance, 3), "), which it can do on an unlucky sample. ",
           "There is no analytic figure to report.")
  }
  if (is.null(reason)) return(a)

  if (variance == "analytic") {
    return(list(variance = NA_real_, method = "none", note = reason))
  }

  jk <- jackknife_variance(design, units, pop, z, pi_i, domain, ctx)
  first <- sub("\n.*", "", reason)
  if (is.na(jk$variance)) {
    # The jackknife declined too. Say so, and keep its reason rather than
    # claiming a fallback that did not happen.
    return(list(variance = NA_real_, method = "none",
                note = paste0("No variance is available for this sample. ",
                              first,
                              if (!identical(jk$note, first)) {
                                paste0(" ", jk$note)
                              })))
  }
  jk$note <- paste0("The design's own estimator declined, so the jackknife ",
                    "was used instead. ", first, " ", jk$note)
  jk
}

# ---- variance context -------------------------------------------------------

#' A store for what the variance estimators read from the frame
#'
#' The stratum allocation, interval and cluster codes, group sizes, the
#' systematic walk, the certainty split, the spread neighbours, the
#' jackknife's strata and the joint probabilities depend on the design, the
#' frame and the sampled rows, and not on the variable being totalled.
#' `estimate_by()` runs an estimator once per domain, and each run worked all
#' of it out again from the whole frame: 200 domains of a stratified sample
#' from 100,000 rows took 3 seconds, nearly all of it rebuilding the strata.
#' One context is made per call and handed to every run, and each estimator
#' stores a piece there the first time it needs it, with `from_context()`.
#'
#' A context holds the pieces for one design, frame and set of sampled rows,
#' with the same inclusion probabilities. A certainty design keeps the context
#' of its `rest`, over the rows below the threshold, inside its own.
#'
#' @noRd
variance_context <- function() new.env(parent = emptyenv())

#' A piece of a variance context, worked out on first use
#'
#' `value` is evaluated only when `ctx` does not hold `name` yet. With no
#' context it is evaluated every time, so an estimator called without one
#' works out the piece itself.
#'
#' @noRd
from_context <- function(ctx, name, value) {
  if (is.null(ctx)) return(value)
  if (!exists(name, envir = ctx, inherits = FALSE)) {
    assign(name, value, envir = ctx)
  }
  get(name, envir = ctx, inherits = FALSE)
}

#' The context of a design nested in another, such as a certainty `rest`
#' @noRd
sub_context <- function(ctx, name) {
  if (is.null(ctx)) return(NULL)
  from_context(ctx, name, variance_context())
}

#' Decline the variance of a Hajek mean that rests on one sampling unit
#'
#' The linearised residual of a ratio, `(y - mean) / N_hat` over the domain,
#' sums to zero inside any primary sampling unit that holds the whole domain.
#' Every variance estimator then returns 0, or a rounding error away from it,
#' however much `y` varies. A domain of one sampled row, a domain inside one
#' sampled cluster and a Poisson sample of one row all reported a zero-width
#' interval. `survey` gives the same 0.
#'
#' Units are counted as the jackknife deletes them, clusters or rows, over the
#' sampled rows in `member`. The mean is exact, and a variance of 0 right,
#' only when every unit of the domain that the design can reach was taken with
#' certainty. That is read from the draw-time probabilities over `frame`, a
#' logical vector marking the domain's frame rows: `TRUE` for the whole
#' sample, `NA` when the frame cannot say. `label` names the domain in the
#' note, and `psu` takes `sample_units()` worked out once for many domains.
#' Returns the decline, or `NULL` when the variance can go ahead.
#'
#' @noRd
ratio_decline <- function(parts, member, frame = TRUE, label = NULL,
                          psu = sample_units(parts$design, parts$units,
                                             parts$pop, parts$rows)) {
  if (length(unique(psu[member])) >= 2L) return(NULL)
  pi_all <- parts$pi_all
  if (!is.null(pi_all) && !anyNA(frame)) {
    reach <- pi_all > 0 & frame
    if (all(is_certain(pi_all[reach]))) return(NULL)
  }
  unit <- unit_noun(parts$design)
  core <- paste0("A Hajek mean's linearised residuals sum to zero over it, so ",
                 "every variance estimator would give a standard error of ",
                 "about 0, whatever the spread of `y`.")
  note <- if (is.null(label)) {
    paste0("The sample has one ", unit, ". ", core, " A standard error needs ",
           "at least two.")
  } else {
    paste0("The domain `", label, "` has ",
           if (unit == "row") "one sampled row" else
             paste0("sampled rows in one ", unit, " only"), ". ", core,
           if (unit == "row") {
             " A standard error needs at least two sampled rows in the domain."
           } else {
             paste0(" A standard error needs the domain's sampled rows to ",
                    "span at least two ", unit, "s.")
           })
  }
  list(variance = NA_real_, method = "none", note = note)
}

#' Each sampled row's primary sampling unit
#'
#' Clusters for a cluster or multistage design, as the jackknife deletes
#' them, the replicates of a systematic design with several starts, and
#' otherwise rows. Under a certainty design each row at or above the
#' threshold is a unit of its own and the rows below it are grouped as
#' `rest` groups them. `pop` is the frame and `rows` the sampled frame rows,
#' which place a systematic sample's rows on its walk.
#'
#' @noRd
sample_units <- function(design, units, pop, rows) {
  if (replicated(design)) {
    return(systematic_replicates(design, systematic_positions(design, pop),
                                 rows))
  }
  if (!inherits(design, "drawn_design_certainty")) {
    return(jackknife_groups(design, units, pop))
  }
  v <- units[[design$above]]
  below <- !is.na(v) & v < design$threshold
  out <- -seq_len(nrow(units))
  if (any(below)) {
    sp <- certainty_split(design, pop)
    out[below] <- sample_units(design$rest, units[below, , drop = FALSE],
                               sp$data[sp$rest, , drop = FALSE],
                               match(rows[below], sp$keep[sp$rest]))
  }
  out
}

#' What a design's primary sampling unit is called in a message
#' @noRd
unit_noun <- function(design) {
  clustered <- c("drawn_design_cluster", "drawn_design_multistage")
  if (inherits(design, clustered)) return("cluster")
  if (replicated(design)) return("replicate")
  if (inherits(design, "drawn_design_certainty") &&
      (inherits(design$rest, clustered) || replicated(design$rest))) {
    return("primary sampling unit")
  }
  "row"
}

#' The reason given for a variance too large to represent
#'
#' Squaring `y / pi` passes the largest double, about `1.8e308`, once the
#' values pass about `1e154`, and the sum comes out infinite or `NaN`. The
#' analytic estimators and the jackknife give the same note, so a fallback
#' from one to the other states it once.
#'
#' @noRd
overflow_note <- function() {
  paste0("The variance overflowed, because squaring `y / pi`, or totals ",
         "built from it, passes the largest number R can hold. Rescale `y`, ",
         "to thousands for example, and estimate again.")
}

#' Degrees of freedom for a design's confidence interval
#'
#' Primary sampling units minus strata, the usual rule (Korn and Graubard
#' 1999), counting only units that could have been left out. A unit taken with
#' certainty is in every sample, so it carries no variance and no degree of
#' freedom: the units `design_weighted()` and `design_spread(size =)` take with
#' certainty, strata and intervals taken whole, and the rows above a
#' certainty design's threshold. A stratum holding none of the other units
#' drops out of the count of strata. Counting them put eight certain units of
#' ten at 9 degrees of freedom where the same selection built with
#' `design_certainty()` had 1, and the t interval covered 79%.
#'
#' A cluster counts when any of its sampled rows is below certainty. A
#' two-stage design that takes every cluster, but not every row in them, has
#' the number of clusters less one. With `size`, a cluster taken with
#' certainty is a stratum whose units are its sampled rows, as
#' `cluster_units()` lays it out.
#'
#' With no unit below certainty the variance is exactly 0 and the result is
#' `Inf`. An empty sample, and a certainty design whose `rest` drew none of the
#' rows below the threshold, have no variance to build an interval on and get
#' 0. A systematic design with `n_starts` of 2 or more has the number of
#' starts less one, the degrees of freedom of its random group variance.
#'
#' `as_svydesign()` hands each group of certain units to `survey` as one
#' primary sampling unit in a stratum of its own, so `survey::degf()` gives the
#' same count. Strata, those of a stratified or spread design or a temporal
#' design's intervals, and clusters are read from `units`, the frame's copy
#' of the sampled rows; `pi_i` holds their inclusion probabilities.
#'
#' @noRd
design_df <- function(design, units, pop, rows, pi_i) {
  if (inherits(design, "drawn_design_certainty")) {
    sp <- certainty_split(design, pop)
    free <- !(rows %in% sp$keep[sp$take])
    if (!any(free)) return(if (length(sp$rest)) 0 else Inf)
    return(design_df(design$rest, units[free, , drop = FALSE],
                     sp$data[sp$rest, , drop = FALSE],
                     match(rows[free], sp$keep[sp$rest]), pi_i[free]))
  }
  n <- nrow(units)
  if (n == 0L) return(0)
  if (replicated(design)) {
    # The random group variance has one degree of freedom fewer than the
    # number of starts, counting a replicate that drew no rows.
    return(if (any(!is_certain(pi_i))) systematic_starts(design) - 1 else Inf)
  }
  if (inherits(design, "drawn_design_multistage") && !is.null(design$size)) {
    cu <- cluster_units(design, units, pop)
    random <- cu$psu %in% cu$psu[!is_certain(pi_i)]
    if (!any(random)) return(Inf)
    return(max(length(unique(cu$psu[random])) -
                 length(unique(cu$stratum[random])), 0))
  }
  psu <- switch(design_type(design),
    cluster = ,
    multistage = cluster_codes(units[[design$clusters]]),
    seq_len(n)
  )
  stratum <- switch(design_type(design),
    stratified = ,
    spread = ,
    weighted = if (!is.null(design$strata)) {
      as.integer(group_key(units, design$strata))
    },
    temporal = temporal_bucket(design, units),
    NULL
  ) %||% rep("all", n)
  random <- psu %in% psu[!is_certain(pi_i)]
  if (!any(random)) return(Inf)
  max(length(unique(psu[random])) - length(unique(stratum[random])), 0)
}

#' The variance of an empty sample
#'
#' A design of random size, such as Poisson sampling or a systematic one whose
#' start falls past the last row, can draw nothing. The Horvitz-Thompson
#' estimate is then 0, which keeps it unbiased over repeated draws, but one
#' empty draw has no spread to measure.
#'
#' @noRd
empty_sample <- function() {
  list(variance = NA_real_, method = "none", note = paste0(
    "The sample is empty. The design could have drawn rows, but this draw ",
    "took none, so the estimate is 0 and has no variance. A design with a ",
    "fixed sample size, or a larger expected one, avoids this."))
}

#' @noRd
ht_variance <- function(design, data, rows, y, pi_i, domain = NULL,
                        ctx = NULL) UseMethod("ht_variance")

#' Sen-Yates-Grundy, for fixed-size designs with closed-form joint
#' probabilities
#'
#' The designs in the package that select by simple random sampling within
#' groups have methods of their own built on `srs_variance()` and
#' `two_stage_variance()`, which give the same figure in linear time and
#' memory. What reaches this method is a design with a joint matrix and no such
#' form, or one with no joint probabilities at all, whose `joint_inclusion()`
#' declines.
#'
#' The joint probabilities and the sums over pairs are n x n matrices of
#' doubles, about eight of them live at once with `syg_variance()`. Past
#' `pair_limit()` bytes projected, about 1 GB or 4,000 sampled rows, the
#' estimator declines with the size it would need, before anything is
#' allocated. Linux can stop R for want of memory before an allocation error
#' is raised, so waiting for one is no guard. The joint matrix is kept in
#' `ctx`, so the domains of one sample share it.
#'
#' @noRd
ht_variance.default <- function(design, data, rows, y, pi_i, domain = NULL,
                                ctx = NULL) {
  n <- length(rows)
  bytes <- 8 * 8 * as.numeric(n)^2
  if (bytes > pair_limit()) {
    stop_declined(paste0(
      "The Sen-Yates-Grundy estimator would need about ", format_gb(bytes),
      " of memory for its ", n, " x ", n,
      " matrices.\nUse variance = \"auto\" to fall back to the jackknife."))
  }
  pij <- from_context(ctx, "joint", joint_inclusion(design, data, rows))
  if (length(rows) < 2L) {
    return(list(variance = NA_real_,
                note = "A variance needs at least two sampled rows."))
  }
  syg_variance(y, pi_i, pij)
}

#' The variance of a sample whose rows were all taken with certainty
#'
#' Every such row is in every possible sample, so the total is known exactly.
#' `method` names the estimator that found this, or is left unset for the
#' design's own exact one.
#'
#' @noRd
all_certain <- function(method = NULL) {
  out <- list(variance = 0, note = paste0(
    "Every sampled row was taken with certainty, so the total is exact ",
    "rather than estimated."))
  out$method <- method
  out
}

#' The Sen-Yates-Grundy sum over pairs of sampled rows
#'
#' Kept as written, over the full joint matrix, as the reference the closed
#' forms below are tested against.
#'
#' @noRd
syg_variance <- function(y, pi_i, pij) {
  yk <- y / pi_i
  d <- outer(yk, yk, "-")^2
  num <- outer(pi_i, pi_i) - pij
  w <- num / pij
  w[!is.finite(w)] <- 0                 # a zero joint probability contributes nothing
  diag(w) <- 0
  list(variance = 0.5 * sum(w * d), note = NULL)
}

#' Variance of a total under simple random sampling within groups
#'
#' For a design that takes a simple random sample of `n_g` of the `N_g` rows
#' in each group, independently across groups, the Sen-Yates-Grundy sum comes
#' to `sum(N_g^2 * (1 - n_g / N_g) * s_g^2 / n_g)`, with `s_g^2` the variance
#' of `z` over the sampled rows of group `g`. This computes it in that form,
#' in time and memory linear in the rows, where the sum over pairs held
#' several n x n matrices: about 1.4 GB for a stratified sample of 4,000 rows.
#'
#' `group` is each sampled row's position in `n_g` and `N_g`. A group taken
#' whole adds exactly 0, a group of one row included. A group with one sampled
#' row from several has no pair of rows to compare and adds 0, as it does in
#' the sum over pairs; callers refuse it first with `lonely_note()` where it
#' matters. `z` is whatever is being totalled: `y`, `y` times a domain's
#' indicator, or a Hajek mean's linearised residual. The probability is the
#' same for every row of a group, so the formula holds for all three. An
#' integer `z` is summed as double, where `rowsum()` would overflow to `NA`.
#'
#' @noRd
srs_variance <- function(z, group, n_g, N_g) {
  z <- as.numeric(z)
  live <- sort(unique(group))
  g <- match(group, live)
  m <- tabulate(g, length(live))
  n <- as.numeric(n_g)[live]
  N <- as.numeric(N_g)[live]
  use <- m >= 2L & n < N
  if (!any(use)) return(0)
  mu <- rowsum(z, g)[, 1L] / m
  ss <- rowsum((z - mu[g])^2, g)[, 1L]
  sum((N^2 * (1 - n / N) * ss / ((m - 1) * n))[use])
}

#' Variance of a total from a two-stage sample with an equal take
#'
#' `a` of `A` clusters by simple random sampling, then `m` of the `sizes[i]`
#' rows of each selected cluster `i`. The unbiased estimator is
#' `A^2 * (1 - a / A) * s_t^2 / a + (A / a) * sum(N_i^2 * (1 - m / N_i) *
#' s_i^2 / m)`, with `s_t^2` the variance of the estimated cluster totals
#' `N_i * mean(z_i)` and `s_i^2` the variance of `z` within cluster `i`. It is
#' the Sen-Yates-Grundy sum for this design, computed without the n x n joint
#' matrix. Selecting every cluster leaves nothing between them to vary, and a
#' cluster taken whole adds nothing within it.
#'
#' `cluster` is each sampled row's cluster code, indexing `sizes`.
#'
#' @noRd
two_stage_variance <- function(z, cluster, a, A, m, sizes) {
  z <- as.numeric(z)
  live <- sort(unique(cluster))
  g <- match(cluster, live)
  N_i <- as.numeric(sizes)[live]
  totals <- N_i * rowsum(z, g)[, 1L] / tabulate(g, length(live))
  between <- if (a < A) A^2 * (1 - a / A) * stats::var(totals) / a else 0
  within <- srs_variance(z, g, rep(m, length(live)), N_i)
  between + A / a * within
}

#' The variance of a simple random sample from one group of `N` rows
#'
#' Shared by the simple, reservoir and spatial designs, which differ only in
#' the `N` they draw `n` from. One sampled row has nothing to compare with
#' unless it is the whole group.
#'
#' @noRd
srs_one_group <- function(z, n, N) {
  if (length(z) < 2L && n < N) {
    return(list(variance = NA_real_,
                note = "A variance needs at least two sampled rows."))
  }
  list(variance = srs_variance(z, rep(1L, length(z)), n, N), note = NULL)
}

#' The groups a design takes a simple random sample within
#'
#' For the simple, reservoir and spatial designs, one group: `n_g` rows drawn
#' from the `N_g` the design reaches. For stratified and temporal designs,
#' `code` gives each frame row's stratum or interval as a position in `n_g`,
#' `N_g` and `labels`. `NULL` for any other design. The analytic estimators of
#' these designs and `srs_layout()`, which estimates every domain in one pass,
#' read their groups from here, so the two cannot disagree. A design drawn
#' with replacement gets its groups too; the callers decide what to do with
#' it.
#'
#' @noRd
srs_groups <- function(design, data, ctx = NULL) {
  from_context(ctx, "srs_groups", switch(design_type(design),
    simple = list(n_g = design$n, N_g = nrow(data)),
    reservoir = {
      # Rows past `max_items` are never read: the sample is a simple random
      # one from the rows before it.
      reach <- reservoir_reach(design, nrow(data))
      list(n_g = min(design$n, reach), N_g = reach)
    },
    spatial = {
      require_suggested("sf", "spatial sampling")
      N <- sum(spatial_inside(design, data))
      list(n_g = min(design$n, N), N_g = N)
    },
    stratified = {
      al <- stratum_allocation(design, data)
      list(code = al$code, labels = al$labels, n_g = al$n_alloc,
           N_g = al$sizes)
    },
    temporal = {
      bucket <- key_from_columns(list(temporal_bucket(design, data)))
      code <- as.integer(bucket)
      sizes <- tabulate(code, nlevels(bucket))
      list(code = code, labels = levels(bucket),
           n_g = pmin(sizes, design$per_interval), N_g = sizes)
    },
    NULL
  ))
}

#' @noRd
ht_variance.drawn_design_simple <- function(design, data, rows, y, pi_i,
                                            domain = NULL, ctx = NULL) {
  if (design$replace) return(NextMethod())
  g <- srs_groups(design, data, ctx)
  srs_one_group(y, g$n_g, g$N_g)
}

#' @noRd
ht_variance.drawn_design_reservoir <- function(design, data, rows, y, pi_i,
                                               domain = NULL, ctx = NULL) {
  g <- srs_groups(design, data, ctx)
  srs_one_group(y, g$n_g, g$N_g)
}

#' @noRd
ht_variance.drawn_design_spatial <- function(design, data, rows, y, pi_i,
                                             domain = NULL, ctx = NULL) {
  g <- srs_groups(design, data, ctx)
  srs_one_group(y, g$n_g, g$N_g)
}

#' Groups with one sampled row out of several
#'
#' `code` is each frame row's group as a position among `k` groups, and
#' `rows` the sampled frame rows. A group taken whole is not lonely. The one
#' rule for `lonely_note()` and `srs_layout()`.
#'
#' @noRd
lonely_groups <- function(code, rows, k) {
  tabulate(code[rows], k) == 1L & tabulate(code, k) > 1L
}

#' Refuse a stratified variance with a stratum that holds one sampled row
#'
#' Sen-Yates-Grundy sums over pairs, and a stratum with one sampled row has no
#' within-stratum pair: its contribution is silently zero. With one row in
#' every stratum the standard error came out as exactly 0 against a true
#' figure in the hundreds. `survey` fails in the same situation unless told
#' otherwise (`options(survey.lonely.psu)`); declining here keeps the two in
#' step and says what to change.
#'
#' Strata taken whole (`n_h == N_h`) are not lonely: they contribute no
#' variance because they have none.
#'
#' For a domain estimate, `domain` marks the frame rows in the domain. A
#' lonely stratum with none of them adds nothing to the domain's variance,
#' because the domain's `y` is 0 on every unit of it, so it is skipped. That
#' is decided from the frame, never from the sampled rows: a stratum whose one
#' sampled row lies outside the domain can still hold domain units that were
#' not drawn. `domain = NA` means the frame cannot say, and every lonely
#' stratum is refused.
#'
#' @noRd
lonely_note <- function(group, rows, noun, fix, domain = NULL) {
  # `noun` is c(singular, plural). `group` is each frame row's group, as a
  # factor or as values, and is counted by position.
  if (!is.factor(group)) {
    group <- key_from_columns(list(group))
  }
  code <- as.integer(group)
  k <- nlevels(group)
  lonely <- lonely_groups(code, rows, k)
  if (known_domain(domain, length(code))) {
    lonely <- lonely & tabulate(code[domain], k) > 0L
  }
  lonely_text(levels(group)[lonely], noun, fix, domain)
}

#' The refusal of lonely groups, naming them
#'
#' `lonely` holds the groups' key levels, in key order. Shared by the analytic
#' estimators and the jackknife, so the two refuse a sample in the same words
#' and `variance = "auto"` says it once.
#'
#' @noRd
lonely_text <- function(lonely, noun, fix, domain = NULL) {
  if (!length(lonely)) return(NULL)
  shown <- group_label(utils::head(lonely, 5L))
  many <- length(lonely) > 1L
  them <- if (many) "them" else "it"
  paste0(if (many) paste(length(lonely), noun[2]) else paste("The", noun[1]),
         " ", paste0("`", shown, "`", collapse = ", "),
         if (length(lonely) > 5L) ", ..." else "",
         " ", if (many) "each have" else "has",
         " a single sampled row, so the variation within ", them,
         if (identical(domain, NA)) {
           paste0(" cannot be measured. If ", if (many) "they hold" else
                    "it holds", " units of this domain, which cannot be told ",
                  "from the frame, leaving ", them, " out would understate ",
                  "the standard error. ")
         } else {
           paste0(" cannot be measured and leaving ", them, " out would ",
                  "understate the standard error. ")
         }, fix)
}

#' What to change so that no stratum or interval holds one sampled row
#' @noRd
lonely_fix <- function(design) {
  if (inherits(design, "drawn_design_multistage")) {
    return(paste0("Draw at least two rows per cluster (`n` of at least ",
                  "twice `n_clusters`)."))
  }
  if (inherits(design, "drawn_design_temporal")) {
    return("Draw at least two rows per interval (per_interval = 2).")
  }
  if (inherits(design, c("drawn_design_spread", "drawn_design_weighted"))) {
    return(paste0("Draw at least two rows below certainty per stratum (",
                  if (!has_counts(design)) "a larger `n`, or ",
                  "a count of 2 or more for each in `allocation`)."))
  }
  if (has_counts(design)) {
    return(paste0("Draw at least two rows per stratum (a count of 2 or more ",
                  "for each in `allocation`)."))
  }
  "Draw at least two rows per stratum (min_per_stratum = 2)."
}

#' Does `domain` mark which of `n` frame rows belong to a domain?
#' @noRd
known_domain <- function(domain, n) {
  is.logical(domain) && length(domain) == n && !anyNA(domain)
}

#' Variance within the groups of `srs_groups()`, refusing a lonely group
#'
#' The stratified and temporal estimators: each group is a simple random
#' sample of its `n_g` rows from its `N_g`, so the variance is
#' `srs_variance()` over groups. `noun` and `fix` word the refusal of a
#' lonely group. A design drawn with replacement is refused for a lonely
#' group too, and otherwise gets `NULL`, for its method to pass on.
#'
#' @noRd
srs_within_groups <- function(design, g, rows, y, domain, noun, fix) {
  group <- structure(g$code, levels = g$labels, class = "factor")
  note <- lonely_note(group, rows, noun, fix, domain)
  if (!is.null(note)) return(list(variance = NA_real_, note = note))
  if (isTRUE(design$replace)) return(NULL)
  list(variance = srs_variance(y, g$code[rows], g$n_g, g$N_g), note = NULL)
}

#' Stratified variance, from the allocation the design used
#'
#' Each stratum is a simple random sample of its `n_alloc` rows from its
#' `sizes` rows, both from `stratum_allocation()`, so the variance is
#' `srs_variance()` over strata. A lonely stratum is refused first.
#'
#' @noRd
ht_variance.drawn_design_stratified <- function(design, data, rows, y, pi_i,
                                                domain = NULL, ctx = NULL) {
  if (walks_strata(design)) {
    return(walk_strata_variance(design, data, rows, y, pi_i, domain, ctx))
  }
  out <- srs_within_groups(
    design, srs_groups(design, data, ctx), rows, y, domain,
    c("stratum", "strata"), lonely_fix(design))
  if (is.null(out)) return(NextMethod())
  out
}

#' Successive-difference variance of a sample walked within strata
#'
#' Each stratum is a systematic sample of its `n_h` rows from its `N_h` at a
#' step of `N_h / n_h`, from a start of its own, so the variance is the sum
#' over strata of `(1 - n_h / N_h) * n_h / (2 * (n_h - 1)) * sum(diff(u)^2)`,
#' `u` being `y / pi` over the stratum's sampled rows in walk order. That is
#' the systematic estimator with `n` applied to each stratum. A lonely
#' stratum is refused first, as for simple random sampling within strata,
#' and a stratum taken whole adds 0.
#'
#' @noRd
walk_strata_variance <- function(design, data, rows, y, pi_i, domain = NULL,
                                 ctx = NULL) {
  g <- srs_groups(design, data, ctx)
  group <- structure(g$code, levels = g$labels, class = "factor")
  note <- lonely_note(group, rows, c("stratum", "strata"), lonely_fix(design),
                      domain)
  if (!is.null(note)) return(list(variance = NA_real_, note = note))
  k <- length(g$N_g)
  pos <- from_context(ctx, "walk", walk_positions(design, data, g$code, k))
  h <- g$code[rows]
  o <- order(h, pos[rows])
  h <- h[o]
  u <- as.numeric(y / pi_i)[o]
  # Successive differences within a stratum, never across two.
  same <- h[-1L] == h[-length(h)]
  ss <- vapply(split(diff(u)[same]^2,
                     factor(h[-1L][same], levels = seq_len(k))),
               sum, numeric(1), USE.NAMES = FALSE)
  m <- tabulate(h, k)
  f <- as.numeric(g$n_g) / as.numeric(g$N_g)
  use <- m >= 2L
  list(variance = sum(((1 - f) * m / (2 * (m - 1)) * ss)[use]),
       method = "successive difference",
       note = paste0("Systematic sampling within strata has no ",
                     "design-unbiased variance estimator; this is the sum ",
                     "over strata of the successive-difference ",
                     "approximation. It would be understated if a stratum ",
                     "cycled with a period matching its step, or if a ",
                     "smooth trend along the sort order swamped all other ",
                     "variation."))
}

#' Temporal variance: intervals are strata
#'
#' Each interval is a simple random sample of `per_interval` of its rows, or
#' all of them when it holds fewer. A lonely interval is refused first.
#'
#' @noRd
ht_variance.drawn_design_temporal <- function(design, data, rows, y, pi_i,
                                              domain = NULL, ctx = NULL) {
  srs_within_groups(
    design, srs_groups(design, data, ctx), rows, y, domain,
    c("interval", "intervals"), lonely_fix(design))
}

#' Each sampled row's group, for every domain's variance in one pass
#'
#' Returns `NULL` when the design's analytic variance is not `srs_variance()`
#' over the groups of `srs_groups()`, as for a stratified design that walks
#' its strata, or when it declines for the whole
#' sample, as with one sampled row from a group of several. Otherwise
#' `group` gives each sampled row's position in `n_g` and `N_g`, `NA` for a
#' row taken with certainty above a certainty design's threshold, and
#' `lonely_rows` the frame rows in a group with one sampled row, whose domains
#' `srs_variance()` cannot measure.
#'
#' @noRd
srs_layout <- function(design, data, rows, ctx = NULL) {
  if (inherits(design, "drawn_design_certainty")) {
    sp <- certainty_parts(design, data, ctx)
    free <- !(rows %in% sp$keep[sp$take])
    pos <- match(rows[free], sp$keep[sp$rest])
    if (!any(free) || anyNA(pos)) return(NULL)
    out <- srs_layout(design$rest, sp$rest_data, pos, sub_context(ctx, "rest"))
    if (is.null(out)) return(NULL)
    group <- rep(NA_integer_, length(rows))
    group[free] <- out$group
    out$group <- group
    out$lonely_rows <- sp$keep[sp$rest][out$lonely_rows]
    return(out)
  }
  if (isTRUE(design$replace) || walks_strata(design)) return(NULL)
  g <- srs_groups(design, data, ctx)
  if (is.null(g)) return(NULL)
  if (is.null(g$code)) {
    if (length(rows) < 2L && g$n_g < g$N_g) return(NULL)
    return(list(group = rep(1L, length(rows)), n_g = g$n_g, N_g = g$N_g,
                lonely_rows = integer(0)))
  }
  lonely <- lonely_groups(g$code, rows, length(g$N_g))
  list(group = g$code[rows], n_g = g$n_g, N_g = g$N_g,
       lonely_rows = which(lonely[g$code]))
}

#' `srs_variance()` for many domains at once
#'
#' `z` holds each sampled row's value in its own domain: `y`, or the Hajek
#' residual or `y / N_d` of that domain. Domain `d`'s variable is `z` on its
#' rows and 0 on every other row, and `domain` gives each row's `d` among
#' `n_dom`. Within group `g` that variable has mean `S_gd / m_g`, with `S_gd`
#' the sum of `z` over the group's rows in `d` and `m_g` the group's sampled
#' rows, and its squared deviations are those of the rows in `d` plus
#' `(m_g - c_gd)` times the squared mean for the `c_gd` rows outside it. So
#' one `rowsum()` over the group and domain cells gives every domain's
#' variance in time linear in the rows, where a pass per domain was linear in
#' rows times domains. The figures equal `srs_variance()` on each domain's
#' variable, up to the order of summation.
#'
#' @noRd
srs_variance_by <- function(z, group, domain, n_dom, n_g, N_g) {
  z <- as.numeric(z)
  out <- numeric(n_dom)
  live <- sort(unique(group))
  g <- match(group, live)
  m <- tabulate(g, length(live))
  n <- as.numeric(n_g)[live]
  N <- as.numeric(N_g)[live]
  use <- m >= 2L & n < N
  rows <- use[g]
  if (!any(rows)) return(out)
  g <- g[rows]
  d <- domain[rows]
  z <- z[rows]
  # One cell per group and domain that share a sampled row. The arithmetic is
  # in doubles, which hold the cell numbers exactly.
  cell_id <- (g - 1) * as.numeric(n_dom) + d
  cells <- unique(cell_id)
  cell <- match(cell_id, cells)
  first <- match(seq_along(cells), cell)
  cg <- g[first]
  mu <- rowsum(z, cell, reorder = FALSE)[, 1L] / m[cg]
  ss <- rowsum((z - mu[cell])^2, cell, reorder = FALSE)[, 1L] +
    (m[cg] - tabulate(cell, length(cells))) * mu^2
  v <- N[cg]^2 * (1 - n[cg] / N[cg]) * ss / ((m[cg] - 1) * n[cg])
  cd <- d[first]
  out[unique(cd)] <- rowsum(v, cd, reorder = FALSE)[, 1L]
  out
}

#' Variance of a single-stage cluster total
#'
#' Whole clusters are taken, so the number of *rows* is random whenever the
#' clusters differ in size. Sen-Yates-Grundy assumes a fixed size, and applied
#' row by row here it understates the variance badly -- by a factor of five on a
#' frame whose clusters vary from 2 to 10 rows -- and can return a negative
#' number or a zero-width interval.
#'
#' The cluster is the sampling unit, so the estimator belongs at that level: the
#' clusters are a simple random sample of `a` from `A`, and the quantity summed
#' over them is each cluster's contribution to the total. That is the textbook
#' form, and it is algebraically identical to the delete-a-cluster jackknife.
#'
#' @noRd
ht_variance.drawn_design_cluster <- function(design, data, rows, y, pi_i,
                                             domain = NULL, ctx = NULL) {
  cl <- from_context(ctx, "clusters", count_clusters(design, data))
  A <- cl$total
  code <- cl$code[rows]
  if (anyNA(code)) {
    return(list(variance = NA_real_,
                note = "Some sampled rows have no cluster label."))
  }
  u <- rowsum(y / pi_i, code, reorder = FALSE)[, 1L]
  m <- length(u)
  if (m < 2L) {
    return(list(variance = NA_real_, note = paste0(
      "A variance needs at least two clusters; this sample has ", m, ".")))
  }
  fpc <- if (is.finite(A) && A > m) 1 - m / A else 0
  list(variance = fpc * m * stats::var(u), note = NULL)
}

#' Two-stage variance, declining where a stage is unmeasurable
#'
#' With an equal take that divides `n`, drawn without replacement, the exact
#' estimator is `two_stage_variance()`. Any other multistage design has no
#' closed-form joint probabilities, and the default method declines.
#'
#' A sample from a single cluster has no variation between clusters to
#' measure. Sen-Yates-Grundy came out negative on every such sample, or 0 for
#' a constant `y`, so this declines first, in the words `design_cluster()`
#' uses. A sample whose every row is certain goes on to the exact estimator,
#' which gives 0.
#'
#' With one row per selected cluster, pairs inside a cluster can never be drawn
#' together, so the exact estimator has no way to see the within-cluster
#' variation. With `f1 = n_clusters / A`, it understates by `f1` times the
#' within-cluster component of the variance. The delete-a-cluster jackknife
#' that `variance = "auto"` falls back to is the ultimate-cluster
#' approximation. It has no first-stage finite population correction, so it
#' overstates by `f1` times the between-cluster component instead, which
#' `fallback_note()` reports once `f1` passes 0.2.
#'
#' @noRd
ht_variance.drawn_design_multistage <- function(design, data, rows, y, pi_i,
                                                domain = NULL, ctx = NULL) {
  if (!is.null(design$size)) {
    if (design$replace) return(NextMethod())
    return(pps_stage_variance(design, data, rows, y, pi_i, ctx))
  }
  cl <- from_context(ctx, "clusters", count_clusters(design, data))
  n_cl <- length(unique(cl$code[rows]))
  if (n_cl < 2L && !all(is_certain(pi_i))) {
    return(list(variance = NA_real_, note = paste0(
      "A variance needs at least two clusters; this sample has ", n_cl, ".")))
  }
  if (!design$replace && design$allocation == "equal" &&
      design$n %/% design$n_clusters == 1L &&
      design$n %% design$n_clusters == 0L) {
    return(list(variance = NA_real_, note = paste0(
      "With one row per selected cluster the variation within clusters ",
      "cannot be measured, and the exact two-stage estimator would understate ",
      "the variance.", fallback_note(design$n_clusters, cl$total))))
  }
  if (design$replace || design$allocation != "equal" ||
      design$n %% design$n_clusters != 0L) {
    return(NextMethod())
  }
  v <- two_stage_variance(y, cl$code[rows], design$n_clusters, cl$total,
                          design$n %/% design$n_clusters,
                          from_context(ctx, "cluster_sizes",
                                       tabulate(cl$code, cl$total)))
  list(variance = v, note = NULL)
}

#' Two-stage variance with clusters selected in proportion to size
#'
#' The unbiased two-stage estimator is the first-stage estimator applied to
#' the estimated cluster totals, plus `sum(V_i / pi_i)` over the selected
#' clusters, where `V_i = N_i^2 * (1 - m / N_i) * s_i^2 / m` is the variance
#' of the simple random sample of `m` rows within cluster `i` (Särndal,
#' Swensson and Wretman 1992). Systematic PPS has no joint
#' probabilities for the first stage, so Deville's approximation takes its
#' place, over the clusters below certainty; a cluster taken with certainty
#' adds only its own `V_i`. With every `pi_i` equal this is
#' `two_stage_variance()`.
#'
#' One row per cluster leaves `V_i` unmeasurable, and one cluster below
#' certainty leaves nothing to compare between clusters, so both decline.
#'
#' @noRd
pps_stage_variance <- function(design, data, rows, y, pi_i, ctx = NULL) {
  if (all(is_certain(pi_i))) return(all_certain())
  st <- from_context(ctx, "stage_one", cluster_stage(design, data))
  code <- st$code[rows]
  live <- unique(code)
  g <- match(code, live)
  k <- length(live)
  p1 <- st$pi1[live]
  N <- as.numeric(st$sizes[live])
  m <- tabulate(g, k)
  if (any(m < 2L & m < N)) {
    return(list(variance = NA_real_, note = paste0(
      "With one row per selected cluster the variation within clusters ",
      "cannot be measured, and leaving it out would understate the ",
      "variance.")))
  }
  open <- which(!is_certain(p1))
  if (length(open) == 1L) {
    return(list(variance = NA_real_, note = paste0(
      "A variance needs at least two clusters below certainty; this sample ",
      "has one, besides the ", k - 1L, " taken with certainty.")))
  }
  z <- as.numeric(y)
  mu <- rowsum(z, g)[, 1L] / m
  ss <- rowsum((z - mu[g])^2, g)[, 1L]
  use <- m >= 2L & m < N
  within <- sum((N^2 * (1 - m / N) * ss / ((m - 1) * m) / p1)[use])
  if (!length(open)) return(list(variance = within, note = NULL))
  # Each cluster's estimated total over its stage-one probability.
  u <- rowsum(z / pi_i, g)[, 1L]
  between <- deville_variance(u[open] * p1[open], p1[open])$variance
  list(variance = between + within, method = "deville",
       note = paste0("Deville's approximation for the clusters selected in ",
                     "proportion to size, which have no closed-form joint ",
                     "inclusion probabilities, plus the variance within the ",
                     "selected clusters."))
}

#' What the ultimate-cluster jackknife does to a multistage variance
#'
#' The jackknife has no first-stage finite population correction, so its
#' between-cluster part is `1 / (1 - f1)` times the true one, with
#' `f1 = a / A`. That is under 1.25 while `f1` is at most 0.2 and is not
#' reported. When every cluster is selected the true between-cluster part is
#' 0, and the clusters are strata. Returns a sentence to append to a note,
#' starting with a space, or `""`.
#'
#' @noRd
fallback_note <- function(a, A) {
  f1 <- a / A
  if (!is.finite(f1) || f1 <= 0.2) return("")
  if (a >= A) {
    return(paste0(
      " Every cluster is selected, so the clusters are strata and the ",
      "variance between them is zero by design, but the delete-a-cluster ",
      "jackknife counts it in full because it has no first-stage finite ",
      "population correction. Draw with design_stratified() using the ",
      "cluster column as `strata` and at least two rows in each."))
  }
  paste0(
    " With ", a, " of the ", A, " clusters selected, the delete-a-cluster ",
    "jackknife, which has no first-stage finite population correction, puts ",
    "the variance between clusters at ", signif(1 / (1 - f1), 2),
    " times its true size.")
}

#' Variance for the weighted designs
#'
#' Poisson sampling is exact: rows are independent, so no joint matrix is
#' needed. Systematic PPS has no closed-form joint probabilities; Deville's
#' approximation is the standard answer for a high-entropy fixed-size design.
#'
#' @noRd
ht_variance.drawn_design_weighted <- function(design, data, rows, y, pi_i,
                                              domain = NULL, ctx = NULL) {
  if (design$method == "poisson") {
    return(list(variance = sum((1 - pi_i) / pi_i^2 * y^2), note = NULL))
  }
  if (design$method == "systematic") {
    if (!is.null(design$strata)) {
      return(deville_strata(design, data, rows, y, pi_i, domain, ctx))
    }
    return(deville_variance(y, pi_i))
  }
  stop_declined(paste0("`design_weighted(method = \"successive\")` has no ",
                       "closed-form inclusion probability."))
}

#' Deville's (1999) variance approximation for fixed-size unequal-probability
#' sampling
#'
#' `sum(c_k * (y_k / pi_k - A)^2) / (1 - sum(a_k^2))` with `c_k = 1 - pi_k`,
#' `a_k = c_k / sum(c)` and `A = sum(a_k * y_k / pi_k)`: the Hajek
#' approximation to the joint probabilities, which Matei and Tillé (2005) found
#' the most reliable of the standard first-order-only estimators. Certainty
#' rows have `c_k = 0` and drop out, as they should. Identical to
#' `sampling::varest()`.
#'
#' @noRd
deville_variance <- function(y, pi_i) {
  c_k <- ifelse(is_certain(pi_i), 0, 1 - pi_i)
  if (sum(c_k) <= 0) return(all_certain("deville"))
  a <- c_k / sum(c_k)
  denom <- 1 - sum(a^2)
  if (denom <= 0) {
    return(list(variance = NA_real_, method = "deville", note = paste0(
      "Deville's approximation needs at least two sampled rows below ",
      "certainty.")))
  }
  u <- y / pi_i
  A <- sum(a * u)
  list(variance = sum(c_k * (u - A)^2) / denom, method = "deville",
       note = paste0("Deville's approximation for unequal-probability ",
                     "sampling: systematic PPS has no closed-form joint ",
                     "inclusion probabilities."))
}

#' Deville's approximation for systematic PPS within strata
#'
#' The strata are drawn independently, so the variance is the sum of
#' `deville_variance()` over each stratum's sampled rows, which is
#' `sampling::varest()` applied stratum by stratum and added up. Rows taken
#' with certainty add nothing. A stratum with one sampled row below certainty
#' has nothing to compare it with, so it is refused by `lonely_note()`,
#' unless `domain` shows it holds no unit of the domain, when it adds 0.
#'
#' @noRd
deville_strata <- function(design, data, rows, y, pi_i, domain = NULL,
                           ctx = NULL) {
  free <- !is_certain(pi_i)
  if (!any(free)) return(all_certain("deville"))
  st <- from_context(ctx, "weighted", weighted_setup(design, data))
  p <- from_context(ctx, "pps_pi", weighted_pi(design, st))
  # Each frame row's stratum, counting only the rows below certainty.
  code <- rep(NA_integer_, nrow(data))
  open <- !is_certain(p)
  code[st$keep[open]] <- st$code[open]
  group <- structure(code, levels = st$labels, class = "factor")
  note <- lonely_note(group, rows, c("stratum", "strata"), lonely_fix(design),
                      domain)
  if (!is.null(note)) {
    return(list(variance = NA_real_, method = "deville", note = note))
  }
  h <- code[rows]
  v <- 0
  for (g in unique(h[free])) {
    i <- which(free & h == g)
    if (length(i) < 2L) next
    v <- v + deville_variance(y[i], pi_i[i])$variance
  }
  list(variance = v, method = "deville",
       note = paste0("Deville's approximation for unequal-probability ",
                     "sampling, summed over strata: systematic PPS has no ",
                     "closed-form joint inclusion probabilities."))
}

#' Successive-difference variance for a systematic sample
#'
#' Each sampled row is compared with its neighbour in the order the design
#' walked, and the squared differences stand in for the population variance:
#' `(1 - f) * n / (2 * (n - 1)) * sum(diff(y / pi)^2)`. A trend along the walk
#' order therefore counts as structure, not noise -- the reason systematic
#' sampling on a sorted frame is efficient, and the reason the simple random
#' formula overstates its variance there.
#'
#' When the interval `k` does not divide the `N` rows the design walks, the
#' sample holds `N %/% k` or `N %/% k + 1` rows, the larger with probability
#' `(N %% k) / k`. The total `sum(y / pi)` moves with that count by about
#' `mean(y / pi)` per row, which successive differences cannot see. A
#' constant `y` got a variance of 0, and a variable with a coefficient of
#' variation of 0.02 an interval that covered none of the time. The term
#' `mean(y / pi)^2 * (r / k) * (1 - r / k)`, with `r` that remainder, is added
#' for it. It is exact for a constant `y`, ignores the covariance between the
#' count and the sample mean, and is 0 for the linearised residuals of a Hajek
#' mean, whose mean is 0.
#'
#' A design given `n` walks at a step of `N / n` rows and always draws `n`,
#' so `f` is `n / N` and no term is added.
#'
#' A design with `n_starts` of 2 or more has the exact random group variance
#' of `replicate_variance()` instead.
#'
#' @noRd
ht_variance.drawn_design_systematic <- function(design, data, rows, y, pi_i,
                                                domain = NULL, ctx = NULL) {
  refuse_fixed_start(design)
  if (replicated(design)) {
    walk <- from_context(ctx, "walk", systematic_positions(design, data))
    return(replicate_variance(design, walk, rows, y, pi_i))
  }
  n <- length(rows)
  if (n < 2L) {
    return(list(variance = NA_real_,
                note = "A variance needs at least two sampled rows."))
  }
  walk <- from_context(ctx, "walk", systematic_positions(design, data))
  pos <- walk[rows]
  u <- (y / pi_i)[order(pos)]
  N <- sum(!is.na(walk))
  if (is.null(design$n)) {
    k <- design$interval
    period <- "`interval`"
    r <- N %% k
  } else {
    # A step of N / n rows, which need not be whole, and a fixed size.
    k <- N / design$n
    period <- paste0("the step of ", format(k, digits = 4), " rows")
    r <- 0
  }
  f <- 1 / k
  v <- (1 - f) * n / (2 * (n - 1)) * sum(diff(u)^2)
  if (r > 0) v <- v + mean(u)^2 * (r / k) * (1 - r / k)
  list(variance = v,
       method = "successive difference",
       note = paste0("Systematic sampling has no design-unbiased variance ",
                     "estimator; this is the successive-difference ",
                     "approximation. It would be understated if the frame ",
                     "cycled with a period matching ", period, ", or if a ",
                     "smooth trend along the sort order swamped all other ",
                     "variation.",
                     if (r > 0) {
                       paste0(" The ", N, " rows walked are not a multiple of ",
                              "`interval` (", k, "), so the sample size is ",
                              "random, ", N %/% k, " or ", N %/% k + 1L,
                              " rows, and the variance includes a term for ",
                              "that.")
                     }))
}

#' Random group variance of a systematic sample with several starts
#'
#' The `m` starts are a simple random sample of the `m * k` residue classes
#' modulo `m * k`, `k` the interval, and each replicate's total `t_r`, `m`
#' times the sum of `y / pi` over its rows, is the Horvitz-Thompson total of
#' a systematic sample at a step of `m * k`. Their mean is the sample's
#' total, and `(1 - 1 / k) * var(t_r) / m` is the textbook unbiased variance
#' of a mean of `m` units drawn from `m * k` without replacement. Over the
#' sampled rows it is the Horvitz-Thompson variance estimator with the exact
#' joint probabilities of `joint_inclusion()`, and, when `m * k` divides the
#' rows walked and the size is fixed, the Sen-Yates-Grundy one too, computed
#' from the `m` replicate totals alone. A replicate whose start fell past the
#' last row walked drew nothing, and its total of 0 counts among the `m`.
#'
#' @noRd
replicate_variance <- function(design, walk, rows, y, pi_i) {
  if (all(is_certain(pi_i))) return(all_certain())
  m <- systematic_starts(design)
  rep_id <- systematic_replicates(design, walk, rows)
  t_r <- m * rowsum(as.numeric(y / pi_i), rep_id, reorder = FALSE)[, 1L]
  t_r <- c(t_r, numeric(m - length(t_r)))
  list(variance = (1 - 1 / design$interval) * stats::var(t_r) / m,
       note = NULL)
}

# ---- jackknife --------------------------------------------------------------

#' Delete-a-group jackknife variance
#'
#' Works from the sampled rows alone, which is what makes it available where
#' the analytic form is not. `units` is the frame's copy of those rows, so the
#' strata and clusters are the ones the sample was drawn with. Primary sampling units are whole clusters where the
#' design has them, otherwise individual rows; where the design has strata --
#' stratified, and temporal, whose intervals are strata -- units are deleted
#' within their own stratum and only that stratum is reweighted (JKn). Deleting
#' across strata instead counts the differences *between* strata as sampling
#' variance, which stratification exists to remove: on a frame with four
#' well-separated strata that overstated the standard error sixty-fold.
#'
#' Two designs are declined rather than approximated. A systematic sample has a
#' single primary sampling unit -- the random start -- so deleting rows does not
#' reflect its randomness at all. Poisson sampling has an exact variance and a
#' random size; deleting rows from it understates by a factor of three.
#'
#' Units taken with certainty are held fixed. They are in every possible
#' sample, so deleting one and inflating the rest invents variance the design
#' does not have. On a frame with six dominant units, systematic PPS samples
#' got 181 times the true variance on average that way. The rule covers the
#' units `design_weighted()` and `design_spread(size =)` take with certainty,
#' strata taken whole and a census. A cluster is held fixed when every one of
#' its sampled rows is certain. The finite population correction counts only
#' the units that were not certain, in the sample and in the frame.
#'
#' A stratum with one unit cannot be deleted from and is refused, unless
#' `domain` shows it holds no unit of the domain being estimated, as in
#' `lonely_note()`.
#'
#' `ctx` is a `variance_context()` for `design`, `pop` and `units`, with
#' `pi_i` the same in every call that shares it.
#'
#' @noRd
jackknife_variance <- function(design, units, pop, y, pi_i, domain = NULL,
                               ctx = NULL) {
  refusal <- jackknife_refusal(design)
  if (!is.null(refusal)) {
    return(list(variance = NA_real_, method = "jackknife", note = refusal))
  }
  if (inherits(design, "drawn_design_certainty")) {
    return(jackknife_certainty(design, units, pop, y, pi_i, domain, ctx))
  }

  psu <- from_context(ctx, "psu", jackknife_groups(design, units, pop, ctx))
  fixed <- !(psu %in% psu[!is_certain(pi_i)])
  if (all(fixed)) return(all_certain("jackknife"))
  n_fixed <- length(unique(psu[fixed]))

  st <- jackknife_strata(design, units, pop, psu, fixed, domain, ctx)
  base <- y / pi_i
  total <- sum(base)
  free <- which(!fixed)
  by_stratum <- split(free, factor(match(st$stratum[free], st$levels),
                                   levels = seq_along(st$levels)))

  v <- 0
  n_psu <- 0L
  lonely <- integer(0)
  for (k in seq_along(st$levels)) {
    in_h <- by_stratum[[k]]
    # Each unit's total in the stratum. The replicates are symmetric in the
    # units, so their order does not matter.
    sums <- rowsum(base[in_h], psu[in_h], reorder = FALSE)[, 1L]
    m <- length(sums)
    f <- st$fpc[k]
    n_psu <- n_psu + m
    if (m < 2L) {
      if (f > 0 && st$in_domain[k]) lonely <- c(lonely, k)
      next
    }
    reps <- (total - sum(sums)) + (sum(sums) - sums) * m / (m - 1)
    v <- v + f * (m - 1) / m * sum((reps - mean(reps))^2)
  }

  n_strata <- length(st$levels)
  if (n_strata == 1L && n_psu < 2L) {
    return(list(variance = NA_real_, method = "jackknife",
                note = paste0("The jackknife needs at least two primary ",
                              "sampling units; this sample has ", n_psu,
                              if (n_fixed > 0L) {
                                paste0(" besides the ", n_fixed, " taken ",
                                       "with certainty")
                              }, ". Draw more units below certainty.")))
  }
  if (length(lonely)) {
    if (!is.null(st$between) && st$between %in% st$levels[lonely]) {
      return(list(variance = NA_real_, method = "jackknife", note = paste0(
        "The jackknife needs at least two clusters below certainty; this ",
        "sample has one.")))
    }
    noun <- if (inherits(design, "drawn_design_temporal")) {
      c("interval", "intervals")
    } else if (!is.null(st$between)) {
      c("cluster taken with certainty", "clusters taken with certainty")
    } else {
      c("stratum", "strata")
    }
    lonely <- lonely[order(st$levels[lonely])]
    return(list(variance = NA_real_, method = "jackknife",
                note = lonely_text(st$labels[lonely], noun,
                                   lonely_fix(design), domain)))
  }
  if (!is.finite(v)) {
    return(list(variance = NA_real_, method = "jackknife",
                note = overflow_note()))
  }
  list(variance = v, method = "jackknife",
       note = paste0("Jackknife over ", n_psu, " primary sampling unit(s)",
                     if (n_strata > 1L) paste0(" in ", n_strata, " strata")
                     else "",
                     if (any(st$fpc < 1)) ", with a finite population correction"
                     else "",
                     if (n_fixed > 0L) {
                       paste0("; the ", n_fixed, " taken with certainty add ",
                              "no variance and were held fixed")
                     }, "."))
}

#' The jackknife for a certainty design
#'
#' Rows at or above the threshold are in every sample, so the whole problem is
#' handed to `rest`, over the rows below the threshold. Those are found from
#' the frame's copy of the sampled rows, which holds the `above` column. Rows
#' that `rest` itself takes with certainty go along and are held fixed there,
#' so its finite population correction counts them correctly.
#'
#' @noRd
jackknife_certainty <- function(design, units, pop, y, pi_i, domain = NULL,
                                ctx = NULL) {
  sp <- certainty_parts(design, pop, ctx)
  v <- units[[design$above]]
  below <- !is.na(v) & v < design$threshold
  if (!any(below)) {
    if (length(sp$rest)) {
      out <- empty_rest(design, length(sp$rest))
      out$method <- "jackknife"
      return(out)
    }
    return(all_certain("jackknife"))
  }
  jackknife_variance(design$rest, units[below, , drop = FALSE], sp$rest_data,
                     y[below], pi_i[below], rest_domain(domain, sp),
                     sub_context(ctx, "rest"))
}

#' Designs the jackknife should decline rather than approximate
#'
#' @noRd
jackknife_refusal <- function(design) {
  if (replicated(design)) {
    return(paste0("A systematic sample with several starts has the starts as ",
                  "its primary sampling units, not the rows, so deleting ",
                  "rows does not reflect the design's randomness and the ",
                  "jackknife is not valid here. variance = \"analytic\" ",
                  "gives the exact variance between the replicates."))
  }
  if (inherits(design, "drawn_design_systematic")) {
    return(paste0("A systematic sample has one primary sampling unit -- the ",
                  "random start -- so deleting rows does not reflect the ",
                  "design's randomness and the jackknife is not valid here."))
  }
  if (walks_strata(design)) {
    return(paste0("Each stratum of a systematic sample within strata has one ",
                  "random start, so deleting rows does not reflect the ",
                  "design's randomness and the jackknife is not valid here."))
  }
  if (inherits(design, "drawn_design_weighted") && design$method == "poisson") {
    return(paste0("Poisson sampling selects rows independently and has an ",
                  "exact variance, sum((1 - pi) / pi^2 * y^2). Deleting rows ",
                  "from a sample whose size is itself random understates it ",
                  "by a factor of three. Use variance = \"analytic\"."))
  }
  NULL
}

#' Strata for the jackknife, and each stratum's finite population correction
#'
#' Stratified and temporal designs, and spread and weighted designs with
#' `strata`, delete within strata, with each stratum's own `1 - n_h / N_h`.
#' The strata of a spread or weighted design count only the rows it can
#' draw. Everything else is a single stratum, corrected by
#' `jackknife_fpc()`. Only strata holding a unit that is not `fixed` are
#' listed, and both counts leave the fixed units out. A unit taken with
#' certainty is in every sample, so the fixed units in the frame are the
#' ones in the sample. A stratum taken whole drops out.
#'
#' A multistage design with `size` deletes clusters below certainty in one
#' stratum, `between`, with no correction, and the rows of each cluster
#' taken with certainty within a stratum of its own, as `cluster_units()`
#' lays them out.
#'
#' Returns each row's `stratum`, the `levels` listed, their `labels`, `fpc`
#' and `in_domain`, in the same order, and `between`. The levels are integer
#' codes in key order, and the labels the key's levels, as the analytic
#' estimators name the same strata and intervals. `in_domain` is `FALSE` for
#' a stratum that `domain`, a logical vector over the frame, shows to hold
#' no unit of the domain, and `TRUE` otherwise.
#'
#' @noRd
jackknife_strata <- function(design, units, pop, psu, fixed, domain = NULL,
                             ctx = NULL) {
  # Integer codes, shared by the sampled rows and the frame.
  grp <- from_context(ctx, "strata", switch(design_type(design),
    stratified = ,
    spread = ,
    weighted = if (!is.null(design$strata)) {
      key <- group_key2(units, pop, design$strata)
      pop_code <- as.integer(key$table)
      # Rows the design drops have no stratum to count in.
      keep <- switch(design_type(design),
        spread = from_context(ctx, "spread", spread_frame(design, pop))$keep,
        weighted = from_context(ctx, "weighted",
                                weighted_setup(design, pop))$keep,
        NULL
      )
      if (!is.null(keep)) {
        out <- rep(TRUE, nrow(pop))
        out[keep] <- FALSE
        pop_code[out] <- NA_integer_
      }
      list(sample = as.integer(key$x), pop = pop_code,
           labels = levels(key$x))
    },
    multistage = if (!is.null(design$size)) {
      # Clusters below certainty share stratum 1, with no finite population
      # correction, as for any multistage design; each cluster taken with
      # certainty is a stratum of its own, of rows drawn without replacement.
      cu <- cluster_units(design, units, pop, ctx)
      st <- from_context(ctx, "stage_one", cluster_stage(design, pop))
      first <- match(seq_along(st$pi1), st$code)
      list(sample = cu$stratum, pop = cu$pop_stratum,
           labels = c("(below certainty)",
                      as.character(pop[[design$clusters]][first])),
           between = 1L)
    },
    temporal = {
      # Keyed as srs_groups() keys the intervals, so the two name them alike.
      key <- key_from_columns(list(c(temporal_bucket(design, units),
                                     temporal_bucket(design, pop))))
      code <- as.integer(key)
      list(sample = code[seq_len(nrow(units))],
           pop = code[nrow(units) + seq_len(nrow(pop))],
           labels = levels(key))
    },
    NULL
  ))
  if (is.null(grp)) {
    m <- length(unique(psu[!fixed]))
    n_fixed <- length(unique(psu[fixed]))
    return(list(stratum = rep("all", nrow(units)), levels = "all",
                labels = "all",
                fpc = jackknife_fpc(design, pop, m, n_fixed, ctx),
                in_domain = TRUE))
  }
  lev <- unique(grp$sample[!fixed])
  N_h <- tabulate(match(grp$pop, lev), length(lev))
  n_h <- tabulate(match(grp$sample[!fixed], lev), length(lev))
  free_h <- N_h - tabulate(match(grp$sample[fixed], lev), length(lev))
  in_domain <- if (known_domain(domain, length(grp$pop))) {
    tabulate(match(grp$pop[domain], lev), length(lev)) > 0L
  } else {
    rep(TRUE, length(lev))
  }
  fpc <- ifelse(free_h > n_h, 1 - n_h / free_h, 0)
  if (!is.null(grp$between)) fpc[lev == grp$between] <- 1
  list(stratum = grp$sample, levels = lev, labels = grp$labels[lev],
       fpc = fpc, in_domain = in_domain, between = grp$between)
}

#' The share of primary sampling units NOT taken
#'
#' Applied for single-stage designs, where deleting a sampling unit accounts for
#' all the randomness there is. That includes probability-proportional-to-size
#' selection, which is still a fixed-size draw without replacement: measured
#' against the empirical sampling variance over four frames and two response
#' variables, the corrected estimator sits between 0.86 and 1.11 of the truth
#' where the uncorrected one runs from 0.98 to 1.77, overstating badly once the
#' sampling fraction gets large.
#'
#' `m` counts the sampled units that were not certain, and `n_fixed` the
#' certain ones, which come off the frame's count. A census therefore gets 0.
#' Without a finite count of units in the frame there is nothing to correct
#' by, and the correction is 1.
#'
#' Not applied to multi-stage designs. Each estimated cluster total carries
#' the second stage's variance as well as the first's, and multiplying by the
#' first-stage correction `1 - f1`, with `f1 = n_clusters / A`, would scale
#' both, understating by `f1` times the within-cluster component. Leaving it
#' out, the ultimate-cluster approximation usual in survey practice, overstates
#' by `f1` times the between-cluster component instead, which errs on the
#' side of a wider interval.
#'
#' @noRd
jackknife_fpc <- function(design, pop, m, n_fixed = 0L, ctx = NULL) {
  if (inherits(design, "drawn_design_multistage")) return(1)
  total <- if (inherits(design, "drawn_design_cluster")) {
    from_context(ctx, "clusters", count_clusters(design, pop))$total
  } else {
    from_context(ctx, "reachable", reachable_rows(design, pop))
  }
  if (!is.finite(total)) return(1)
  free <- total - n_fixed
  if (free <= m) return(0)
  1 - m / free
}

#' Rows of the frame a design can select
#'
#' The population a finite population correction counts. Rows a design drops,
#' through `na_rm`, past `max_items` or outside its region or window, have
#' probability 0 and are not part of it; counting them put the correction
#' for a systematic or PPS design with `na_rm` too close to 1. The count is
#' read from the design's own columns, so it needs no inclusion probabilities
#' and is defined for a design without a closed form, such as
#' `design_weighted(method = "successive")`. A stratum allocated no rows is
#' not seen here; `sum(pi > 0)` over the draw-time probabilities counts those
#' out too, where a sample stores them.
#'
#' @noRd
reachable_rows <- function(design, pop) {
  n <- nrow(pop)
  keyed <- function(cols) {
    sum(!Reduce(`|`, lapply(cols, function(col) is.na(pop[[col]]))))
  }
  switch(design_type(design),
    reservoir = reservoir_reach(design, n),
    spatial = sum(spatial_inside(design, pop)),
    temporal = sum(!is.na(temporal_bucket(design, pop))),
    systematic = sum(!is.na(systematic_positions(design, pop))),
    weighted = keyed(c(design$weights, design$strata)),
    spread = {
      bad <- spread_missing(
        coordinate_columns(pop, design$across, "across", time = TRUE),
        if (!is.null(design$size)) pop[[design$size]])
      if (!is.null(design$strata)) {
        bad <- bad | is.na(group_key(pop, design$strata))
      }
      sum(!bad)
    },
    stratified = keyed(c(design$strata, design$within$order_by)),
    cluster = ,
    multistage = keyed(c(design$clusters, design$size)),
    certainty = {
      sp <- certainty_split(design, pop)
      length(sp$take) +
        reachable_rows(design$rest, sp$data[sp$rest, , drop = FALSE])
    },
    n
  )
}

#' Primary sampling units for the jackknife
#'
#' Whole clusters for a cluster or multistage design, read from the frame's
#' copy of the sampled rows, which always has the cluster column, as integer
#' codes; otherwise each row is its own unit. A multistage design with
#' `size` reads which clusters were taken with certainty from the frame
#' `pop`, and makes each of their sampled rows a unit, as `cluster_units()`
#' does.
#'
#' @noRd
jackknife_groups <- function(design, units, pop = NULL, ctx = NULL) {
  if (inherits(design, "drawn_design_multistage") && !is.null(design$size)) {
    return(cluster_units(design, units, pop, ctx)$psu)
  }
  if (inherits(design, c("drawn_design_cluster", "drawn_design_multistage"))) {
    return(cluster_codes(units[[design$clusters]]))
  }
  seq_len(nrow(units))
}

#' @export
print.drawn_ht <- function(x, ...) {
  print_estimate(x, "total", "Horvitz-Thompson total")
}
