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
#'     known joint probabilities use the Sen-Yates-Grundy estimator; Poisson
#'     sampling, whose rows are independent, uses `sum((1 - pi) / pi^2 * y^2)`;
#'     cluster designs apply the estimator at the cluster level, because taking
#'     whole clusters makes the row count random; certainty designs hand the
#'     problem to `rest`, since certainty rows are in every sample and add
#'     nothing. Sen-Yates-Grundy works on matrices with a row and a column per
#'     sampled row, and past about 4,700 sampled rows they would need more
#'     than 2 GB, so it declines and says how much.}
#'   \item{`"deville"`}{Systematic probability-proportional-to-size
#'     ([design_weighted()] with `method = "systematic"`). Its joint
#'     probabilities have no closed form, so this uses Deville's (1999)
#'     approximation, which needs only the first-order probabilities and is the
#'     best-performing of the standard approximations for high-entropy designs
#'     (Matei and Tillé 2005). Randomised systematic selection is close to high
#'     entropy (Hartley and Rao 1962). Unlike a plain jackknife, each row's
#'     contribution is scaled by its own `1 - pi`, so rows near certainty stop
#'     adding variance they do not have.}
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
#'     matching `interval`, and when a smooth trend so dominates `y` that the
#'     random start is almost the only source of variation, it runs low (0.3 to
#'     0.7 of the truth when the trend spanned 25 noise standard deviations).
#'     When `interval` does not divide the number of rows walked, the sample
#'     size is random, one of two values, and the total moves with it in a way
#'     successive differences cannot see. A term for that,
#'     `mean(y / pi)^2 * (r / k) * (1 - r / k)` with `k` the interval and `r`
#'     the remainder, is added, and the note says so. Without it a constant
#'     `y` had a standard error of 0, and a variable with a coefficient of
#'     variation of 0.02 an interval that never covered the truth.}
#'   \item{`"local mean"`}{[design_spread()]. A spatially balanced sample is
#'     well spread precisely because neighbours are rarely taken together, so
#'     this compares each sampled row with the local mean of its nearest
#'     sampled neighbours (Stevens and Olsen 2003; Grafström and Schelin 2014).}
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
#' @seealso [ht_mean()], [deff()], [joint_prob()], [inclusion_prob()]
#' @export
ht_total <- function(sample, y, variance = c("auto", "analytic", "jackknife",
                                            "none"), level = 0.95,
                     by = NULL, df = NULL,
                     alternative = c("two.sided", "less", "greater")) {
  variance <- check_choice(variance,
                           c("auto", "analytic", "jackknife", "none"),
                           "variance")
  alternative <- check_choice(alternative, c("two.sided", "less", "greater"),
                              "alternative")
  parts <- ht_prepare(sample, y, level, df)
  if (!is.null(by)) {
    return(estimate_by(parts, by, variance, level, df, what = "total",
                       alternative = alternative))
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
#' @noRd
ht_variance_dispatch <- function(design, units, pop, rows, z, pi_i, variance,
                                 domain = NULL) {
  if (variance == "none") {
    return(list(variance = NA_real_, method = "none",
                note = "Variance not requested."))
  }
  if (variance == "jackknife") {
    jk <- jackknife_variance(design, units, pop, z, pi_i, domain)
    if (is.na(jk$variance)) jk$method <- "none"
    return(jk)
  }

  # A negative Sen-Yates-Grundy estimate, a declining estimator and an
  # overflow are failures rather than variances. Each gets a reason instead
  # of a number that sqrt() and deff() would then have to guess about.
  a <- tryCatch(ht_variance(design, pop, rows, z, pi_i, domain),
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

  jk <- jackknife_variance(design, units, pop, z, pi_i, domain)
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
                          psu = sample_units(parts$design, parts$units)) {
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

#' Each sampled row's primary sampling unit, as the jackknife deletes them
#'
#' Clusters for a cluster or multistage design, otherwise rows. Under a
#' certainty design each row at or above the threshold is a unit of its own
#' and the rows below it are grouped as `rest` groups them.
#'
#' @noRd
sample_units <- function(design, units) {
  if (!inherits(design, "drawn_design_certainty")) {
    return(jackknife_groups(design, units))
  }
  v <- units[[design$above]]
  below <- !is.na(v) & v < design$threshold
  out <- -seq_len(nrow(units))
  out[below] <- sample_units(design$rest, units[below, , drop = FALSE])
  out
}

#' What a design's primary sampling unit is called in a message
#' @noRd
unit_noun <- function(design) {
  clustered <- c("drawn_design_cluster", "drawn_design_multistage")
  if (inherits(design, clustered)) return("cluster")
  if (inherits(design, "drawn_design_certainty") &&
      inherits(design$rest, clustered)) {
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
#' the number of clusters less one.
#'
#' With no unit below certainty the variance is exactly 0 and the result is
#' `Inf`. An empty sample, and a certainty design whose `rest` drew none of the
#' rows below the threshold, have no variance to build an interval on and get
#' 0.
#'
#' `as_svydesign()` hands each group of certain units to `survey` as one
#' primary sampling unit in a stratum of its own, so `survey::degf()` gives the
#' same count. Strata and clusters are read from `units`, the frame's copy of
#' the sampled rows; `pi_i` holds their inclusion probabilities.
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
  psu <- switch(design_type(design),
    cluster = ,
    multistage = cluster_codes(units[[design$clusters]]),
    seq_len(n)
  )
  stratum <- switch(design_type(design),
    stratified = as.integer(group_key(units, design$strata)),
    temporal = temporal_bucket(design, units),
    rep("all", n)
  )
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
ht_variance <- function(design, data, rows, y, pi_i,
                        domain = NULL) UseMethod("ht_variance")

#' Sen-Yates-Grundy, for fixed-size designs with closed-form joint
#' probabilities
#'
#' The joint probabilities and the sums over pairs are n x n matrices of
#' doubles. A stratified sample held about eleven of them at once, 1.4 GB at
#' 4,000 rows, so past 2 GiB projected the estimator declines with the size it
#' would need. `variance = "auto"` then uses the jackknife, which gives the
#' same figure for a simple random or stratified sample in time linear in the
#' rows.
#'
#' @noRd
ht_variance.default <- function(design, data, rows, y, pi_i, domain = NULL) {
  n <- length(rows)
  bytes <- 12 * 8 * as.numeric(n)^2
  if (bytes > 2^31) {
    stop_declined(paste0(
      "The Sen-Yates-Grundy estimator would need about ",
      format(signif(bytes / 1e9, 2)), " GB of memory for its ", n, " x ", n,
      " matrices.\nUse variance = \"auto\" to fall back to the jackknife."))
  }
  pij <- joint_inclusion(design, data, rows)
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
  N_h <- tabulate(code, k)
  n_h <- tabulate(code[rows], k)
  lonely <- n_h == 1L & N_h > 1L
  if (known_domain(domain, length(code))) {
    lonely <- lonely & tabulate(code[domain], k) > 0L
  }
  lonely <- levels(group)[lonely]
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

#' Does `domain` mark which of `n` frame rows belong to a domain?
#' @noRd
known_domain <- function(domain, n) {
  is.logical(domain) && length(domain) == n && !anyNA(domain)
}

#' @noRd
ht_variance.drawn_design_stratified <- function(design, data, rows, y, pi_i,
                                                domain = NULL) {
  group <- group_key(data, design$strata)
  note <- lonely_note(group, rows, c("stratum", "strata"),
                      "Draw at least two rows per stratum (min_per_stratum = 2).",
                      domain)
  if (!is.null(note)) return(list(variance = NA_real_, note = note))
  NextMethod()
}

#' @noRd
ht_variance.drawn_design_temporal <- function(design, data, rows, y, pi_i,
                                              domain = NULL) {
  bucket <- temporal_bucket(design, data)
  note <- lonely_note(bucket, rows, c("interval", "intervals"),
                      "Draw at least two rows per interval (per_interval = 2).",
                      domain)
  if (!is.null(note)) return(list(variance = NA_real_, note = note))
  NextMethod()
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
                                             domain = NULL) {
  cl <- count_clusters(design, data)
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
#' A sample from a single cluster has no variation between clusters to
#' measure. Sen-Yates-Grundy comes out negative on every such sample, or 0 for
#' a constant `y`, so this declines first, in the words `design_cluster()`
#' uses. A sample whose every row is certain goes on to the exact estimator,
#' which gives 0.
#'
#' With one row per selected cluster, pairs inside a cluster can never be drawn
#' together, so Sen-Yates-Grundy has no way to see the within-cluster
#' variation. With `f1 = n_clusters / A`, it understates by `f1` times the
#' within-cluster component of the variance. The delete-a-cluster jackknife
#' that `variance = "auto"` falls back to is the ultimate-cluster
#' approximation. It has no first-stage finite population correction, so it
#' overstates by `f1` times the between-cluster component instead, which
#' `fallback_note()` reports once `f1` passes 0.2.
#'
#' @noRd
ht_variance.drawn_design_multistage <- function(design, data, rows, y, pi_i,
                                                domain = NULL) {
  cl <- count_clusters(design, data)
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
  NextMethod()
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
                                              domain = NULL) {
  if (design$method == "poisson") {
    return(list(variance = sum((1 - pi_i) / pi_i^2 * y^2), note = NULL))
  }
  if (design$method == "systematic") {
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
#' @noRd
ht_variance.drawn_design_systematic <- function(design, data, rows, y, pi_i,
                                                domain = NULL) {
  if (!is.null(design$start)) {
    stop_declined(paste0("A systematic design with a fixed `start` is not a ",
                         "probability sample: each row is selected with ",
                         "probability 0 or 1."))
  }
  n <- length(rows)
  if (n < 2L) {
    return(list(variance = NA_real_,
                note = "A variance needs at least two sampled rows."))
  }
  walk <- systematic_positions(design, data)
  pos <- walk[rows]
  u <- (y / pi_i)[order(pos)]
  k <- design$interval
  f <- 1 / k
  v <- (1 - f) * n / (2 * (n - 1)) * sum(diff(u)^2)
  N <- sum(!is.na(walk))
  r <- N %% k
  if (r > 0) v <- v + mean(u)^2 * (r / k) * (1 - r / k)
  list(variance = v,
       method = "successive difference",
       note = paste0("Systematic sampling has no design-unbiased variance ",
                     "estimator; this is the successive-difference ",
                     "approximation. It would be understated if the frame ",
                     "cycled with a period matching `interval`, or if a smooth ",
                     "trend along the sort order swamped all other variation.",
                     if (r > 0) {
                       paste0(" The ", N, " rows walked are not a multiple of ",
                              "`interval` (", k, "), so the sample size is ",
                              "random, ", N %/% k, " or ", N %/% k + 1L,
                              " rows, and the variance includes a term for ",
                              "that.")
                     }))
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
#' @noRd
jackknife_variance <- function(design, units, pop, y, pi_i, domain = NULL) {
  refusal <- jackknife_refusal(design)
  if (!is.null(refusal)) {
    return(list(variance = NA_real_, method = "jackknife", note = refusal))
  }
  if (inherits(design, "drawn_design_certainty")) {
    return(jackknife_certainty(design, units, pop, y, pi_i, domain))
  }

  psu <- jackknife_groups(design, units)
  fixed <- !(psu %in% psu[!is_certain(pi_i)])
  if (all(fixed)) return(all_certain("jackknife"))
  n_fixed <- length(unique(psu[fixed]))

  st <- jackknife_strata(design, units, pop, psu, fixed, domain)
  base <- y / pi_i
  total <- sum(base)
  free <- which(!fixed)
  by_stratum <- split(free, factor(match(st$stratum[free], st$levels),
                                   levels = seq_along(st$levels)))

  v <- 0
  n_psu <- 0L
  lonely <- 0L
  for (k in seq_along(st$levels)) {
    in_h <- by_stratum[[k]]
    # Each unit's total in the stratum. The replicates are symmetric in the
    # units, so their order does not matter.
    sums <- rowsum(base[in_h], psu[in_h], reorder = FALSE)[, 1L]
    m <- length(sums)
    f <- st$fpc[k]
    n_psu <- n_psu + m
    if (m < 2L) {
      if (f > 0 && st$in_domain[k]) lonely <- lonely + 1L
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
  if (lonely > 0L) {
    noun <- if (inherits(design, "drawn_design_temporal")) {
      c("interval", "intervals")
    } else {
      c("stratum", "strata")
    }
    return(list(variance = NA_real_, method = "jackknife", note = paste0(
      "The stratified jackknife needs at least two units in each ", noun[1],
      ", and ", lonely, " ", noun[if (lonely == 1L) 1L else 2L],
      if (lonely == 1L) " has" else " have", " one.")))
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
jackknife_certainty <- function(design, units, pop, y, pi_i, domain = NULL) {
  sp <- certainty_split(design, pop)
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
  jackknife_variance(design$rest, units[below, , drop = FALSE],
                     sp$data[sp$rest, , drop = FALSE],
                     y[below], pi_i[below], rest_domain(domain, sp))
}

#' Designs the jackknife should decline rather than approximate
#'
#' @noRd
jackknife_refusal <- function(design) {
  if (inherits(design, "drawn_design_systematic")) {
    return(paste0("A systematic sample has one primary sampling unit -- the ",
                  "random start -- so deleting rows does not reflect the ",
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
#' Stratified and temporal designs delete within strata, with each stratum's
#' own `1 - n_h / N_h`. Everything else is a single stratum, corrected by
#' `jackknife_fpc()`. Only strata holding a unit that is not `fixed` are
#' listed, and both counts leave the fixed units out. A unit taken with
#' certainty is in every sample, so the fixed units in the frame are the ones
#' in the sample. A stratum taken whole drops out.
#'
#' Returns each row's `stratum`, the `levels` listed, their `fpc` and
#' `in_domain`, in the same order. `in_domain` is `FALSE` for a stratum that
#' `domain`, a logical vector over the frame, shows to hold no unit of the
#' domain, and `TRUE` otherwise.
#'
#' @noRd
jackknife_strata <- function(design, units, pop, psu, fixed, domain = NULL) {
  grp <- switch(design_type(design),
    stratified = {
      key <- group_key2(units, pop, design$strata)
      list(sample = as.integer(key$x), pop = as.integer(key$table))
    },
    temporal = list(
      sample = temporal_bucket(design, units),
      pop = temporal_bucket(design, pop)
    ),
    NULL
  )
  if (is.null(grp)) {
    m <- length(unique(psu[!fixed]))
    n_fixed <- length(unique(psu[fixed]))
    return(list(stratum = rep("all", nrow(units)), levels = "all",
                fpc = jackknife_fpc(design, pop, m, n_fixed),
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
  list(stratum = grp$sample, levels = lev,
       fpc = ifelse(free_h > n_h, 1 - n_h / free_h, 0),
       in_domain = in_domain)
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
jackknife_fpc <- function(design, pop, m, n_fixed = 0L) {
  if (inherits(design, "drawn_design_multistage")) return(1)
  total <- if (inherits(design, "drawn_design_cluster")) {
    count_clusters(design, pop)$total
  } else {
    reachable_rows(design, pop)
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
    weighted = keyed(design$weights),
    spread = keyed(c(design$across, design$size)),
    stratified = keyed(design$strata),
    cluster = ,
    multistage = keyed(design$clusters),
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
#' codes; otherwise each row is its own unit.
#'
#' @noRd
jackknife_groups <- function(design, units) {
  if (inherits(design, c("drawn_design_cluster", "drawn_design_multistage"))) {
    return(cluster_codes(units[[design$clusters]]))
  }
  seq_len(nrow(units))
}

#' @export
print.drawn_ht <- function(x, ...) {
  print_estimate(x, "total", "Horvitz-Thompson total")
}
