#' Estimate a population mean
#'
#' The Hajek estimator, `sum(y / pi) / sum(1 / pi)`, is the default. It divides
#' by the *estimated* population size rather than the known one, so a sample
#' that happens to over-represent heavy-weight rows inflates numerator and
#' denominator together and they partly cancel. `estimator = "ht"` divides by
#' the true `N` instead.
#'
#' @section Proportions:
#' For a 0/1 or logical variable either one estimates a proportion, and the
#' confidence interval is then formed on the logit scale and transformed back
#' (Korn and Graubard 1999, section 5.3). A symmetric interval around a small
#' proportion runs below zero; this one cannot, and matches
#' `survey::svyciprop(method = "xlogit")`.
#'
#' Whether `y` is a proportion is read from what `y` is, not from the values
#' the sample happened to draw. A logical `y` is one, and so is a column of
#' the frame whose values are all 0 or 1. A count that came out 0 or 1 in
#' every sampled row is not: a logit interval cannot pass 1, so it would never
#' cover a mean of 1.2. A numeric vector, or a column added to the sample after
#' the draw, such as a 0/1 result recorded in fieldwork, gets the symmetric
#' interval unless `proportion = TRUE` says it is a proportion.
#' `proportion = FALSE` asks for the symmetric interval whatever `y` is.
#'
#' `ci_method = "beta"` gives the Korn and Graubard (1998) interval instead,
#' a Clopper-Pearson interval at an effective sample size of `p(1 - p) / v`,
#' with `p` the estimate and `v` its variance, capped at `n` and then scaled
#' by `(qt(a, n - 1) / qt(a, df))^2` for the design's degrees of freedom.
#' Here `n` is the number of sampled rows below certainty behind the estimate
#' and `a` the tail probability of the interval. With no rows taken with
#' certainty it matches `survey::svyciprop(method = "beta")` whenever the cap
#' does not bind. It ignores the finite population correction, so it errs
#' wide when the sample is a large share of the frame.
#'
#' An estimate of exactly 0 or 1, such as no deviations found in 60 sampled
#' invoices, has a variance of 0, and a logit interval around it has no
#' width. It gets the beta interval whatever `ci_method` says, with the
#' effective sample size set to `n`, and `note` says so. For 0 of 60 from a
#' simple random sample, the 95% upper bound (`alternative = "less"`) is
#' 4.87%. When every sampled row behind the estimate was taken with
#' certainty, the estimate is exact and the interval has no width.
#'
#' @section When the two differ, and which to use:
#' They coincide **exactly** whenever the design weights of the rows you drew
#' sum to `N` — which covers every fixed-size equal-probability design and every
#' [design_stratified()] without replacement, because each stratum contributes
#' `n_h * N_h / n_h = N_h`. There is nothing to choose between them there.
#'
#' They differ when that sum is random or uneven:
#'
#' * the sample size itself is random — [design_weighted()] with
#'   `method = "poisson"`, [design_cluster()] over clusters of unequal size, or
#'   [design_systematic()] with an `interval` that does not divide `N`;
#' * the weights vary within a fixed-size sample — probability-proportional-to-
#'   size selection.
#'
#' Hajek is usually the steadier of the two and is the default for that reason
#' (Särndal, Swensson and Wretman 1992, section 5.7). The exception is worth
#' knowing: when `y` is close to proportional to the size measure that drove
#' selection, `y / pi` is nearly constant, the Horvitz-Thompson numerator
#' barely moves, and dividing it by the known `N` beats dividing by an
#' estimate.
#'
#' One caveat on `"ht"`. It is unbiased for the frame mean *provided every row
#' could have been selected*. Rows with inclusion probability 0 — outside a time
#' window, outside a region, dropped under `na_rm = TRUE` — sit inside the `N`
#' it divides by but can never enter the numerator, so the estimate is biased
#' low by exactly their share of the frame. [sample_summary()] reports how many
#' such rows there are.
#' The Hajek mean is unaffected, because it estimates the mean of the part of
#' the frame the design can actually reach.
#'
#' @section Variance and domains:
#' The Hajek mean is a ratio, and its variance is the variance of the total of
#' the linearised residuals `(y - mean) / N_hat` (Deville 1999), computed with
#' whichever estimator suits the design — see [ht_total()], which also explains
#' the confidence interval's degrees of freedom, how `by` estimates domain
#' means from the whole sample rather than a subset of it, and what an empty
#' sample gives.
#'
#' The residuals sum to zero inside a primary sampling unit that holds every
#' sampled row of the domain, so a Hajek mean over one unit, such as a domain
#' with one sampled row or a Poisson sample of one row, gets no standard
#' error. The "Domains" section of [ht_total()] says more.
#'
#' `contrast` gives a weighted sum of domain means, such as the difference
#' between two, with its standard error; see "Contrasts" in [ht_total()].
#' Its interval is symmetric whatever `proportion` and `ci_method` say.
#'
#' When `y` has the same value in every sampled row, the residuals are 0 and
#' the Hajek mean's variance is set to exactly 0, with a note, rather than
#' left to rounding error, which gave standard errors near `1e-16` and design
#' effects that changed from seed to seed. The design effect is then `NA`.
#'
#' @section Missing values:
#' A missing `y` is refused unless `na_rm = TRUE`. Then the mean is the
#' Hajek mean over the rows with `y` present, treated as a domain of the
#' whole sample, with the design's own variance estimator, and equals
#' `survey::svymean(na.rm = TRUE)` on [as_svydesign()]. It assumes the values
#' are missing completely at random, and the printed result says so. `n`
#' counts the rows with `y` present, `n_missing` the others, and the design
#' effect is `NA`, as for a domain. With `by`, each domain is the rows of
#' its `by` values with `y` present. `estimator = "ht"` is refused, since
#' the number of frame rows that would have a `y` is unknown. See "Missing
#' values and nonresponse" in [ht_total()] for when the assumption fails and
#' what to do then.
#'
#' @param sample A data frame returned by [draw()] with `weights = TRUE`. It is
#'   read as [ht_total()] reads it: `.prob` for the probabilities, and the
#'   stored frame for the design's columns. See "Which columns are read" there.
#' @param y The variable to average: a column name, several column names, or
#'   a numeric vector as long as `sample`. See "Several variables" in
#'   [ht_total()].
#' @param estimator `"hajek"` or `"ht"`. See above.
#' @param variance,level,by,df,alternative,contrast As for [ht_total()].
#' @param proportion Whether `y` is a proportion, which decides the shape of
#'   the interval. `NULL` decides from `y`: a logical `y`, or a column of the
#'   frame whose values are all 0 or 1, is a proportion. `TRUE` treats `y` as
#'   one, and needs it to be 0 or 1 in every sampled row; `FALSE` gives the
#'   symmetric interval. See "Proportions".
#' @param ci_method The interval for a proportion: `"logit"`, or `"beta"` for
#'   Korn and Graubard's. A proportion of exactly 0 or 1 gets `"beta"` either
#'   way. See "Proportions".
#' @param na_rm `FALSE` refuses a missing `y`. `TRUE` estimates the mean over
#'   the rows with `y` present, as a domain of the sample. See "Missing
#'   values".
#'
#' @return A list with a `print()` method, holding `mean`, `estimator`,
#'   `variable`, `variance`, `se`, `ci`, `level`, `alternative`, `df`, `n`,
#'   `n_random`, `design`, `deff`, `method`, `note`, `ci_scale` (`"linear"`,
#'   `"logit"` or `"beta"`), `unreachable` and `n_missing`, the sampled rows
#'   with `y` missing that `na_rm = TRUE` left out. See [ht_total()] for what
#'   `variable`, `n_random`, `method`, `df`, `note` and `unreachable` say,
#'   and [deff()] for the design effect. With `by`, a data frame with one row
#'   per domain, and with `contrast` a `drawn_contrast`, as for [ht_total()].
#'   With several column names as `y`, a `drawn_by` table with a row per
#'   variable, as for [ht_total()]. Each has `coef()`, `vcov()`, `confint()`,
#'   `as.data.frame()` and `tidy()` methods; see "Working with the result" in
#'   [ht_total()].
#'
#' @references
#' Korn, E. L. and Graubard, B. I. (1998). Confidence intervals for
#' proportions with small expected number of positive counts estimated from
#' survey data. *Survey Methodology*, 24, 193–201.
#'
#' Korn, E. L. and Graubard, B. I. (1999). *Analysis of Health Surveys*.
#' Wiley.
#'
#' Hájek, J. (1971). Comment on "An essay on the logical foundations of survey
#' sampling, part one" by D. Basu. In V. P. Godambe and D. A. Sprott (eds.),
#' *Foundations of Statistical Inference*, p. 236. Holt, Rinehart and Winston.
#'
#' Särndal, C.-E., Swensson, B. and Wretman, J. (1992). *Model Assisted Survey
#' Sampling*. Springer.
#'
#' Deville, J.-C. (1999). Variance estimation for complex statistics and
#' estimators: linearization and residual techniques. *Survey Methodology*,
#' 25, 193–203.
#'
#' @examples
#' set.seed(1)
#' pop <- data.frame(
#'   id = 1:200,
#'   site = rep(c("a", "b"), times = c(150, 50)),
#'   spend = round(stats::runif(200, 10, 500))
#' )
#' s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
#'
#' ht_mean(s, "spend")
#' mean(pop$spend)   # the truth
#'
#' # A mean for each site, estimated from the whole sample
#' ht_mean(s, "spend", by = "site")
#'
#' # The difference between the two sites' means
#' ht_mean(s, "spend", by = "site", contrast = c(a = 1, b = -1))
#'
#' # A proportion is the mean of a 0/1 variable
#' ht_mean(s, s$spend > 250)
#'
#' # Stratified without replacement: the weights sum to N, so the two agree
#' all.equal(ht_mean(s, "spend", estimator = "ht")$mean, ht_mean(s, "spend")$mean)
#'
#' # Poisson sampling has a random size, so they part company
#' p <- draw(pop, design_weighted("spend", n = 40, method = "poisson"),
#'           seed = 2, weights = TRUE)
#' c(hajek = ht_mean(p, "spend", variance = "none")$mean,
#'   ht    = ht_mean(p, "spend", "ht", variance = "none")$mean)
#'
#' @seealso [ht_total()], [ht_ratio()], [deff()], [sample_summary()]
#' @export
ht_mean <- function(sample, y, estimator = c("hajek", "ht"),
                    variance = c("auto", "analytic", "jackknife", "none"),
                    level = 0.95, by = NULL, df = NULL, proportion = NULL,
                    ci_method = c("logit", "beta"),
                    alternative = c("two.sided", "less", "greater"),
                    contrast = NULL, na_rm = FALSE) {
  estimator <- check_choice(estimator, c("hajek", "ht"), "estimator")
  variance <- check_choice(variance,
                           c("auto", "analytic", "jackknife", "none"),
                           "variance")
  ci_method <- check_choice(ci_method, c("logit", "beta"), "ci_method")
  alternative <- check_choice(alternative, c("two.sided", "less", "greater"),
                              "alternative")
  check_contrast_by(contrast, by)
  check_flag(na_rm, "na_rm")
  if (several_names(y)) {
    return(estimate_several(sample, y, "mean", variance, level, by, df,
                            alternative, contrast, estimator = estimator,
                            proportion = proportion, ci_method = ci_method,
                            na_rm = na_rm))
  }
  parts <- ht_prepare(sample, y, level, df, na_rm = na_rm)
  parts <- mean_setup(parts, estimator, proportion, ci_method)
  if (!is.null(by)) {
    return(estimate_by(parts, by, variance, level, df, what = "mean",
                       estimator = estimator, ci_method = ci_method,
                       alternative = alternative, contrast = contrast))
  }
  warn_grouped(sample, "ht_mean")
  mean_one(parts, estimator, variance, level, df, ci_method,
           alternative)$result
}

#' Checks a mean makes of its `y`, and whether it is a proportion
#'
#' `label` names `y` in the messages, as the column it is when `y` names
#' several.
#'
#' @noRd
mean_setup <- function(parts, estimator, proportion, ci_method,
                       label = "`y`") {
  n_missing <- sum(!parts$respond)
  if (n_missing && estimator == "ht") {
    stop("`estimator = \"ht\"` divides by the number of frame rows, and ",
         "with ", label, " missing in ", n_missing, " sampled row",
         if (n_missing > 1L) "s", " the number of frame rows the mean is ",
         "over is unknown. Use the Hajek estimator, which estimates it from ",
         "the sample.", call. = FALSE)
  }
  parts$proportion <- mean_is_proportion(parts, proportion)
  if (ci_method == "beta" && !parts$proportion) {
    stop("`ci_method = \"beta\"` gives an interval for a proportion, and ",
         "the mean of ", label, " is not one. Pass `proportion = TRUE` if ",
         "it is a 0/1 indicator; see \"Proportions\" in ?ht_mean.",
         call. = FALSE)
  }
  parts
}

#' One mean, as ht_mean() gives it without `by`
#'
#' Returns the estimate as `result`, with `z`, the variable whose total has
#' its variance, and `domain`, the `domain` its variance was run with, which
#' estimate_several() needs for the covariances between variables.
#'
#' @noRd
mean_one <- function(parts, estimator, variance, level, df, ci_method,
                     alternative, ctx = NULL) {
  if (any(!parts$respond)) {
    return(respondent_mean(parts, variance, level, df, ci_method,
                           alternative, ctx))
  }
  if (!nrow(parts$sample) && estimator == "hajek") {
    stop("`sample` has no rows, so the Hajek mean is 0 / 0. Use ",
         "estimator = \"ht\", which divides the total, 0, by the size of the ",
         "frame.", call. = FALSE)
  }

  pi_i <- parts$pi
  yv <- parts$y
  n_hat <- sum(1 / pi_i)
  N <- nrow(parts$pop)

  denom <- if (estimator == "hajek") n_hat else N
  est <- sum(yv / pi_i) / denom

  # Linearise: the variance of a ratio is the variance of a total of residuals.
  # For the HT mean there is no denominator to vary, so the residual is just y.
  z <- if (estimator == "hajek") (yv - est) / denom else yv / N
  # A constant `y` has residuals of exactly 0. Computed, they come out a
  # rounding error away from it, and so did the standard error.
  constant <- estimator == "hajek" && length(yv) > 0L && all(yv == yv[1])
  if (constant) {
    est <- yv[1]
    z <- numeric(length(yv))
  }
  declined <- if (estimator == "hajek" && variance != "none") {
    ratio_decline(parts, rep(TRUE, length(yv)))
  }

  out <- estimate_one(parts, est, z, variance, level, df, what = "mean",
                      class = "drawn_mean", deff_y = yv,
                      extra = list(estimator = estimator),
                      proportion = parts$proportion, ci_method = ci_method,
                      alternative = alternative, var_out = declined,
                      divide_by = denom_rows(parts, estimator), ctx = ctx)
  if (constant && identical(out$variance, 0) &&
      !identical(out$ci_scale, "beta") && any(!is_certain(pi_i))) {
    out$note <- paste(c(out$note, constant_note(yv[1], parts$proportion)),
                      collapse = " ")
  }
  out$n_missing <- 0L
  list(result = out, z = z, domain = NULL)
}

#' The Hajek mean over the rows with `y` present
#'
#' Those rows are a domain of the whole sample, as a `by` domain is. The
#' estimate is the Hajek mean over them, and its variance that of the total
#' of the linearised residuals `(y - mean) / N_hat` over every sampled row,
#' 0 for a row with `y` missing, from the design's own estimator. Which frame
#' rows would have a `y` cannot be read from the frame, so `domain` is `NA`.
#' This is `survey::svymean(na.rm = TRUE)`. `n`, `n_random` and the beta
#' interval count the rows with `y` present. A domain's design effect is not
#' estimated, so `deff` is `NA`. Returns what mean_one() returns.
#'
#' @noRd
respondent_mean <- function(parts, variance, level, df, ci_method,
                            alternative, ctx = NULL) {
  r <- parts$respond
  pi_i <- parts$pi
  yv <- parts$y
  n_hat <- sum(1 / pi_i[r])
  est <- sum(yv[r] / pi_i[r]) / n_hat
  z <- numeric(length(yv))
  z[r] <- (yv[r] - est) / n_hat
  constant <- all(yv[r] == yv[r][1])
  if (constant) {
    est <- yv[r][1]
    z[] <- 0
  }
  var_out <- if (variance != "none") {
    ratio_decline(parts, which(r), NA,
                  subject = "The domain of rows with `y` present")
  }
  var_out <- var_out %||%
    ht_variance_dispatch(parts$design, parts$units, parts$pop, parts$rows, z,
                         pi_i, variance, domain = NA, ctx = ctx)
  df <- df %||% design_df(parts$design, parts$units, parts$pop, parts$rows,
                          pi_i)
  # The interval's sample size and `n_random` count the rows behind the
  # estimate.
  held <- parts
  held$pi <- pi_i[r]
  out <- finish_estimate(est, var_out, held, level, df, what = "mean",
                         class = "drawn_mean", proportion = parts$proportion,
                         ci_method = ci_method, alternative = alternative)
  out$estimator <- "hajek"
  out$n <- sum(r)
  out$n_missing <- sum(!r)
  notes <- c(out$note, missing_note(sum(!r), length(r)))
  if (constant && identical(out$variance, 0) &&
      !identical(out$ci_scale, "beta") && any(!is_certain(pi_i[r]))) {
    notes <- c(notes, constant_note(est, parts$proportion))
  }
  out$note <- paste(notes, collapse = " ")
  list(result = out, z = z, domain = NA)
}

#' What a mean over the rows with `y` present assumes
#' @noRd
missing_note <- function(n_missing, n) {
  one <- n_missing == 1L
  paste0(n_missing, " of the ", n, " sampled rows ", if (one) "has" else
           "have", " `y` missing and ", if (one) "is" else "are", " left ",
         "out. The mean is over the rows with `y` present, estimated as a ",
         "domain of the whole sample, which assumes the values are missing ",
         "completely at random, with no relation between whether `y` is ",
         "missing and `y` itself. When response differs between groups ",
         "whose `y` differs, adjust the weights within those groups and ",
         "estimate with as_svydesign().")
}

#' What the mean divides its Horvitz-Thompson total by, for the design effect
#'
#' The Hajek mean estimates the mean over the rows the design can reach. The
#' Horvitz-Thompson mean divides by every row of the frame.
#'
#' @noRd
denom_rows <- function(parts, estimator) {
  if (estimator == "ht") nrow(parts$pop) else parts$reach
}

#' Why a Hajek mean of a constant `y` has a variance of 0
#' @noRd
constant_note <- function(value, proportion, domains = NULL) {
  hint <- if (!proportion && value %in% c(0, 1)) {
    paste0(" If `y` is a 0/1 indicator, `proportion = TRUE` gives a ",
           "proportion's interval, which allows for rows with the other value.")
  }
  if (is.null(domains)) {
    return(paste0(
      "`y` is ", format(value), " in every sampled row, so the Hajek mean's ",
      "residuals are 0 and its estimated variance is 0. The interval has no ",
      "width, though rows the sample did not draw may differ.", hint))
  }
  paste0("`y` takes one value in every sampled row of the domain",
         if (length(domains) > 1L) "s", " ",
         paste0("`", domains, "`", collapse = ", "), ", so ",
         if (length(domains) > 1L) "their Hajek means have" else
           "its Hajek mean has",
         " an estimated variance of 0 and an interval with no width, though ",
         "rows the sample did not draw may differ.")
}

#' Whether a mean is a proportion, for the shape of its interval
#'
#' Read from what `y` is rather than from the values the sample happened to
#' draw: a logical `y`, or a column of the frame whose values are all 0 or 1
#' (a logical column of the frame counts too). A count that came out 0 or 1
#' in every sampled row is not one, since a logit interval cannot pass 1 and
#' would never cover a mean above it. A numeric vector, and a column added to
#' the sample after the draw, are proportions only when `proportion` says so.
#'
#' @noRd
mean_is_proportion <- function(parts, proportion) {
  y <- parts$y[parts$respond]
  binary <- all(y == 0 | y == 1)
  if (!is.null(proportion)) {
    if (!is.logical(proportion) || length(proportion) != 1L ||
        is.na(proportion)) {
      stop("`proportion` must be NULL, TRUE or FALSE, not ",
           format_bad(proportion), ".", call. = FALSE)
    }
    if (proportion && !binary) {
      bad <- unique(y[!(y == 0 | y == 1)])
      stop("`proportion = TRUE` needs `y` to be 0 or 1 in every sampled row, ",
           "and it takes ", length(bad), " other value",
           if (length(bad) > 1L) "s", ", such as ", format(bad[1]), ". Leave ",
           "`proportion` at NULL for the mean of a count or a measurement.",
           call. = FALSE)
    }
    return(proportion)
  }
  if (!binary) return(FALSE)
  if (parts$y_logical) return(TRUE)
  col <- parts$y_name
  if (is.null(col)) return(FALSE)
  # compact_sample() records which of the columns it dropped were 0/1.
  if (!col %in% names(parts$pop)) return(col %in% parts$compact$binary)
  binary_column(parts$pop[[col]])
}

#' Estimate a ratio of two population totals
#'
#' `ht_ratio()` estimates `R = Y / X`, the ratio of the population totals of
#' `y` and `x`, by `sum(y / pi) / sum(x / pi)`: audited value per unit of
#' book value, spend per visit, or yield per hectare. A Hajek mean is the
#' ratio with `x` equal to 1 in every row. [ht_total()] with `ratio`
#' multiplies the ratio by the frame's total of `x` to estimate the total of
#' `y`.
#'
#' @section Variance:
#' The ratio is linearised. Its variance is that of the Horvitz-Thompson
#' total of the residuals `(y - R x) / X_hat`, with `X_hat = sum(x / pi)`
#' (Särndal, Swensson and Wretman 1992; Deville 1999), and is computed with
#' the design's own estimator, as for a total: see "Variance" and
#' "Confidence intervals" in [ht_total()] for the estimators, the fallback
#' to the jackknife, the degrees of freedom and the notes. For the designs
#' [as_svydesign()] hands to `survey` with identical standard errors, the
#' standard error equals that of `survey::svyratio()` on its result. For
#' systematic, spread and systematic PPS samples it comes from the
#' estimator this package uses for their totals, which `survey` does not
#' have.
#'
#' The residuals sum to zero over the sample, weighted by `1 / pi`, so a
#' ratio whose sampled rows lie in one primary sampling unit, one row or one
#' cluster, gets no standard error, as a Hajek mean does. When `y` is
#' proportional to `x` in every sampled row the residuals are 0, the
#' variance is set to exactly 0 rather than left to rounding error, and
#' `note` says so.
#'
#' The ratio of two estimated totals is biased, by an amount of order
#' `1 / n`, which is small against its standard error except in small
#' samples. The design effect compares its variance with that of the same
#' ratio from a simple random sample of the same size, as
#' `survey::svyratio(deff = TRUE)` does.
#'
#' @section Domains:
#' `by` gives a ratio for each domain, `R_d = Y_d / X_d` over the rows of
#' the domain, estimated from the whole sample. Each domain's residuals,
#' `(y - R_d x) / X_hat_d`, are 0 outside it, as for the domain totals and
#' means described under "Domains" in [ht_total()], and the standard errors
#' are those of
#' `survey::svyby(~y, ~d, design, survey::svyratio, denominator = ~x)`. A
#' domain whose estimated total of `x` is 0 is refused by name.
#'
#' A missing `y` or `x` is refused. `survey::svyratio(na.rm = TRUE)` on
#' [as_svydesign()] estimates the ratio over the rows that have both.
#'
#' @inheritParams ht_total
#' @param y,x The numerator and the denominator: each a column name of
#'   `sample`, or a numeric vector as long as `sample`.
#'
#' @return A list of class `drawn_ratio` with a `print()` method, holding
#'   `ratio`, the estimate; `variable` and `x`, the names of `y` and `x`
#'   when they are column names, else `NULL`; `total_y` and `total_x`,
#'   their estimated totals; and `variance`, `se`, `ci`, `level`,
#'   `alternative`, `df`, `n`, `n_random`, `ci_scale`, `design`, `deff`,
#'   `method`, `note` and `unreachable` as for [ht_total()]. With `by`, a
#'   `drawn_by` table with a row per domain and a `ratio` column, as for
#'   [ht_total()]. Both have `coef()`, `vcov()`, `confint()`,
#'   `as.data.frame()` and `tidy()` methods; see "Working with the result"
#'   in [ht_total()].
#'
#' @references
#' Deville, J.-C. (1999). Variance estimation for complex statistics and
#' estimators: linearization and residual techniques. *Survey Methodology*,
#' 25, 193–203.
#'
#' Särndal, C.-E., Swensson, B. and Wretman, J. (1992). *Model Assisted Survey
#' Sampling*. Springer.
#'
#' @examples
#' set.seed(1)
#' pop <- data.frame(
#'   id = 1:300,
#'   site = rep(c("a", "b", "c"), each = 100),
#'   book = round(stats::rlnorm(300, 6, 1))
#' )
#' pop$audited <- round(pop$book * stats::runif(300, 0.9, 1.05))
#' s <- draw(pop, design_stratified("site", n = 45), seed = 1, weights = TRUE)
#'
#' # Audited value per unit of book value
#' ht_ratio(s, "audited", "book")
#' sum(pop$audited) / sum(pop$book)   # the truth
#'
#' # For each site, from the whole sample
#' ht_ratio(s, "audited", "book", by = "site")
#'
#' # The ratio estimator of the total uses the frame's total of `book`, and
#' # is far more precise than the Horvitz-Thompson total here
#' ht_total(s, "audited", ratio = "book")
#' ht_total(s, "audited")
#' sum(pop$audited)
#'
#' @seealso [ht_total()], [ht_mean()], [as_svydesign()]
#' @export
ht_ratio <- function(sample, y, x, by = NULL,
                     variance = c("auto", "analytic", "jackknife", "none"),
                     level = 0.95, df = NULL,
                     alternative = c("two.sided", "less", "greater")) {
  variance <- check_choice(variance,
                           c("auto", "analytic", "jackknife", "none"),
                           "variance")
  alternative <- check_choice(alternative, c("two.sided", "less", "greater"),
                              "alternative")
  parts <- prepare_sample(sample, level, df)
  parts <- add_ratio_variables(parts, y, x)
  parts <- finish_prepare(parts)
  if (!is.null(by)) {
    return(estimate_by(parts, by, variance, level, df, what = "ratio",
                       alternative = alternative))
  }
  warn_grouped(sample, "ht_ratio")
  ratio_one(parts, variance, level, df, alternative)$result
}

#' Read the `y` and `x` of a ratio into the parts of a sample
#'
#' Both must be present and finite in every sampled row.
#'
#' @noRd
add_ratio_variables <- function(parts, y, x) {
  sample <- parts$sample
  read <- function(v, arg) {
    out <- read_variable(sample, v, arg)
    gap <- is.na(out$values) & !is.nan(out$values)
    if (any(gap)) {
      stop("`", arg, "` is missing in ", sum(gap), " of the ", length(gap),
           " sampled rows. A ratio needs `y` and `x` in every sampled row. ",
           "survey::svyratio(na.rm = TRUE) on as_svydesign() estimates it ",
           "over the rows that have both.", call. = FALSE)
    }
    need <- "A ratio needs a finite `y` and `x` for every sampled row."
    check_finite_values(out$values, gap, paste0("`", arg, "`"), need)
    out
  }
  yv <- read(y, "y")
  xv <- read(x, "x")
  parts$y <- yv$values
  parts$respond <- rep(TRUE, length(yv$values))
  parts$y_name <- yv$name
  parts$y_logical <- yv$logical
  parts$x <- xv$values
  parts$x_name <- xv$name
  parts
}

#' The ratio of `y` to `x`, and the residuals that carry its variance
#'
#' Returns `ratio`, `sum(y / pi) / sum(x / pi)`; `total_y` and `total_x`;
#' `z`, the residuals `(y - ratio * x) / total_x`; and `proportional`,
#' whether `y` is `ratio * x` in every row up to rounding. Then `z` is set
#' to exactly 0: the subtraction leaves a rounding error, which would give a
#' standard error a little above 0 and a design effect made of that error.
#' The caller refuses a `total_x` of 0.
#'
#' @noRd
ratio_fit <- function(y, x, pi_i) {
  total_y <- sum(y / pi_i)
  total_x <- sum(x / pi_i)
  ratio <- total_y / total_x
  e <- y - ratio * x
  proportional <- length(y) > 0L && is.finite(ratio) &&
    all(abs(e) <= 1e-9 * pmax(abs(y), abs(ratio * x)))
  z <- if (proportional) numeric(length(y)) else e / total_x
  list(ratio = ratio, total_y = total_y, total_x = total_x, z = z,
       proportional = proportional)
}

#' A variable as a message names it: its column name, or `arg` for a vector
#' @noRd
ratio_name <- function(name, arg) paste0("`", name %||% arg, "`")

#' Refuse a ratio whose estimated denominator is 0
#' @noRd
check_ratio_denominator <- function(parts, total_x, name) {
  if (!nrow(parts$sample)) {
    stop("`sample` has no rows, so the ratio is 0 / 0.", call. = FALSE)
  }
  if (total_x == 0) {
    stop("The estimated total of ", name, ", `sum(x / pi)` over the ",
         "sampled rows, is 0, so the ratio is undefined. Use a denominator ",
         "whose sampled values do not cancel out.", call. = FALSE)
  }
  invisible(NULL)
}

#' One ratio, as ht_ratio() gives it without `by`
#'
#' Returns what mean_one() returns.
#'
#' @noRd
ratio_one <- function(parts, variance, level, df, alternative, ctx = NULL) {
  fit <- ratio_fit(parts$y, parts$x, parts$pi)
  check_ratio_denominator(parts, fit$total_x, ratio_name(parts$x_name, "x"))
  declined <- if (variance != "none") {
    ratio_decline(parts, rep(TRUE, length(parts$y)), estimate = "ratio")
  }
  out <- estimate_one(parts, fit$ratio, fit$z, variance, level, df,
                      what = "ratio", class = "drawn_ratio", deff_y = fit$z,
                      alternative = alternative, var_out = declined,
                      ctx = ctx)
  out$x <- parts$x_name
  out$total_y <- fit$total_y
  out$total_x <- fit$total_x
  out <- proportional_result(out, fit, parts)
  list(result = out, z = fit$z, domain = NULL)
}

#' The ratio estimator of a total, `R * X`, with `X` the frame's total of `x`
#'
#' Its variance is `X^2` times the ratio's, the variance of the total of
#' the residuals `X (y - R x) / X_hat`.
#'
#' @noRd
ratio_total <- function(parts, variance, level, df, alternative, ctx = NULL) {
  aux <- parts$ratio
  fit <- ratio_fit(parts$y, aux$x, parts$pi)
  check_ratio_denominator(parts, fit$total_x, ratio_name(aux$name, "x"))
  z <- aux$total * fit$z
  declined <- if (variance != "none") {
    ratio_decline(parts, rep(TRUE, length(parts$y)),
                  estimate = "ratio estimator")
  }
  out <- estimate_one(parts, fit$ratio * aux$total, z, variance, level, df,
                      what = "total", class = "drawn_ht", deff_y = z,
                      alternative = alternative, var_out = declined,
                      ctx = ctx)
  out$ratio <- fit$ratio
  out$x <- aux$name
  out$frame_total <- aux$total
  out <- proportional_result(out, fit, parts)
  list(result = out, z = z, domain = NULL)
}

#' Add the note for a `y` proportional to `x`, when its variance is 0
#' @noRd
proportional_result <- function(out, fit, parts) {
  if (fit$proportional && identical(out$variance, 0) &&
      any(!is_certain(parts$pi))) {
    out$note <- paste(c(out$note, proportional_note(fit$ratio)),
                      collapse = " ")
  }
  out
}

#' Why a ratio of a `y` proportional to `x` has a variance of 0
#' @noRd
proportional_note <- function(value = NULL, domains = NULL) {
  tail <- "though rows the sample did not draw may differ."
  if (is.null(domains)) {
    return(paste0(
      "`y / x` is ", format(value), " in every sampled row, so the ",
      "ratio's residuals are 0 and its estimated variance is 0. The ",
      "interval has no width, ", tail))
  }
  several <- length(domains) > 1L
  paste0("`y / x` takes one value in every sampled row of the domain",
         if (several) "s", " ", paste0("`", domains, "`", collapse = ", "),
         ", so ", if (several) "their ratios have" else "its ratio has",
         " an estimated variance of 0 and an interval with no width, ", tail)
}

#' A ratio estimator of a total needs no domains
#' @noRd
check_ratio_by <- function(ratio, by) {
  if (!is.null(ratio) && !is.null(by)) {
    stop("`ratio` gives one ratio estimator of the total, and `by` asks for ",
         "domains. For each domain's ratio use ht_ratio() with `by`, and ",
         "multiply it by the domain's total of `", ratio[1], "` in the ",
         "frame.", call. = FALSE)
  }
  invisible(NULL)
}

#' The column `ratio` names, and its total over the frame
#'
#' The total is over the frame rows the design can reach, the population a
#' total estimates, read from the frame stored with the sample, as
#' frame_totals() reads it. Returns `NULL` without `ratio`, else its
#' `name`, the frame's value `x` for each sampled row and `total`.
#'
#' @noRd
ratio_parts <- function(parts, ratio) {
  if (is.null(ratio)) return(NULL)
  if (!is.character(ratio) || length(ratio) != 1L || is.na(ratio)) {
    stop("`ratio` must be the name of one column of the frame, not ",
         format_bad(ratio), ".", call. = FALSE)
  }
  pop <- parts$pop
  if (!ratio %in% names(pop)) {
    if (ratio %in% parts$compact$dropped) {
      stop("`ratio` names `", ratio, "`, which compact_sample() dropped ",
           "from the frame stored with `sample`. Name it in `keep` when ",
           "compacting, or rebuild the sample from the full frame with ",
           "as_drawn().", call. = FALSE)
    }
    stop("`ratio` names `", ratio, "`, which is not a column of the frame ",
         "`sample` was drawn from. The ratio estimator needs its total over ",
         "the frame, so name a column of that frame.", call. = FALSE)
  }
  x_all <- pop[[ratio]]
  if (!is.numeric(x_all) && !is.logical(x_all)) {
    stop("`ratio` names `", ratio, "`, a ", class(x_all)[1], " column of ",
         "the frame. Name a numeric column.", call. = FALSE)
  }
  x_all <- as.numeric(x_all)
  # A sample saved by drawn 0.2.0 records no probabilities, so they are
  # worked out again from the design.
  pi_all <- parts$pi_all %||% exact_inclusion(parts$design, pop)
  reach <- reachable(pi_all)
  vals <- x_all[reach]
  gap <- is.na(vals)
  if (any(gap)) {
    stop("`", ratio, "` is missing in ", sum(gap), " of the ", length(vals),
         " frame rows the design could reach, so its frame total is ",
         "unknown. Name a column known for every row the design can reach, ",
         "or estimate the ratio with ht_ratio() and multiply it by a total ",
         "you know.", call. = FALSE)
  }
  if (any(!is.finite(vals))) {
    stop("`", ratio, "` is infinite in ", sum(!is.finite(vals)), " of the ",
         "frame rows the design could reach, so its frame total is not a ",
         "number. Name a column with a finite value in every row.",
         call. = FALSE)
  }
  list(name = ratio, x = x_all[parts$rows], total = sum(vals))
}

#' @export
print.drawn_ratio <- function(x, ...) {
  shown <- function(v) format(v, big.mark = ",")
  print_estimate(x, "ratio", "Ratio of Horvitz-Thompson totals", lines = c(
    ratio = paste(ratio_name(x$variable, "y"), "/", ratio_name(x$x, "x")),
    totals = paste(shown(x$total_y), "/", shown(x$total_x))
  ))
}

# ---- several variables at once ---------------------------------------------

#' Whether `y` names several columns
#' @noRd
several_names <- function(y) is.character(y) && length(y) > 1L

#' Whether a `by` table holds the estimates of several variables
#' @noRd
several_y <- function(x) length(attr(x, "variable")) > 1L

#' Check the names in a `y` that names several columns
#'
#' A character vector that names no column of `sample`, or holds a missing
#' or empty name, is taken for data passed as `y`, and gets the message that
#' lists the forms `y` takes.
#'
#' @noRd
check_several_y <- function(sample, y) {
  if (anyNA(y) || !all(nzchar(y)) || !any(y %in% names(sample))) {
    stop(variable_forms(y, "y", TRUE, sample), call. = FALSE)
  }
  twice <- unique(y[duplicated(y)])
  if (length(twice)) {
    stop("`y` names ", quote_names(twice), " more than once. Name each ",
         "column once.", call. = FALSE)
  }
  absent <- setdiff(y, names(sample))
  if (length(absent)) {
    stop("`sample` has no column", if (length(absent) > 1L) "s", " ",
         quote_names(absent), ".", call. = FALSE)
  }
  invisible(NULL)
}

#' Estimates of several variables, as one table
#'
#' Each variable is estimated as a call for it alone estimates it, with one
#' `variance_context()` for them all, and the rows are stacked under a
#' `variable` column, which the table's `by` attribute leads, so that
#' coef(), confint() and the other methods name each row by its variable and
#' domain.
#'
#' The covariances between the variables come from polarisation. Every
#' variance estimator here is a quadratic form in the variable it totals,
#' `V(z) = z' A z` with `A` set by the design and the sampled rows, so
#' `Cov(a, b) = (V(z_a + z_b) - V(z_a) - V(z_b)) / 2`, one more run of the
#' same estimator for each pair. They are kept in the `vcov` attribute,
#' which vcov() reads. With `by` they are the covariances between the
#' variables within each domain, and 0 between domains, as in the vcov() of
#' a `by` table of one variable.
#'
#' @noRd
estimate_several <- function(sample, y, what, variance, level, by, df,
                             alternative, contrast, estimator = "hajek",
                             proportion = NULL, ci_method = "logit",
                             na_rm = NULL, ratio = NULL) {
  check_several_y(sample, y)
  if (!is.null(contrast)) {
    stop("`contrast` weighs the domains of one variable, and `y` names ",
         length(y), " columns. Estimate each contrast with a call of its ",
         "own.", call. = FALSE)
  }
  base <- prepare_sample(sample, level, df)
  parts <- lapply(y, function(v) {
    add_y(base, v, na_rm, label = paste0("`", v, "`, named in `y`,"))
  })
  base <- finish_prepare(base)
  if (what == "mean") {
    parts <- lapply(seq_along(y), function(i) {
      mean_setup(parts[[i]], estimator, proportion, ci_method,
                 label = paste0("`", y[i], "`"))
    })
  }
  aux <- ratio_parts(base, ratio)
  for (i in seq_along(parts)) parts[[i]]["ratio"] <- list(aux)
  ctx <- variance_context()
  if (is.null(by)) {
    warn_grouped(sample, paste0("ht_", what))
    fits <- lapply(parts, function(p) {
      if (what == "total") {
        total_one(p, variance, level, df, alternative, ctx)
      } else {
        mean_one(p, estimator, variance, level, df, ci_method, alternative,
                 ctx)
      }
    })
    results <- lapply(fits, `[[`, "result")
    pieces <- lapply(results, function(r) {
      row <- data.frame(n = as.integer(r$n))
      row[[what]] <- r[[what]]
      row$variance <- r$variance
      row$se <- r$se
      row$ci_lower <- r$ci[1]
      row$ci_upper <- r$ci[2]
      row$method <- r$method
      row
    })
    notes <- lapply(results, `[[`, "note")
    n_random <- lapply(results, `[[`, "n_random")
    df_use <- results[[1]]$df
  } else {
    fits <- lapply(parts, function(p) {
      estimate_by(p, by, variance, level, df, what, estimator = estimator,
                  ci_method = ci_method, alternative = alternative,
                  ctx = ctx, reserved = "variable", keep = TRUE)
    })
    tables <- lapply(fits, `[[`, "table")
    pieces <- lapply(tables, plain_frame)
    notes <- lapply(tables, attr, "notes")
    n_random <- lapply(tables, function(t) unname(attr(t, "n_random")))
    df_use <- attr(tables[[1]], "df")
  }
  sizes <- vapply(pieces, nrow, integer(1))
  out <- do.call(rbind, lapply(seq_along(y), function(i) {
    cbind(data.frame(variable = rep(y[i], sizes[i]),
                     stringsAsFactors = FALSE), pieces[[i]])
  }))
  rownames(out) <- NULL
  out <- structure(out, class = c("drawn_by", "data.frame"), what = what,
                   estimator = if (what == "mean") estimator, level = level,
                   alternative = alternative, df = df_use,
                   design = design_type(base$design),
                   notes = several_notes(y, notes),
                   unreachable = base$unreachable, variable = y,
                   by = c("variable", by))
  if (!is.null(aux)) attr(out, "x") <- aux$name
  nm <- domain_names(out)
  # A proportion's interval is on the scale of `ci_method`, any other
  # variable's on the linear scale; one scale for all is kept as one value.
  scale <- rep(vapply(parts, function(p) {
    if (isTRUE(p$proportion)) ci_method else "linear"
  }, character(1)), sizes)
  attr(out, "ci_scale") <- if (length(unique(scale)) > 1L) {
    stats::setNames(scale, nm)
  } else {
    c(scale, "linear")[1]
  }
  attr(out, "n_random") <- stats::setNames(as.integer(unlist(n_random)), nm)
  attr(out, "vcov") <- if (is.null(by)) {
    several_vcov(base, fits, variance, nm, ctx)
  } else {
    several_vcov_by(base, by, fits, out, variance, nm, ctx)
  }
  out
}

#' The notes of several variables: once if every variable has it, else
#' after the variable's name
#' @noRd
several_notes <- function(y, notes) {
  notes <- lapply(notes, function(n) n[!is.na(n) & nzchar(n)])
  every <- unique(unlist(notes))
  shared <- every[vapply(every, function(n) {
    all(vapply(notes, function(v) n %in% v, logical(1)))
  }, logical(1))]
  own <- unlist(lapply(seq_along(y), function(i) {
    mine <- setdiff(notes[[i]], shared)
    if (length(mine)) paste0("`", y[i], "`: ", mine)
  }))
  c(shared, own)
}

#' A covariance by polarisation, or `NA`
#'
#' `va` and `vb` are the two variances and `ma` and `mb` the estimators
#' that gave them. `run(use)` gives the variance of the sum with `variance =
#' use`. The covariance is `NA` when either variance is missing or the
#' three did not come from one estimator, as when the design's own declined
#' for one of them and the jackknife stood in.
#'
#' @noRd
polarised <- function(va, vb, ma, mb, run) {
  if (!is.finite(va) || !is.finite(vb) || !identical(ma, mb) ||
      identical(ma, "none")) {
    return(NA_real_)
  }
  s <- run(if (ma == "jackknife") "jackknife" else "analytic")
  if (!identical(s$method %||% "analytic", ma) || !is.finite(s$variance)) {
    return(NA_real_)
  }
  (s$variance - va - vb) / 2
}

#' The covariance matrix of several variables' estimates, without `by`
#' @noRd
several_vcov <- function(base, fits, variance, nm, ctx) {
  res <- lapply(fits, `[[`, "result")
  v <- vapply(res, function(r) r$variance, numeric(1))
  m <- vapply(res, function(r) r$method %||% "analytic", character(1))
  k <- length(fits)
  out <- diag(v, k)
  for (a in seq_len(k - 1L)) {
    for (b in seq.int(a + 1L, k)) {
      # A mean over the rows with `y` present was run on a domain the frame
      # cannot place, and the sum is run the same way.
      domain <- if (anyNA(c(fits[[a]]$domain, fits[[b]]$domain))) NA
      z <- fits[[a]]$z + fits[[b]]$z
      out[a, b] <- out[b, a] <- polarised(v[a], v[b], m[a], m[b],
                                          function(use) {
        ht_variance_dispatch(base$design, base$units, base$pop, base$rows, z,
                             base$pi, use, domain, ctx)
      })
    }
  }
  dimnames(out) <- list(nm, nm)
  out
}

#' The covariance matrix of several variables' domain estimates
#'
#' Each pair's sum is run over every domain at once, with the rows of each
#' domain whether or not they have each `y`; a row without one has a `z` of
#' 0 for it.
#'
#' @noRd
several_vcov_by <- function(base, by, fits, out, variance, nm, ctx) {
  v <- out$variance
  V <- diag(v, length(v))
  dimnames(V) <- list(nm, nm)
  if (variance == "none") return(V)
  every <- base
  every$respond <- rep(TRUE, nrow(base$sample))
  dom <- domain_keys(every, by)
  start <- cumsum(c(0L, vapply(fits, function(f) nrow(f$table),
                               integer(1))))
  k <- length(fits)
  for (a in seq_len(k - 1L)) {
    for (b in seq.int(a + 1L, k)) {
      vs <- domain_variances(every, dom, list(z = fits[[a]]$z + fits[[b]]$z),
                             variance, NULL, ctx)
      for (j in seq_along(dom$levels_in)) {
        ra <- start[a] + match(dom$levels_in[j], fits[[a]]$dom$levels_in)
        rb <- start[b] + match(dom$levels_in[j], fits[[b]]$dom$levels_in)
        if (is.na(ra) || is.na(rb)) next
        V[ra, rb] <- V[rb, ra] <- polarised(v[ra], v[rb], out$method[ra],
                                            out$method[rb],
                                            function(use) vs[[j]])
      }
    }
  }
  V
}

#' Estimates for each domain, from the whole sample
#'
#' A domain estimate uses the full sample with `y` zeroed outside the domain.
#' Subsetting first would treat the number of rows that happened to land in the
#' domain as fixed by design, which it is not, and would hand the variance
#' estimator a sample the design never drew.
#'
#' Every domain's variance is computed against one `variance_context()`, so
#' the strata, clusters and other pieces read from the frame are worked out
#' once per call rather than once per domain. With `contrast`, the domains'
#' estimates are combined by `estimate_contrast()` instead.
#'
#' `what` is `"total"`, `"mean"` or `"ratio"`. `ctx` is the context to share
#' with other calls on the same sample, `reserved` names further columns the
#' table uses, which `by` must not name, and `keep = TRUE` returns the table
#' with the domains and their variable `z`, as estimate_several() needs them.
#'
#' @noRd
estimate_by <- function(parts, by, variance, level, df, what,
                        estimator = "hajek", ci_method = "logit",
                        alternative = "two.sided", contrast = NULL,
                        ctx = variance_context(), reserved = NULL,
                        keep = FALSE) {
  dom <- domain_keys(parts, by)
  est <- domain_estimates(parts, dom, what, estimator)
  pi_i <- parts$pi
  df_use <- df %||% design_df(parts$design, parts$units, parts$pop,
                              parts$rows, pi_i)
  if (!is.null(contrast)) {
    return(estimate_contrast(parts, dom, est, contrast, variance, level,
                             df_use, what, estimator, alternative, ctx))
  }
  # The table's own columns would overwrite a `by` column of the same name.
  own <- c(reserved, "n", what, "variance", "se", "ci_lower", "ci_upper",
           "method")
  clash <- intersect(by, own)
  if (length(clash)) {
    stop("`by` names ", paste0("`", clash, "`", collapse = ", "), ", which ",
         "the table of domain estimates uses for its own column",
         if (length(clash) > 1L) "s", ". Rename ",
         if (length(clash) > 1L) "them" else "it", " in `sample` before ",
         "estimating.", call. = FALSE)
  }
  linearised <- if (what == "ratio") {
    "ratio"
  } else if (what == "mean" && estimator == "hajek") {
    "Hajek mean"
  }
  vars <- domain_variances(parts, dom, est, variance, linearised, ctx = ctx)
  proportion <- what == "mean" && isTRUE(parts$proportion)
  labels <- dom$labels
  levels_in <- dom$levels_in

  rows <- lapply(seq_along(levels_in), function(j) {
    v <- vars[[j]]
    ci <- ci_for(est$est[j], v$variance, level, df_use, alternative,
                 proportion, ci_method, pi_i[dom$idx[[j]]])
    list(variance = v$variance,
         se = if (is.finite(v$variance) && v$variance >= 0)
           sqrt(v$variance) else NA_real_,
         lo = ci$ci[1], hi = ci$ci[2], method = v$method %||% "analytic",
         note = v$note, edge = ci$edge,
         constant = est$constant[j] && identical(v$variance, 0) &&
           !ci$edge && !ci$exact)
  })

  first <- match(levels_in, dom$key)
  out <- data.frame(stats::setNames(lapply(by, function(b) {
    parts$sample[[b]][first]
  }), by), check.names = FALSE)
  out$n <- lengths(dom$idx)
  out[[what]] <- est$est
  out$variance <- vapply(rows, `[[`, numeric(1), "variance")
  out$se <- vapply(rows, `[[`, numeric(1), "se")
  out$ci_lower <- vapply(rows, `[[`, numeric(1), "lo")
  out$ci_upper <- vapply(rows, `[[`, numeric(1), "hi")
  out$method <- vapply(rows, `[[`, character(1), "method")
  notes <- unique(unlist(lapply(rows, `[[`, "note")))
  if (any(is.finite(out$se) & (is.na(out$ci_lower) | is.na(out$ci_upper)))) {
    notes <- c(notes, df_note(df_use))
  }
  edge <- vapply(rows, `[[`, logical(1), "edge")
  if (any(edge)) {
    notes <- c(notes, edge_note(domains = labels[levels_in[edge]],
                                ci_method = ci_method,
                                certain = any(is_certain(pi_i))))
  }
  constant <- vapply(rows, `[[`, logical(1), "constant")
  if (any(constant)) {
    shown <- labels[levels_in[constant]]
    notes <- c(notes, if (what == "ratio") {
      proportional_note(domains = shown)
    } else {
      constant_note(NA, TRUE, domains = shown)
    })
  }
  if (any(!parts$respond)) {
    notes <- c(notes, missing_note(sum(!parts$respond),
                                   length(parts$respond)))
  }
  out <- structure(out, class = c("drawn_by", "data.frame"),
                   what = what, estimator = if (what == "mean") estimator,
                   level = level, alternative = alternative, df = df_use,
                   design = design_type(parts$design),
                   ci_scale = if (proportion) ci_method else "linear",
                   notes = notes, unreachable = parts$unreachable,
                   variable = parts$y_name, by = by)
  if (what == "ratio") attr(out, "x") <- parts$x_name
  # The beta interval's sample size, kept by domain name rather than by
  # position so that it still matches after rows are selected.
  attr(out, "n_random") <- stats::setNames(
    vapply(dom$idx, function(i) sum(!is_certain(pi_i[i])), integer(1)),
    domain_names(out))
  if (keep) return(list(table = out, z = est$z, dom = dom))
  out
}

#' The domains of `by`, keyed in the sample and, where it can be, the frame
#'
#' Returns `key`, each sampled row's domain as a code into `labels`;
#' `pop_key`, each frame row's, or `NULL` when the frame lacks a `by` column;
#' `levels_in`, the codes of the domains with a sampled row, in key order;
#' `pos`, each sampled row's position in `levels_in`, `NA` for a row with `y`
#' missing under `na_rm = TRUE`; `idx`, the sampled rows of each domain in
#' `levels_in`, without those; and `aligned`, whether the frame's `by`
#' columns agree with the sample's on the sampled rows, so that the frame can
#' say which of its rows each domain holds.
#'
#' @noRd
domain_keys <- function(parts, by) {
  if (!is.character(by) || !length(by) || anyNA(by)) {
    stop("`by` must name one or more columns of `sample`.", call. = FALSE)
  }
  absent <- setdiff(by, names(parts$sample))
  if (length(absent)) {
    stop("`by` names ", paste0("`", absent, "`", collapse = ", "), ", which ",
         "`sample` does not have.", call. = FALSE)
  }
  check_key_columns(parts$sample, by, "by")
  # Columns are read one at a time, so an sf sample's geometry stays out.
  if (any(vapply(by, function(b) anyNA(parts$sample[[b]]), logical(1)))) {
    stop("`by` columns must not contain missing values.", call. = FALSE)
  }
  # Domains come out in key order: factor levels, else sorted values, with
  # text in byte order whatever the locale.
  if (all(by %in% names(parts$pop))) {
    keys <- group_key2(parts$sample, parts$pop, by)
    key <- as.integer(keys$x)
    pop_key <- as.integer(keys$table)
    labels <- group_label(levels(keys$x))
  } else {
    k <- group_key(parts$sample, by)
    key <- as.integer(k)
    pop_key <- NULL
    labels <- group_label(levels(k))
  }
  # Rows with `y` missing, left out by `na_rm = TRUE`, belong to no domain,
  # so each domain is the rows of its `by` values with `y` present.
  use <- parts$respond %||% rep(TRUE, length(key))
  levels_in <- which(tabulate(key[use], length(labels)) > 0L)
  pos <- match(key, levels_in)
  pos[!use] <- NA_integer_
  idx <- split(seq_along(key)[use],
               structure(pos[use], levels = as.character(seq_along(levels_in)),
                         class = "factor"))
  list(by = by, key = key, pop_key = pop_key, labels = labels,
       levels_in = levels_in, pos = pos, idx = unname(idx),
       aligned = !is.null(pop_key) && identical(pop_key[parts$rows], key))
}

#' Which frame rows the domains at positions `j` of `levels_in` hold
#'
#' Read from the frame's own `by` columns. Where those are missing, or
#' disagree with the sample's values for the sampled rows, the frame cannot
#' say and the result is `NA`.
#'
#' @noRd
domain_frame <- function(dom, j) {
  if (!dom$aligned) return(NA)
  dom$pop_key %in% dom$levels_in[j]
}

#' Each domain's estimate, and the variable whose total has its variance
#'
#' For a total that variable is `y` on the domain's rows; for a Hajek mean,
#' the linearised residual `(y - mean) / N_hat` of the domain; for a
#' Horvitz-Thompson mean, `y / N_d` with `N_d` the domain's rows in the frame;
#' and for a ratio, `(y - R_d x) / X_hat_d` with `R_d` the domain's ratio and
#' `X_hat_d` its estimated total of `x`.
#' Each is 0 outside its domain, so one vector `z` holds them all: each
#' sampled row's value in its own domain. A Hajek mean of a `y` constant over
#' the domain, and a ratio of a `y` proportional to `x` there, have residuals
#' of exactly 0, not the rounding error the subtraction leaves, and
#' `constant` marks them.
#'
#' @noRd
domain_estimates <- function(parts, dom, what, estimator) {
  yv <- parts$y
  pi_i <- parts$pi
  k <- length(dom$levels_in)
  est <- numeric(k)
  constant <- logical(k)
  z <- numeric(length(yv))
  if (what == "ratio") {
    zero <- logical(k)
    for (j in seq_len(k)) {
      i <- dom$idx[[j]]
      fit <- ratio_fit(yv[i], parts$x[i], pi_i[i])
      zero[j] <- fit$total_x == 0
      est[j] <- fit$ratio
      constant[j] <- fit$proportional
      z[i] <- fit$z
    }
    if (any(zero)) {
      shown <- dom$labels[dom$levels_in[zero]]
      several <- length(shown) > 1L
      stop("The estimated total of ", ratio_name(parts$x_name, "x"),
           " is 0 in the domain", if (several) "s", " ",
           quote_labels(shown), ", so ", if (several) "their ratios are" else
             "its ratio is", " undefined. Merge ",
           if (several) "those domains" else "that domain",
           " with another in the `by` columns.", call. = FALSE)
    }
    return(list(est = est, z = z, constant = constant))
  }
  size <- NULL
  if (what == "mean" && estimator == "ht") {
    if (is.null(dom$pop_key)) {
      gone <- intersect(dom$by, parts$compact$dropped)
      if (length(gone)) {
        stop("estimator = \"ht\" needs each domain's size in the frame, ",
             "and compact_sample() dropped ", quote_names(gone), " from the ",
             "frame stored with `sample`. Name the `by` columns in `keep` ",
             "when compacting, or rebuild the sample from the full frame ",
             "with as_drawn().", call. = FALSE)
      }
      stop("estimator = \"ht\" needs the domain's population size, and ",
           "the `by` columns are not in the population.", call. = FALSE)
    }
    size <- tabulate(dom$pop_key, length(dom$labels))[dom$levels_in]
  }
  for (j in seq_len(k)) {
    i <- dom$idx[[j]]
    total <- sum(yv[i] / pi_i[i])
    if (what == "total") {
      est[j] <- total
      z[i] <- yv[i]
      next
    }
    denom <- if (estimator == "hajek") sum(1 / pi_i[i]) else size[j]
    est[j] <- total / denom
    in_d <- yv[i]
    if (estimator == "hajek" && all(in_d == in_d[1])) {
      constant[j] <- TRUE
      est[j] <- in_d[1]
    } else {
      z[i] <- if (estimator == "hajek") (in_d - est[j]) / denom else
        in_d / denom
    }
  }
  list(est = est, z = z, constant = constant)
}

#' The variance of every domain's estimate
#'
#' A Hajek mean or a ratio over one primary sampling unit is declined first,
#' with `ratio_decline()`; `linearised` names the estimate in its note, and
#' is `NULL` for an estimate that has no such decline. Where the design's
#' analytic estimator is
#' `srs_variance()` over strata, intervals or one group, every other domain is
#' estimated in one pass with `srs_variance_by()`. A domain that pass cannot
#' settle, because a lonely group holds some of its units or its variance
#' overflowed, goes through `ht_variance_dispatch()` as every domain of the
#' other designs does, with the shared `ctx`, and gets the same note and
#' fallback as an estimate of its own.
#'
#' @noRd
domain_variances <- function(parts, dom, est, variance, linearised, ctx) {
  k <- length(dom$levels_in)
  out <- vector("list", k)
  if (!is.null(linearised) && variance != "none") {
    psu <- sample_units(parts$design, parts$units, parts$pop, parts$rows)
    for (j in seq_len(k)) {
      out[j] <- list(ratio_decline(parts, dom$idx[[j]], domain_frame(dom, j),
                                   dom$labels[dom$levels_in[j]], psu,
                                   estimate = linearised))
    }
  }
  open <- which(vapply(out, is.null, logical(1)))
  if (length(open) && variance %in% c("auto", "analytic")) {
    v <- srs_domains(parts, dom, est$z, ctx)
    for (j in open[!is.na(v[open])]) {
      out[[j]] <- list(variance = v[[j]], note = NULL, method = "analytic")
    }
    open <- which(vapply(out, is.null, logical(1)))
  }
  for (j in open) {
    i <- dom$idx[[j]]
    z <- numeric(length(est$z))
    z[i] <- est$z[i]
    out[[j]] <- ht_variance_dispatch(parts$design, parts$units, parts$pop,
                                     parts$rows, z, parts$pi, variance,
                                     domain_frame(dom, j), ctx)
  }
  out
}

#' Every domain's `srs_variance()` in one pass, or `NA` where it cannot say
#'
#' One variance per domain. Every domain gets `NA` when the design's analytic
#' estimator is not `srs_variance()` over groups (see `srs_layout()`), or when
#' a lonely group exists and the frame cannot say which domains it holds. A
#' domain with units in a lonely group, and a variance that overflowed, get
#' `NA` too.
#'
#' @noRd
srs_domains <- function(parts, dom, z, ctx) {
  lay <- srs_layout(parts$design, parts$pop, parts$rows, ctx)
  k <- length(dom$levels_in)
  if (is.null(lay)) return(rep(NA_real_, k))
  keep <- !is.na(lay$group)
  # A row with `y` missing is in no domain. It still counts among its
  # group's sampled rows, so it goes in a domain of its own, left out of the
  # result, where its `z` of 0 adds nothing.
  pos <- dom$pos
  pos[is.na(pos)] <- k + 1L
  v <- srs_variance_by(z[keep], lay$group[keep], pos[keep], k + 1L, lay$n_g,
                       lay$N_g)[seq_len(k)]
  if (length(lay$lonely_rows)) {
    hit <- if (dom$aligned) {
      match(dom$pop_key[lay$lonely_rows], dom$levels_in)
    } else {
      seq_len(k)
    }
    v[hit[!is.na(hit)]] <- NA_real_
  }
  v[!is.finite(v)] <- NA_real_
  v
}

#' A weighted sum of domain estimates, with its standard error
#'
#' Each domain's estimate comes from the whole sample, with the variable
#' `domain_estimates()` gives it, so a weighted sum of the estimates is
#' estimated by the same sum of those variables, and its variance is that of
#' one total. One run of the design's estimator gives it, covariances between
#' domains included. The interval is symmetric, since a difference of
#' proportions can be negative. A Hajek mean over one primary sampling unit
#' has no measurable variance, so a contrast that weighs one has none either.
#'
#' @noRd
estimate_contrast <- function(parts, dom, est, contrast, variance, level, df,
                              what, estimator, alternative, ctx) {
  w <- contrast_weights(contrast, dom)
  used <- which(w != 0)
  value <- sum(w[used] * est$est[used])
  w_row <- w[dom$pos]
  w_row[is.na(w_row)] <- 0
  z <- w_row * est$z
  labels <- dom$labels[dom$levels_in[used]]
  var_out <- NULL
  if (what == "mean" && estimator == "hajek" && variance != "none") {
    psu <- sample_units(parts$design, parts$units, parts$pop, parts$rows)
    notes <- unlist(lapply(seq_along(used), function(u) {
      ratio_decline(parts, dom$idx[[used[u]]], domain_frame(dom, used[u]),
                    labels[u], psu)$note
    }))
    if (length(notes)) {
      var_out <- list(variance = NA_real_, method = "none", note = paste(
        c(notes, paste0("The contrast weighs ",
                        if (length(notes) > 1L) "these domains" else
                          "this domain",
                        ", so it has no standard error either.")),
        collapse = " "))
    }
  }
  var_out <- var_out %||%
    ht_variance_dispatch(parts$design, parts$units, parts$pop, parts$rows, z,
                         parts$pi, variance, domain_frame(dom, used), ctx)
  out <- finish_estimate(value, var_out, parts, level, df, what = "contrast",
                         class = "drawn_contrast", alternative = alternative)
  if (any(!parts$respond)) {
    out$note <- paste(c(out$note, missing_note(sum(!parts$respond),
                                               length(parts$respond))),
                      collapse = " ")
  }
  out$weights <- stats::setNames(w[used], labels)
  out$estimates <- stats::setNames(est$est[used], labels)
  out$what <- what
  out$estimator <- if (what == "mean") estimator
  out$by <- dom$by
  out
}

#' A contrast weighs `by` domains, so it needs `by`
#' @noRd
check_contrast_by <- function(contrast, by) {
  if (!is.null(contrast) && is.null(by)) {
    stop("`contrast` weighs the domains of `by`, and `by` is NULL. Pass the ",
         "column or columns that define the domains as `by`.", call. = FALSE)
  }
  invisible(NULL)
}

#' Check a contrast against the domains, and spread it over them
#'
#' Returns one weight per domain in `dom$levels_in`, 0 for a domain the
#' contrast leaves out.
#'
#' @noRd
contrast_weights <- function(contrast, dom) {
  ticks <- function(x) paste0("`", x, "`", collapse = ", ")
  if (!is.numeric(contrast) || !length(contrast) || is.null(names(contrast))) {
    stop("`contrast` must be a named numeric vector of weights, one per ",
         "`by` domain, not ", format_bad(contrast), ". Name each weight ",
         "after its domain, as in `contrast = c(a = 1, b = -1)`.",
         call. = FALSE)
  }
  nm <- names(contrast)
  if (anyNA(nm) || !all(nzchar(nm))) {
    stop("`contrast` has a weight with no name. Name each weight after the ",
         "`by` domain it applies to.", call. = FALSE)
  }
  twice <- unique(nm[duplicated(nm)])
  if (length(twice)) {
    stop("`contrast` names ", ticks(twice), " more than once. Give each ",
         "domain one weight.", call. = FALSE)
  }
  bad <- which(!is.finite(contrast))
  if (length(bad)) {
    stop("`contrast` gives ", ticks(nm[bad[1]]), " a weight of ",
         format(contrast[[bad[1]]]), ". Give every domain a finite weight.",
         call. = FALSE)
  }
  shown <- dom$labels[dom$levels_in]
  at <- match(nm, shown)
  if (anyNA(at)) {
    miss <- nm[is.na(at)]
    several <- length(miss) > 1L
    empty <- miss[miss %in% dom$labels]
    if (length(empty)) {
      stop("`contrast` names ", ticks(empty), ", which ",
           if (length(empty) > 1L) "have" else "has", " no sampled row and ",
           "so no estimate. Leave ", if (length(empty) > 1L) "them" else "it",
           " out of `contrast`.", call. = FALSE)
    }
    stop("`contrast` names ", ticks(miss), ", which ",
         if (several) "are not domains" else "is not a domain", " of `by`. ",
         "Name the weights after the domains, which are ",
         ticks(utils::head(shown, 5L)),
         if (length(shown) > 5L) ", ..." else "", ".", call. = FALSE)
  }
  if (all(contrast == 0)) {
    stop("`contrast` gives every domain a weight of 0. Give at least one ",
         "domain a weight other than 0.", call. = FALSE)
  }
  w <- numeric(length(shown))
  w[at] <- as.numeric(contrast)
  w
}

#' @export
print.drawn_by <- function(x, ...) {
  what <- attr(x, "what")
  if (is.null(what)) {
    # Selecting columns with `[` or subset() keeps the class but drops the
    # attributes the heading reads, and the columns it describes may be gone.
    print(plain_frame(x), ...)
    return(invisible(x))
  }
  label <- if (what == "ratio") {
    paste("Ratio", ratio_name(attr(x, "variable"), "y"), "/",
          ratio_name(attr(x, "x"), "x"))
  } else if (what == "total") {
    if (is.null(attr(x, "x"))) "Horvitz-Thompson total" else
      paste0("Ratio estimate of the total, to `", attr(x, "x"), "`,")
  } else if (identical(attr(x, "estimator"), "ht")) {
    "Horvitz-Thompson mean"
  } else {
    "Hajek mean"
  }
  rows <- if (!several_y(x)) "domain" else
    if (length(attr(x, "by")) > 1L) "variable and domain" else "variable"
  scale <- unique(attr(x, "ci_scale") %||% "linear")
  shown <- if (length(scale) > 1L) {
    paste0(setdiff(scale, "linear")[1], " for proportions, ")
  } else {
    scale_label(scale)
  }
  cat(label, " by ", rows, "  (", attr(x, "design"), " design, ",
      ci_label(attr(x, "level"), attr(x, "alternative")), ", ",
      shown, df_label(attr(x, "df")), ")\n", sep = "")
  # The standard error is shown; its square adds width and nothing else.
  y <- plain_frame(x)
  y$variance <- NULL
  print(y, digits = 4, row.names = FALSE)
  print_unreachable(attr(x, "unreachable"), "these estimates")
  notes <- attr(x, "notes")
  notes <- notes[nzchar(notes) & !grepl("^Variance not requested", notes)]
  if (length(notes)) cat_note(notes)
  invisible(x)
}

#' @export
print.drawn_contrast <- function(x, ...) {
  what <- if (identical(x$what, "total")) "Horvitz-Thompson totals" else
    if (identical(x$estimator, "ht")) "Horvitz-Thompson means" else
      "Hajek means"
  est <- vapply(x$estimates, format, character(1), big.mark = ",")
  print_estimate(x, "contrast", paste("Contrast of", what), lines = c(
    contrast = contrast_label(x$weights),
    domains = paste0("`", names(est), "` ", est, collapse = ", ")
  ))
}

#' A contrast's weights as a sum, such as "`a` - `b`", or "a - b" without
#' `quote`
#' @noRd
contrast_label <- function(w, quote = TRUE) {
  tick <- if (quote) "`" else ""
  terms <- vapply(seq_along(w), function(i) {
    a <- abs(w[[i]])
    paste0(if (a != 1) paste0(format(a), " "), tick, names(w)[i], tick)
  }, character(1))
  out <- paste0(ifelse(w < 0, "- ", "+ "), terms, collapse = " ")
  sub("^- ", "-", sub("^\\+ ", "", out))
}

#' "95% CI", or "95% upper bound" for a one-sided interval
#' @noRd
ci_label <- function(level, alternative) {
  paste0(format(100 * level), "% ", switch(alternative %||% "two.sided",
    less = "upper bound",
    greater = "lower bound",
    "CI"
  ))
}

#' The interval's scale, as it leads the bracket after an interval
#' @noRd
scale_label <- function(scale) {
  switch(scale %||% "linear", logit = "logit, ", beta = "beta, ", "")
}

#' A line for frame rows the design could never select
#'
#' `draw()` warns about them by name; this keeps them in view on the
#' estimate, which says nothing about them.
#'
#' @noRd
print_unreachable <- function(n, what) {
  if (!isTRUE(n > 0)) return(invisible(NULL))
  cat("  unreachable  ", fmt_n(n), " frame row", if (n > 1) "s",
      " at probability 0, not covered by ", what, "\n", sep = "")
}

#' Design effect
#'
#' How much precision the design costs against simple random sampling of the
#' same size (Kish 1965). `deff = 1` means the design is doing as well as a
#' simple random sample from the rows it can reach; above 1 it is doing worse,
#' which is the usual price of clustering; below 1 it is doing better, which is
#' what stratification and probability-proportional-to-size buy you.
#'
#' Read it as an exchange rate on sample size: at `deff = 2`, a sample of 400
#' carries about as much information as 200 drawn at random.
#'
#' The simple random sample it is measured against is drawn from the rows the
#' design gives a probability above 0, `sum(inclusion_prob(design, frame) > 0)`
#' of them, because those are the rows a total or a Hajek mean estimates
#' over. Rows past `max_items`, outside a region or time window, dropped by
#' `na_rm` or in a stratum allocated no rows are left out. Counting them set
#' an equal-probability sample of 20 from the first 100 rows of a 400-row
#' frame at a design effect of 0.05 for a total, "better than simple random
#' sampling", where the right figure is 1.
#' `survey` measures against the sum of the weights, which is that count for
#' an equal-probability design, and gives the same design effect there. A
#' Horvitz-Thompson mean
#' (`estimator = "ht"`) divides the total by the frame's row count, so its
#' design effect is the total's.
#'
#' The design effect is `NA` when there is no variance, and when `y` takes
#' one value in every sampled row, where there is no spread for a simple
#' random sample to be measured by.
#'
#' A table of domain estimates from `by`, or of several variables, carries
#' no design effect, and `deff()` refuses it with a message saying so. A
#' contrast's design effect is `NA`. A ratio from [ht_ratio()], and a total
#' from the ratio estimator, are measured against the same estimator from a
#' simple random sample, through the variance of its linearised residuals,
#' as `survey::svyratio(deff = TRUE)` measures a ratio.
#'
#' `survey` has a `deff()` of its own, and whichever of the two packages is
#' attached last masks the other's. Either one answers for both. This
#' `deff()` is a generic that passes an object it has no method for, such as
#' `survey::svymean(..., deff = TRUE)`, to `survey::deff()`, and `survey`'s
#' generic finds the methods for the estimates made here.
#'
#' @param object A result from [ht_total()], [ht_mean()] or [ht_ratio()], or
#'   an object from `survey` that holds a design effect.
#' @param ... Passed to `survey::deff()` for an object from `survey`, such as
#'   its `quietly`. Not used for an estimate made here.
#'
#' @return A single number, or `NA` when the design has no variance estimate
#'   or `y` is constant in the sample. For an object from `survey`, what
#'   `survey::deff()` returns.
#'
#' @references
#' Kish, L. (1965). *Survey Sampling*. Wiley.
#'
#' @examples
#' set.seed(1)
#' # Sites differ from each other, and rows within a cluster are alike --
#' # exactly the structure that makes stratifying pay and clustering cost.
#' pop <- data.frame(
#'   id = 1:400,
#'   site = rep(c("a", "b", "c", "d"), each = 100),
#'   cl = rep(paste0("c", 1:40), each = 10)
#' )
#' pop$y <- rep(c(20, 60, 120, 200), each = 100) + round(stats::rnorm(400, 0, 8))
#'
#' # Stratifying on something that matters buys precision (deff below 1)
#' deff(ht_total(draw(pop, design_stratified("site", n = 40), seed = 1,
#'                    weights = TRUE), "y"))
#'
#' # Clustering usually costs it (deff above 1)
#' deff(ht_total(draw(pop, design_cluster("cl", n_clusters = 4), seed = 1,
#'                    weights = TRUE), "y"))
#'
#' @seealso [ht_total()], [ht_mean()], [plan_size()], which takes a design
#'   effect as an input.
#' @export
deff <- function(object, ...) UseMethod("deff")

# Each method is registered twice: for the generic above, and for `survey`'s,
# so whichever `deff()` is found first answers for an estimate made here.

#' @export
#' @exportS3Method survey::deff
deff.drawn_ht <- function(object, quietly = FALSE, ...) object$deff

#' @export
#' @exportS3Method survey::deff
deff.drawn_mean <- function(object, quietly = FALSE, ...) object$deff

#' @export
#' @exportS3Method survey::deff
deff.drawn_contrast <- function(object, quietly = FALSE, ...) object$deff

#' @export
#' @exportS3Method survey::deff
deff.drawn_ratio <- function(object, quietly = FALSE, ...) object$deff

#' @export
#' @exportS3Method survey::deff
deff.drawn_by <- function(object, quietly = FALSE, ...) {
  if (several_y(object)) {
    stop("`object` is a table of estimates for several variables, which ",
         "carries no design effect. Call deff() on ht_total() or ht_mean() ",
         "for one variable without `by`.", call. = FALSE)
  }
  stop("`object` is a table of domain estimates from `by`, which carries ",
       "no design effect. Call deff() on ht_total() or ht_mean() without ",
       "`by` for the design effect of the whole sample.", call. = FALSE)
}

#' @export
deff.default <- function(object, ...) {
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("`object` has class `", class(object)[1], "`, and deff() reads ",
         "the design effect of a result from ht_total() or ht_mean(). For an ",
         "object from `survey`, install it with install.packages(\"survey\").",
         call. = FALSE)
  }
  # survey::deff() dispatches from where it is called. Called here, it would
  # find this method before its own default and call it again, without end;
  # called from the base environment, it finds only its own methods.
  do.call(survey::deff, c(list(object), list(...)), envir = baseenv())
}

#' Shared argument handling for the estimators
#'
#' `na_rm` is `NULL` for ht_total(), which refuses a missing `y`, and the
#' `na_rm` of ht_mean() otherwise. With `na_rm = TRUE` a missing `y` is
#' allowed: `respond` marks the rows with `y` present, and `y` holds 0 in
#' the others, so a sum over every row is a sum over those rows.
#'
#' The work is split in three so that estimate_several() and ht_ratio() can
#' read the sample's state once and then one or more variables:
#' `prepare_sample()` reads the state and checks `level` and `df`, `add_y()`
#' reads one `y`, and `finish_prepare()` checks `.prob` and `.weight`.
#'
#' @noRd
ht_prepare <- function(sample, y, level, df = NULL, na_rm = FALSE) {
  parts <- prepare_sample(sample, level, df)
  finish_prepare(add_y(parts, y, na_rm))
}

#' The sample's state, before any variable is read
#' @noRd
prepare_sample <- function(sample, level, df = NULL) {
  state <- drawn_state(sample)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  if (nrow(sample) == 0L) {
    # Possible under Poisson sampling or a systematic start past the end, and
    # then a total of 0 (see empty_sample()). A design that reaches no row of
    # the frame, and a sample that records no probabilities, give nothing to
    # estimate from.
    if (is.null(state$pi_all)) {
      stop("`sample` has no rows, so there is nothing to estimate from.",
           call. = FALSE)
    }
    if (!any(reachable(state$pi_all))) {
      stop("`sample` has no rows, and its design gives every row of the ",
           "frame a probability of 0, so there is nothing to estimate from. ",
           "Give the design a sample size above 0.", call. = FALSE)
    }
  }
  if (!is.null(df) && (!is.numeric(df) || length(df) != 1L || is.na(df) ||
                       df < 1)) {
    stop("`df` must be NULL, Inf or a single number of 1 or more, not ",
         format_bad(df), ". A t interval needs at least one degree of ",
         "freedom.", call. = FALSE)
  }
  # The rows the design can reach: those with a probability above 0 at draw
  # time, or none recorded because the sample is conditional on clusters
  # that left them out. A sample saved by drawn 0.2.0 has no such record, so
  # its count is read from the design's columns and the unreachable rows are
  # unknown.
  pi_all <- state$pi_all
  list(design = state$design, rows = state$rows, pop = state$pop,
       units = state$units, pi_all = pi_all, pi = sample$.prob,
       sample = sample, compact = state$compact,
       reach = if (is.null(pi_all)) {
         reachable_rows(state$design, state$pop)
       } else {
         sum(reachable(pi_all))
       },
       unreachable = if (is.null(pi_all)) {
         NA_integer_
       } else {
         sum(!reachable(pi_all))
       })
}

#' Read one `y` into the parts of a sample
#'
#' `label` names `y` in the messages: "`y`", or the column with the `y` it
#' came from when `y` names several.
#'
#' @noRd
add_y <- function(parts, y, na_rm, label = "`y`") {
  v <- read_variable(parts$sample, y, "y", several = TRUE)
  yv <- v$values
  missing_y <- is.na(yv) & !is.nan(yv)
  if (any(missing_y) && !isTRUE(na_rm)) {
    stop(missing_y_message(parts$sample, missing_y, na_rm, label),
         call. = FALSE)
  }
  if (length(missing_y) && all(missing_y)) {
    stop(label, " is missing in every one of the ", length(missing_y), " ",
         "sampled rows, so there is no mean to estimate.", call. = FALSE)
  }
  need <- "A total or mean needs a finite `y` for every sampled row."
  check_finite_values(yv, missing_y, label, need)
  yv[missing_y] <- 0
  parts$y <- yv
  parts$respond <- !missing_y
  parts$y_name <- v$name
  parts$y_logical <- v$logical
  parts
}

#' Check `.prob` and `.weight` once the variables are read
#' @noRd
finish_prepare <- function(parts) {
  # A sample saved by drawn 0.2.0 records no probabilities for drawn_state()
  # to compare `.prob` with, so only its type and range are checked.
  if (is.null(parts$pi_all)) check_prob(parts$sample)
  warn_weight_edited(parts$sample)
  parts
}

#' A variable given as a column name or as a vector
#'
#' Returns its `values`, numeric with a logical turned to 0 and 1, whether
#' it was `logical`, and its `name`, `NULL` for a vector. `several` says
#' whether `arg` may also name several columns, which the caller reads one
#' at a time.
#'
#' @noRd
read_variable <- function(sample, y, arg, several = FALSE) {
  name <- NULL
  v <- y
  if (is.character(y) && length(y) == 1L && !is.na(y)) {
    if (!y %in% names(sample)) {
      stop("`sample` has no column `", y, "`.", call. = FALSE)
    }
    name <- y
    v <- sample[[y]]
    if (!is.numeric(v) && !is.logical(v)) {
      stop("`", arg, "` names `", y, "`, a ", class(v)[1], " column of ",
           "`sample`. Name a numeric or logical column.", call. = FALSE)
    }
  }
  if (!(is.numeric(v) || is.logical(v)) || length(v) != nrow(sample)) {
    stop(variable_forms(y, arg, several, sample), call. = FALSE)
  }
  logical <- is.logical(v)
  if (logical) v <- as.numeric(v)
  list(values = v, logical = logical, name = name)
}

#' What `arg` must be, and what it was
#' @noRd
variable_forms <- function(y, arg, several, sample) {
  forms <- if (several) {
    "one column name, several column names, or a numeric vector"
  } else {
    "one column name or a numeric vector"
  }
  a <- function(word) paste(if (grepl("^[aeiou]", word)) "an" else "a", word)
  what <- if (is.null(y)) {
    "NULL"
  } else if (is.character(y) && length(y) > 1L) {
    paste0("a character vector of ", length(y), " values, ",
           if (anyNA(y) || !all(nzchar(y))) {
             "some of them missing or empty"
           } else {
             "none of them a column of `sample`"
           })
  } else if (is.atomic(y)) {
    paste0(a(class(y)[1]), " vector of length ", length(y),
           if (length(y) != nrow(sample)) {
             paste0(", and `sample` has ", nrow(sample), " rows")
           })
  } else {
    a(class(y)[1])
  }
  paste0("`", arg, "` must be ", forms, " as long as `sample`. It is ",
         what, ".")
}

#' Refuse a variable that is infinite or `NaN` in some sampled row
#'
#' `label` names the variable and `need` is the sentence saying what needs
#' it finite.
#'
#' @noRd
check_finite_values <- function(v, missing, label, need) {
  bad <- which(!is.finite(v) & !missing)
  if (length(bad)) {
    stop(label, " is ", paste(unique(as.character(v[bad])), collapse = " or "),
         " in ", length(bad), " row(s) of `sample`: ", format_rows(bad),
         ". ", need, call. = FALSE)
  }
  invisible(NULL)
}

#' Why a missing `y` is refused, and what to do instead
#'
#' A total over the rows with `y` present counts the others as 0, so it is
#' low by about their share of the weight. ht_total() (`na_rm = NULL`) says
#' so; ht_mean() points to its `na_rm`.
#'
#' @noRd
missing_y_message <- function(sample, missing, na_rm, label = "`y`") {
  head <- paste0(label, " is missing in ", sum(missing), " of the ",
                 length(missing), " sampled rows. ")
  see <- "See \"Missing values and nonresponse\" in ?ht_total."
  if (!is.null(na_rm)) {
    return(paste0(head, "Set `na_rm = TRUE` to estimate the mean over the ",
                  "rows with `y` present, as a domain of the sample, which ",
                  "assumes the values are missing completely at random. ",
                  see))
  }
  w <- 1 / sample$.prob
  share <- sum(w[missing]) / sum(w)
  paste0(head, "A total over the rows with `y` present counts the others ",
         "as 0, so it is low by about their share of the weight",
         if (is.finite(share)) {
           paste0(", ", format(signif(100 * share, 2)), "% here")
         }, ". ht_mean(na_rm = TRUE) estimates the mean over the rows with ",
         "`y` present, with `by` for domains, as survey's `na.rm = TRUE` ",
         "does after as_svydesign(). For a total, adjust the weights for ",
         "nonresponse within weighting classes and estimate from them with ",
         "as_svydesign(). ", see)
}

#' Warn when `.weight` is no longer `1 / .prob`
#'
#' ht_total() and ht_mean() estimate from the design probabilities in `.prob`,
#' and as_svydesign() hands `survey` the `.weight` column. After a nonresponse
#' adjustment or trimming written into `.weight`, the two give different
#' answers, and the reason should be said once rather than found by comparing
#' them.
#'
#' @noRd
warn_weight_edited <- function(sample) {
  w <- sample[[".weight"]]
  if (is.null(w) || !nrow(sample)) return(invisible(NULL))
  design_w <- 1 / sample$.prob
  edited <- if (is.numeric(w)) sum(!same_prob(w, design_w)) else nrow(sample)
  if (edited == 0L) return(invisible(NULL))
  warn_drawn(paste0(
    "`sample` has `.weight` values that differ from 1 / `.prob` in ", edited,
    " row(s). ht_total() and ht_mean() use the design probabilities in ",
    "`.prob` and ignore `.weight`; estimate from adjusted weights with ",
    "as_svydesign()."), "drawn_weight_edited")
}

#' Warn that dplyr groups play no part in an estimate
#'
#' ht_total() and ht_mean() estimate over the whole sample, and so would a
#' grouped sample without this: the groups are ignored and one overall figure
#' comes back where the caller may expect one per group.
#'
#' @noRd
warn_grouped <- function(sample, fun) {
  if (!inherits(sample, "grouped_df")) return(invisible(NULL))
  vars <- setdiff(names(attr(sample, "groups", exact = TRUE)), ".rows")
  hint <- if (length(vars)) {
    paste0("`by = ", if (length(vars) > 1L) "c(",
           paste0("\"", vars, "\"", collapse = ", "),
           if (length(vars) > 1L) ")", "`")
  } else {
    "`by`"
  }
  warn_drawn(paste0(
    "`sample` is grouped with dplyr::group_by(), which ", fun, "() ignores, ",
    "so this is one estimate over the whole sample. Use ", hint, " for an ",
    "estimate per group."), "drawn_grouped_df")
}

#' A confidence interval from an estimate, a variance and degrees of freedom
#'
#' `scale` is `"linear"`, `"logit"` or `"beta"`; `logit = TRUE` is short for
#' the second. `alternative` gives a two-sided interval or a one-sided bound,
#' `"less"` an upper bound and `"greater"` a lower one, as in
#' [stats::t.test()]. The open end of a bound is `-Inf` or `Inf`, or 0 or 1
#' for a proportion. A one-sided bound at `level` is the end of a two-sided
#' interval at `1 - 2 * (1 - level)`.
#'
#' The beta interval (Korn and Graubard 1998) is Clopper-Pearson's at an
#' effective sample size of `p(1 - p) / v`, capped at `n` and scaled by
#' `(qt(a, n - 1) / qt(a, df))^2`, the design's degrees of freedom against
#' those of a simple random sample of `n`. With `a` the tail probability this
#' is `survey::svyciprop(method = "beta")` when the cap does not bind. At an
#' estimate of 0 or 1 the variance is 0 and says nothing, and the effective
#' sample size is `n`. `n` counts the sampled rows below certainty behind the
#' estimate; with none, the estimate is exact.
#'
#' @noRd
interval <- function(est, v, level, df, logit = FALSE,
                     alternative = "two.sided",
                     scale = if (logit) "logit" else "linear", n = NA_real_) {
  if (!is.finite(v) || v < 0) return(c(NA_real_, NA_real_))
  proportion <- scale %in% c("logit", "beta")
  open <- if (proportion) c(0, 1) else c(-Inf, Inf)
  bound <- function(lo, hi) {
    switch(alternative, less = c(open[1], hi), greater = c(lo, open[2]),
           c(lo, hi))
  }
  edge <- est <= 0 || est >= 1
  # A variance of 0 is an exact estimate, whatever the degrees of freedom,
  # except for a proportion of 0 or 1, whose variance is 0 because no row of
  # the other kind was drawn.
  exact <- v == 0 && !(scale == "beta" && edge && isTRUE(n >= 1))
  if (exact) return(bound(est, est))
  if (is.finite(df) && df < 1) return(c(NA_real_, NA_real_))
  # The tail probability on each side that has one.
  a <- if (alternative == "two.sided") (1 - level) / 2 else 1 - level
  quantile <- function(p, df) {
    if (is.finite(df)) stats::qt(p, df) else stats::qnorm(p)
  }
  q <- quantile(1 - a, df)
  if (scale == "beta" && est >= 0 && est <= 1) {
    n_eff <- if (edge) n else min(n, est * (1 - est) / v)
    if (n >= 2) n_eff <- n_eff * (stats::qt(a, n - 1) / quantile(a, df))^2
    x <- n_eff * est
    lo <- if (x <= 0) 0 else stats::qbeta(a, x, n_eff - x + 1)
    hi <- if (x >= n_eff) 1 else stats::qbeta(1 - a, x + 1, n_eff - x)
    return(bound(lo, hi))
  }
  if (proportion && !edge) {
    # A proportion's interval on the logit scale, by the delta method, then
    # back: it cannot leave [0, 1], and near 0 or 1 it is asymmetric the way
    # the sampling distribution is. This is survey::svyciprop()'s "xlogit".
    se_logit <- sqrt(v) / (est * (1 - est))
    ends <- stats::plogis(stats::qlogis(est) + c(-1, 1) * q * se_logit)
    return(bound(ends[1], ends[2]))
  }
  if (proportion) open <- c(-Inf, Inf)
  bound(est - q * sqrt(v), est + q * sqrt(v))
}

#' The interval for a total or mean, and the scale it was formed on
#'
#' A mean that is not a proportion gets the symmetric interval. A proportion
#' gets `ci_method`'s, except that an estimate of exactly 0 or 1 gets the beta
#' interval whatever `ci_method` says: its variance is 0, and a logit interval
#' around it has no width. `pi_i` holds the probabilities of the sampled rows
#' behind the estimate, the domain's for a domain; those below certainty are
#' the beta interval's `n`, and with none of them the estimate is exact.
#'
#' `n` is that count, for a caller that holds it in place of `pi_i`.
#'
#' Returns the interval, its `scale`, `edge` (the beta interval was used for
#' an estimate of 0 or 1) and `exact`.
#'
#' @noRd
ci_for <- function(est, v, level, df, alternative = "two.sided",
                   proportion = FALSE, ci_method = "logit", pi_i = NULL,
                   n = sum(!is_certain(pi_i))) {
  if (!proportion) {
    return(list(ci = interval(est, v, level, df, alternative = alternative),
                scale = "linear", edge = FALSE, exact = FALSE))
  }
  edge <- is.finite(v) && v >= 0 && (est == 0 || est == 1) && n >= 1
  scale <- if (edge) "beta" else ci_method
  # A Horvitz-Thompson mean can pass 1. Its interval is then symmetric.
  if (est < 0 || est > 1) scale <- "linear"
  list(ci = interval(est, v, level, df, alternative = alternative,
                     scale = scale, n = n),
       scale = scale, edge = edge, exact = n == 0)
}

#' Why a proportion of 0 or 1 has a beta interval
#' @noRd
edge_note <- function(est = NULL, n = NULL, domains = NULL,
                      ci_method = "logit", certain = FALSE) {
  why <- if (ci_method == "logit") {
    ", and a logit interval around it would have no width"
  } else {
    ""
  }
  below <- if (certain) " below certainty" else ""
  if (is.null(domains)) {
    return(paste0(
      "The estimate is ", format(est), ", with an estimated variance of 0",
      why, ". The interval is Korn and Graubard's beta interval, with the ",
      fmt_n(n), " sampled row", if (n > 1) "s", below, " as its effective ",
      "sample size."))
  }
  several <- length(domains) > 1L
  paste0("The domain", if (several) "s", " ",
         paste0("`", domains, "`", collapse = ", "),
         if (several) " have" else " has",
         " a proportion of 0 or 1, with an estimated variance of 0. ",
         if (several) "Their intervals are" else "Its interval is",
         " Korn and Graubard's beta interval", if (several) "s",
         ", with the domain's sampled rows", below, " as the effective ",
         "sample size.")
}

#' Assemble an estimate, its interval, and its design effect
#'
#' `divide_by` is what a mean divides its Horvitz-Thompson total by, for the
#' design effect: see `deff_value()`. `variable` is the name of `y`, `NULL`
#' for a vector, and `n_random` the sampled rows below certainty, which
#' `confint()` needs to form a beta interval at another level.
#'
#' @noRd
finish_estimate <- function(est, var_out, parts, level, df, what, class,
                            y = NULL, proportion = FALSE, ci_method = "logit",
                            alternative = "two.sided",
                            divide_by = parts$reach) {
  v <- var_out$variance
  se <- if (is.na(v) || v < 0) NA_real_ else sqrt(v)
  ci <- ci_for(est, v, level, df, alternative, proportion, ci_method,
               parts$pi)
  note <- var_out$note
  if (!is.na(se) && anyNA(ci$ci)) {
    note <- paste(c(note, df_note(df)), collapse = " ")
  }
  if (ci$edge) {
    certain <- is_certain(parts$pi)
    note <- paste(c(note, edge_note(est, sum(!certain), ci_method = ci_method,
                                    certain = any(certain))),
                  collapse = " ")
  }
  n <- nrow(parts$sample)

  out <- list(variable = parts$y_name, variance = v, se = se, ci = ci$ci,
              level = level, alternative = alternative, df = df, n = n,
              n_random = sum(!is_certain(parts$pi)),
              ci_scale = ci$scale, design = design_type(parts$design),
              deff = if (is.null(y)) NA_real_ else
                deff_value(v, y, parts$pi, parts$reach, n, what, divide_by),
              method = var_out$method %||% "analytic", note = note,
              unreachable = parts$unreachable)
  out[[what]] <- est
  structure(out, class = class)
}

#' Why an estimate with a standard error has no interval
#'
#' `ht_prepare()` refuses a `df` below 1 from the caller, so a missing
#' interval with a finite standard error comes from the design's degrees of
#' freedom.
#'
#' @noRd
df_note <- function(df) {
  paste0("The design leaves ", format(df), " degrees of freedom, the number ",
         "of primary sampling units below certainty less the number of strata ",
         "holding them, so there is no t interval. Pass `df` to set the ",
         "degrees of freedom, or `df = Inf` for a normal interval.")
}

#' Variance under simple random sampling of the same size, for comparison
#'
#' The population variance is estimated from the sample with design weights, so
#' this works from a sample of any design.
#'
#' The denominator is `sum(w) * (n - 1) / n`, not `sum(w) - 1`. Under simple
#' random sampling `sum(w * (y - mu)^2)` has expectation `(N/n)(n-1)S^2`, so
#' dividing by `N - 1` leaves the estimate short by a factor of `n(N-1)/(N(n-1))`
#' -- and the design effect, which divides by it, long by the reciprocal. That
#' put a plain simple random sample of 10 at a design effect of 1.11, printed as
#' "worse than simple random sampling". This form is exactly unbiased there.
#'
#' `N` is the number of rows the design can reach, those with a probability
#' above 0. The simple random sample is drawn from them, since they are the
#' population a total or a Hajek mean estimates over; `survey` uses the sum of
#' the weights, which is the same count for an equal-probability design.
#' `divide_by` is what a mean divides its Horvitz-Thompson total by: `N` for
#' the Hajek mean, and the frame's row count for the Horvitz-Thompson mean,
#' whose design effect is then the total's.
#'
#' A constant `y` gives `NA`. Its spread is 0, and computed it came out as a
#' rounding error, which put a total's design effect at `1e32`.
#'
#' @noRd
deff_value <- function(v, y, pi_i, N, n, what, divide_by = N) {
  if (!is.finite(v) || v <= 0 || n < 2L || !is.finite(N) || N < 2L) {
    return(NA_real_)
  }
  if (all(y == y[1])) return(NA_real_)
  w <- 1 / pi_i
  mu <- sum(w * y) / sum(w)
  s2 <- sum(w * (y - mu)^2) * n / (sum(w) * (n - 1))
  if (!is.finite(s2) || s2 <= 0) return(NA_real_)
  fpc <- max(0, 1 - n / N)
  v_srs <- if (what == "mean") {
    (N / divide_by)^2 * s2 / n * fpc
  } else {
    N^2 * s2 / n * fpc
  }
  if (!is.finite(v_srs) || v_srs <= 0) return(NA_real_)
  v / v_srs
}

#' @export
print.drawn_mean <- function(x, ...) {
  lines <- if (isTRUE(x$n_missing > 0)) {
    c(missing = paste0(fmt_n(x$n_missing), " sampled row",
                       if (x$n_missing > 1) "s", " left out, assumed missing ",
                       "completely at random"))
  }
  print_estimate(x, "mean", if (identical(x$estimator, "ht"))
    "Horvitz-Thompson mean" else "Hajek mean", lines = lines)
}

#' @noRd
df_label <- function(df) {
  if (is.null(df) || !is.finite(df)) "normal" else
    paste0("t, ", format(signif(df, 4)), " df")
}

#' The variance method as printed: the value `method` holds, which a `by`
#' table shows in its `method` column, and whether it is an approximation
#' @noRd
method_label <- function(m) {
  m <- m %||% "analytic"
  if (m %in% c("deville", "successive difference", "local mean")) {
    paste0(m, ", an approximation")
  } else {
    m
  }
}

#' Print a single estimate
#'
#' `lines` are printed between the heading and the estimate, each after its
#' name.
#'
#' @noRd
print_estimate <- function(x, field, label, lines = NULL) {
  cat(label, "  (", x$design, " design, n = ", x$n, ")\n", sep = "")
  for (nm in names(lines)) {
    cat("  ", formatC(nm, width = -8), " ", lines[[nm]], "\n", sep = "")
  }
  cat("  estimate ", format(x[[field]], big.mark = ","), "\n", sep = "")
  if (is.na(x$se)) {
    cat("  se       NA\n")
  } else {
    cat("  se       ", format(x$se, big.mark = ","), "  (",
        method_label(x$method), ")\n", sep = "")
    if (!anyNA(x$ci)) {
      alt <- x$alternative %||% "two.sided"
      cat("  ", ci_label(x$level, alt), "  ",
          switch(alt,
            less = format(x$ci[2], big.mark = ","),
            greater = format(x$ci[1], big.mark = ","),
            paste0(format(x$ci[1], big.mark = ","), " to ",
                   format(x$ci[2], big.mark = ","))
          ), "  (", scale_label(x$ci_scale), df_label(x$df), ")\n", sep = "")
    }
    if (is.finite(x$deff)) {
      cat("  deff     ", signif(x$deff, 3), "  (",
          if (x$deff > 1.05) "worse than" else if (x$deff < 0.95) "better than"
          else "about the same as", " simple random sampling)\n", sep = "")
    }
  }
  print_unreachable(x$unreachable, "this estimate")
  # A standard error of 0 always has a reason worth stating: rows taken with
  # certainty, a constant `y`, or a proportion of 0 or 1.
  if (!is.null(x$note) && (is.na(x$se) || x$se == 0 ||
                            !identical(x$method, "analytic") ||
                            anyNA(x$ci))) {
    cat_note(x$note)
  }
  invisible(x)
}

#' Print a note, wrapped and indented, after one blank line
#' @noRd
cat_note <- function(text) {
  cat("\n", paste0(strwrap(text, prefix = "  "), collapse = "\n"), "\n",
      sep = "")
}

#' A `by` table without its classes and attributes
#' @noRd
plain_frame <- function(x) {
  attributes(x) <- attributes(x)[c("names", "row.names")]
  class(x) <- "data.frame"
  x
}

# ---- the estimate as numbers ----------------------------------------------
#
# coef(), vcov(), confint(), as.data.frame() and tidy() for the four kinds of
# estimate. survey::SE() and survey::cv() call vcov() and coef(), so they work
# through these.

#' The field that holds a single estimate
#' @noRd
estimate_field <- function(x) {
  if (inherits(x, "drawn_ht")) "total" else
    if (inherits(x, "drawn_mean")) "mean" else
      if (inherits(x, "drawn_ratio")) "ratio" else "contrast"
}

#' What coef() names each estimate
#'
#' A total or mean is named after `y`, or "total" or "mean" when `y` was a
#' vector; a ratio as "y/x", as `survey::svyratio()` names it; a contrast by
#' its weights, as in "a - b"; and each domain of a `by` table as `contrast`
#' names it, after its variable when the table holds several.
#'
#' @noRd
estimate_names <- function(x) {
  if (inherits(x, "drawn_by")) return(domain_names(x))
  if (inherits(x, "drawn_contrast")) {
    return(contrast_label(x$weights, quote = FALSE))
  }
  if (inherits(x, "drawn_ratio")) {
    return(paste0(x$variable %||% "y", "/", x$x %||% "x"))
  }
  x$variable %||% estimate_field(x)
}

#' Each row's domain, named as `contrast` names it, such as "a . east"
#' @noRd
domain_names <- function(x) {
  key <- group_key(x, attr(x, "by"))
  group_label(levels(key))[as.integer(key)]
}

#' Refuse a `by` table that lost the attributes the methods read
#'
#' Selecting columns with `[` or subset() keeps the class and drops them;
#' selecting rows keeps them.
#'
#' @noRd
check_by_table <- function(x) {
  what <- attr(x, "what")
  by <- attr(x, "by")
  if (is.null(what) || is.null(by) || !all(c(by, what) %in% names(x))) {
    stop("`object` does not carry the attributes of a `by` table, which ",
         "selecting columns with `[` or subset() drops. Use the table ",
         "ht_total() or ht_mean() returned, or select rows only.",
         call. = FALSE)
  }
  invisible(x)
}

#' @noRd
check_ci_level <- function(level, arg = "level") {
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("`", arg, "` must be a single number strictly between 0 and 1, not ",
         format_bad(level), ".", call. = FALSE)
  }
  level
}

#' Whether `level` is the level an interval was formed at
#' @noRd
same_level <- function(level, stored) {
  !is.null(stored) && abs(level - stored) < 1e-12
}

#' Intervals as confint() returns them: a row per estimate, and columns named
#' by the probability at each end
#' @noRd
ci_matrix <- function(lo, hi, names, level, alternative, parm = NULL) {
  p <- switch(alternative %||% "two.sided",
    less = c(0, level),
    greater = c(1 - level, 1),
    c((1 - level) / 2, 1 - (1 - level) / 2)
  )
  out <- matrix(c(lo, hi), ncol = 2L, dimnames = list(names, paste(
    format(100 * p, trim = TRUE, scientific = FALSE, digits = 3), "%")))
  if (is.null(parm)) return(out)
  if (is.numeric(parm)) parm <- names[parm]
  out[parm, , drop = FALSE]
}

#' The interval of a total, mean or contrast at `level`
#'
#' The stored one at its own level. At another, it is formed again by
#' `interval()` on the estimate's scale, with the design's degrees of freedom
#' and, for a beta interval, the sampled rows below certainty.
#'
#' @noRd
single_ci <- function(x, level) {
  if (same_level(level, x$level)) return(x$ci)
  scale <- x$ci_scale %||% "linear"
  if (scale == "beta" && is.null(x$n_random)) {
    stop("`object` does not hold the sample size its beta interval needs at ",
         "another level. Estimate again with `level = ", format(level), "`.",
         call. = FALSE)
  }
  interval(x[[estimate_field(x)]], x$variance, level, x$df,
           alternative = x$alternative %||% "two.sided", scale = scale,
           n = x$n_random %||% NA_real_)
}

#' The intervals of a `by` table at `level`, a row per domain
#' @noRd
by_ci <- function(x, level) {
  if (same_level(level, attr(x, "level"))) {
    return(cbind(x$ci_lower, x$ci_upper))
  }
  est <- x[[attr(x, "what")]]
  v <- if ("variance" %in% names(x)) x$variance else x$se^2
  scale <- row_scales(x)
  proportion <- scale %in% c("logit", "beta")
  n <- if (any(proportion)) unname(attr(x, "n_random")[domain_names(x)])
  if (any(proportion) && (is.null(n) || anyNA(n[proportion]))) {
    stop("`object` does not hold the sample sizes its intervals for a ",
         "proportion need at another level. Estimate again with `level = ",
         format(level), "`.", call. = FALSE)
  }
  out <- vapply(seq_along(est), function(i) {
    ci_for(est[i], v[i], level, attr(x, "df"),
           attr(x, "alternative") %||% "two.sided", proportion[i], scale[i],
           n = n[i])$ci
  }, numeric(2))
  t(out)
}

#' The interval scale of each row of a `by` table
#'
#' One value for the whole table, or, for a table of several variables of
#' which some are proportions, one per row, named by domain_names().
#'
#' @noRd
row_scales <- function(x) {
  scale <- attr(x, "ci_scale") %||% "linear"
  if (length(scale) == 1L) return(rep(scale, nrow(x)))
  unname(scale[domain_names(x)])
}

#' @export
coef.drawn_ht <- function(object, ...) {
  stats::setNames(object[[estimate_field(object)]], estimate_names(object))
}

#' @export
coef.drawn_mean <- function(object, ...) coef.drawn_ht(object)

#' @export
coef.drawn_contrast <- function(object, ...) coef.drawn_ht(object)

#' @export
coef.drawn_ratio <- function(object, ...) coef.drawn_ht(object)

#' @export
coef.drawn_by <- function(object, ...) {
  check_by_table(object)
  stats::setNames(object[[attr(object, "what")]], domain_names(object))
}

#' @export
vcov.drawn_ht <- function(object, ...) {
  nm <- estimate_names(object)
  matrix(object$variance %||% NA_real_, 1L, 1L, dimnames = list(nm, nm))
}

#' @export
vcov.drawn_mean <- function(object, ...) vcov.drawn_ht(object)

#' @export
vcov.drawn_contrast <- function(object, ...) vcov.drawn_ht(object)

#' @export
vcov.drawn_ratio <- function(object, ...) vcov.drawn_ht(object)

# The covariances between domains are not estimated, so the off-diagonal
# entries are 0. A difference of domains needs `contrast`. A table of
# several variables holds the covariances between them in its `vcov`
# attribute, matched to the rows by name, so that it follows a selection of
# rows.

#' @export
vcov.drawn_by <- function(object, ...) {
  check_by_table(object)
  nm <- domain_names(object)
  full <- attr(object, "vcov")
  if (!is.null(full) && all(nm %in% rownames(full))) {
    return(full[nm, nm, drop = FALSE])
  }
  v <- if ("variance" %in% names(object)) object$variance else object$se^2
  out <- diag(v, nrow = length(v))
  dimnames(out) <- list(nm, nm)
  out
}

#' @export
confint.drawn_ht <- function(object, parm, level = object$level, ...) {
  check_ci_level(level)
  ci <- single_ci(object, level)
  ci_matrix(ci[1], ci[2], estimate_names(object), level, object$alternative,
            if (!missing(parm)) parm)
}

#' @export
confint.drawn_mean <- function(object, parm, level = object$level, ...) {
  confint.drawn_ht(object, parm, level)
}

#' @export
confint.drawn_contrast <- function(object, parm, level = object$level, ...) {
  confint.drawn_ht(object, parm, level)
}

#' @export
confint.drawn_ratio <- function(object, parm, level = object$level, ...) {
  confint.drawn_ht(object, parm, level)
}

#' @export
confint.drawn_by <- function(object, parm, level = attr(object, "level"),
                             ...) {
  check_by_table(object)
  check_ci_level(level)
  ci <- by_ci(object, level)
  ci_matrix(ci[, 1], ci[, 2], domain_names(object), level,
            attr(object, "alternative"), if (!missing(parm)) parm)
}

#' One row for a total, mean or contrast
#' @noRd
single_frame <- function(x, row.names = NULL) {
  out <- data.frame(variable = x$variable %||% NA_character_,
                    stringsAsFactors = FALSE)
  if (inherits(x, "drawn_contrast")) {
    out$domains <- contrast_label(x$weights, quote = FALSE)
  }
  if (inherits(x, "drawn_ratio")) out$x <- x$x %||% NA_character_
  out$n <- x$n
  out[[estimate_field(x)]] <- x[[estimate_field(x)]]
  out$variance <- x$variance
  out$se <- x$se
  out$ci_lower <- x$ci[1]
  out$ci_upper <- x$ci[2]
  out$level <- x$level
  out$df <- x$df
  out$deff <- x$deff
  out$method <- x$method
  if (!is.null(row.names)) row.names(out) <- row.names
  out
}

#' @export
as.data.frame.drawn_ht <- function(x, row.names = NULL, optional = FALSE,
                                   ...) {
  single_frame(x, row.names)
}

#' @export
as.data.frame.drawn_mean <- function(x, row.names = NULL, optional = FALSE,
                                     ...) {
  single_frame(x, row.names)
}

#' @export
as.data.frame.drawn_contrast <- function(x, row.names = NULL,
                                         optional = FALSE, ...) {
  single_frame(x, row.names)
}

#' @export
as.data.frame.drawn_ratio <- function(x, row.names = NULL, optional = FALSE,
                                      ...) {
  single_frame(x, row.names)
}

#' @export
as.data.frame.drawn_by <- function(x, row.names = NULL, optional = FALSE,
                                   ...) {
  out <- plain_frame(x)
  if (!is.null(row.names)) row.names(out) <- row.names
  out
}

#' The columns of tidy(): `term`, `estimate`, `std.error`, and the interval
#' with `conf.int = TRUE`
#' @noRd
tidy_frame <- function(term, estimate, se, ci = NULL) {
  out <- data.frame(term = term, estimate = estimate, std.error = se,
                    stringsAsFactors = FALSE)
  if (!is.null(ci)) {
    out$conf.low <- ci[, 1]
    out$conf.high <- ci[, 2]
  }
  out
}

#' @exportS3Method generics::tidy
tidy.drawn_ht <- function(x, conf.int = FALSE, conf.level = x$level, ...) {
  ci <- if (isTRUE(conf.int)) {
    check_ci_level(conf.level, "conf.level")
    matrix(single_ci(x, conf.level), 1L)
  }
  tidy_frame(estimate_names(x), x[[estimate_field(x)]], x$se, ci)
}

#' @exportS3Method generics::tidy
tidy.drawn_mean <- function(x, conf.int = FALSE, conf.level = x$level, ...) {
  tidy.drawn_ht(x, conf.int, conf.level)
}

#' @exportS3Method generics::tidy
tidy.drawn_contrast <- function(x, conf.int = FALSE, conf.level = x$level,
                                ...) {
  tidy.drawn_ht(x, conf.int, conf.level)
}

#' @exportS3Method generics::tidy
tidy.drawn_ratio <- function(x, conf.int = FALSE, conf.level = x$level, ...) {
  tidy.drawn_ht(x, conf.int, conf.level)
}

#' @exportS3Method generics::tidy
tidy.drawn_by <- function(x, conf.int = FALSE,
                          conf.level = attr(x, "level"), ...) {
  check_by_table(x)
  what <- attr(x, "what")
  ci <- if (isTRUE(conf.int)) {
    check_ci_level(conf.level, "conf.level")
    by_ci(x, conf.level)
  }
  by <- attr(x, "by")
  term <- if (several_y(x)) {
    by <- setdiff(by, "variable")
    x$variable
  } else if (what == "ratio") {
    rep(paste0(attr(x, "variable") %||% "y", "/", attr(x, "x") %||% "x"),
        nrow(x))
  } else {
    rep(attr(x, "variable") %||% what, nrow(x))
  }
  out <- tidy_frame(term, x[[what]], x$se, ci)
  if (!length(by)) return(out)
  cbind(plain_frame(x)[by], out)
}
