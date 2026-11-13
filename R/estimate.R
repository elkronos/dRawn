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
#' @param sample A data frame returned by [draw()] with `weights = TRUE`. It is
#'   read as [ht_total()] reads it: `.prob` for the probabilities, and the
#'   stored frame for the design's columns. See "Which columns are read" there.
#' @param y The variable to average: a column name, or a numeric vector as long
#'   as `sample`.
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
#'
#' @return A list with a `print()` method, holding `mean`, `estimator`,
#'   `variance`, `se`, `ci`, `level`, `alternative`, `df`, `n`, `design`,
#'   `deff`, `method`, `note`, `ci_scale` (`"linear"`, `"logit"` or `"beta"`)
#'   and `unreachable`. See [ht_total()] for what `method`, `df`, `note` and
#'   `unreachable` say, and [deff()] for the design effect. With `by`, a data
#'   frame with one row per domain, and with `contrast` a `drawn_contrast`, as
#'   for [ht_total()].
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
#' @seealso [ht_total()], [deff()], [sample_summary()]
#' @export
ht_mean <- function(sample, y, estimator = c("hajek", "ht"),
                    variance = c("auto", "analytic", "jackknife", "none"),
                    level = 0.95, by = NULL, df = NULL, proportion = NULL,
                    ci_method = c("logit", "beta"),
                    alternative = c("two.sided", "less", "greater"),
                    contrast = NULL) {
  estimator <- check_choice(estimator, c("hajek", "ht"), "estimator")
  variance <- check_choice(variance,
                           c("auto", "analytic", "jackknife", "none"),
                           "variance")
  ci_method <- check_choice(ci_method, c("logit", "beta"), "ci_method")
  alternative <- check_choice(alternative, c("two.sided", "less", "greater"),
                              "alternative")
  check_contrast_by(contrast, by)
  parts <- ht_prepare(sample, y, level, df)
  parts$proportion <- mean_is_proportion(parts, proportion)
  if (ci_method == "beta" && !parts$proportion) {
    stop("`ci_method = \"beta\"` gives an interval for a proportion, and ",
         "this mean is not one. Pass `proportion = TRUE` if `y` is a 0/1 ",
         "indicator; see \"Proportions\" in ?ht_mean.", call. = FALSE)
  }
  if (!is.null(by)) {
    return(estimate_by(parts, by, variance, level, df, what = "mean",
                       estimator = estimator, ci_method = ci_method,
                       alternative = alternative, contrast = contrast))
  }
  warn_grouped(sample, "ht_mean")
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
                      divide_by = denom_rows(parts, estimator))
  if (constant && identical(out$variance, 0) &&
      !identical(out$ci_scale, "beta") && any(!is_certain(pi_i))) {
    out$note <- paste(c(out$note, constant_note(yv[1], parts$proportion)),
                      collapse = " ")
  }
  out
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
  y <- parts$y
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
  if (is.null(col) || !col %in% names(parts$pop)) return(FALSE)
  x <- parts$pop[[col]]
  if (is.logical(x)) return(TRUE)
  if (!is.numeric(x)) return(FALSE)
  x <- x[!is.na(x)]
  length(x) > 0L && all(x == 0 | x == 1)
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
#' @noRd
estimate_by <- function(parts, by, variance, level, df, what,
                        estimator = "hajek", ci_method = "logit",
                        alternative = "two.sided", contrast = NULL) {
  dom <- domain_keys(parts, by)
  est <- domain_estimates(parts, dom, what, estimator)
  pi_i <- parts$pi
  df_use <- df %||% design_df(parts$design, parts$units, parts$pop,
                              parts$rows, pi_i)
  ctx <- variance_context()
  if (!is.null(contrast)) {
    return(estimate_contrast(parts, dom, est, contrast, variance, level,
                             df_use, what, estimator, alternative, ctx))
  }
  vars <- domain_variances(parts, dom, est, variance,
                           ratio = what == "mean" && estimator == "hajek",
                           ctx = ctx)
  proportion <- what == "mean" && isTRUE(parts$proportion)
  labels <- dom$labels
  levels_in <- dom$levels_in

  rows <- lapply(seq_along(levels_in), function(j) {
    v <- vars[[j]]
    ci <- ci_for(est$est[j], v$variance, level, df_use, alternative,
                 proportion, ci_method, pi_i[dom$idx[[j]]])
    list(se = if (is.finite(v$variance) && v$variance >= 0)
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
  out$n <- as.numeric(lengths(dom$idx))
  out[[what]] <- est$est
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
    notes <- c(notes, constant_note(NA, TRUE,
                                    domains = labels[levels_in[constant]]))
  }
  structure(out, class = c("drawn_by", "data.frame"),
            what = what, estimator = if (what == "mean") estimator,
            level = level, alternative = alternative, df = df_use,
            design = design_type(parts$design),
            ci_scale = if (proportion) ci_method else "linear",
            notes = notes, unreachable = parts$unreachable)
}

#' The domains of `by`, keyed in the sample and, where it can be, the frame
#'
#' Returns `key`, each sampled row's domain as a code into `labels`;
#' `pop_key`, each frame row's, or `NULL` when the frame lacks a `by` column;
#' `levels_in`, the codes of the domains with a sampled row, in key order;
#' `pos`, each sampled row's position in `levels_in`; `idx`, the sampled rows
#' of each domain in `levels_in`; and `aligned`, whether the frame's `by`
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
  levels_in <- which(tabulate(key, length(labels)) > 0L)
  pos <- match(key, levels_in)
  idx <- split(seq_along(key),
               structure(pos, levels = as.character(seq_along(levels_in)),
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
#' Horvitz-Thompson mean, `y / N_d` with `N_d` the domain's rows in the frame.
#' Each is 0 outside its domain, so one vector `z` holds them all: each
#' sampled row's value in its own domain. A Hajek mean of a `y` constant over
#' the domain has residuals of exactly 0, not the rounding error the
#' subtraction leaves, and `constant` marks it.
#'
#' @noRd
domain_estimates <- function(parts, dom, what, estimator) {
  yv <- parts$y
  pi_i <- parts$pi
  k <- length(dom$levels_in)
  est <- numeric(k)
  constant <- logical(k)
  z <- numeric(length(yv))
  size <- NULL
  if (what == "mean" && estimator == "ht") {
    if (is.null(dom$pop_key)) {
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
#' A Hajek mean (`ratio`) over one primary sampling unit is declined first,
#' with `ratio_decline()`. Where the design's analytic estimator is
#' `srs_variance()` over strata, intervals or one group, every other domain is
#' estimated in one pass with `srs_variance_by()`. A domain that pass cannot
#' settle, because a lonely group holds some of its units or its variance
#' overflowed, goes through `ht_variance_dispatch()` as every domain of the
#' other designs does, with the shared `ctx`, and gets the same note and
#' fallback as an estimate of its own.
#'
#' @noRd
domain_variances <- function(parts, dom, est, variance, ratio, ctx) {
  k <- length(dom$levels_in)
  out <- vector("list", k)
  if (ratio && variance != "none") {
    psu <- sample_units(parts$design, parts$units, parts$pop, parts$rows)
    for (j in seq_len(k)) {
      out[j] <- list(ratio_decline(parts, dom$idx[[j]], domain_frame(dom, j),
                                   dom$labels[dom$levels_in[j]], psu))
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
  v <- srs_variance_by(z[keep], lay$group[keep], dom$pos[keep], k, lay$n_g,
                       lay$N_g)
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
  z <- w[dom$pos] * est$z
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
  label <- if (what == "total") "Horvitz-Thompson total" else
    if (identical(attr(x, "estimator"), "ht")) "Horvitz-Thompson mean" else
      "Hajek mean"
  cat(label, " by domain  (", attr(x, "design"), " design, ",
      ci_label(attr(x, "level"), attr(x, "alternative")), ", ",
      scale_label(attr(x, "ci_scale")), df_label(attr(x, "df")), ")\n",
      sep = "")
  y <- x
  attributes(y) <- attributes(x)[c("names", "row.names")]
  class(y) <- "data.frame"
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

#' A contrast's weights as a sum, such as "`a` - `b`"
#' @noRd
contrast_label <- function(w) {
  terms <- vapply(seq_along(w), function(i) {
    a <- abs(w[[i]])
    paste0(if (a != 1) paste0(format(a), " "), "`", names(w)[i], "`")
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
#' `survey` has a `deff()` of its own, and attaching `survey` after this
#' package masks this one. Both work on an estimate from [ht_total()] or
#' [ht_mean()], so the masking is harmless.
#'
#' @param x A result from [ht_total()] or [ht_mean()].
#'
#' @return A single number, or `NA` when the design has no variance estimate
#'   or `y` is constant in the sample.
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
deff <- function(x) {
  if (!inherits(x, c("drawn_ht", "drawn_mean"))) {
    stop("`x` must come from ht_total() or ht_mean().", call. = FALSE)
  }
  x$deff
}

# `survey` exports a `deff()` generic too, and attaching `survey` after this
# package masks ours. Registering methods for its generic means whichever
# `deff()` is found first still answers for an estimate made here, instead of
# failing with "non-numeric argument to mathematical function".

#' @exportS3Method survey::deff
deff.drawn_ht <- function(object, quietly = FALSE, ...) object$deff

#' @exportS3Method survey::deff
deff.drawn_mean <- function(object, quietly = FALSE, ...) object$deff

#' Shared argument handling for the estimators
#' @noRd
ht_prepare <- function(sample, y, level, df = NULL) {
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

  yv <- if (is.character(y) && length(y) == 1L) {
    if (!y %in% names(sample)) {
      stop("`sample` has no column `", y, "`.", call. = FALSE)
    }
    sample[[y]]
  } else {
    y
  }
  y_logical <- is.logical(yv)
  if (y_logical) yv <- as.numeric(yv)
  if (!is.numeric(yv) || length(yv) != nrow(sample)) {
    stop("`y` must be a numeric column name or a numeric vector as long as ",
         "`sample`.", call. = FALSE)
  }
  missing_y <- is.na(yv) & !is.nan(yv)
  if (any(missing_y)) {
    stop(sum(missing_y), " value(s) of `y` are missing.", call. = FALSE)
  }
  if (!all(is.finite(yv))) {
    bad <- which(!is.finite(yv))
    stop("`y` is ", paste(unique(as.character(yv[bad])), collapse = " or "),
         " in ", length(bad), " row(s) of `sample`: ", format_rows(bad),
         ". A total or mean needs a finite `y` for every sampled row.",
         call. = FALSE)
  }
  # A sample saved by drawn 0.2.0 records no probabilities for drawn_state()
  # to compare `.prob` with, so only its type and range are checked.
  if (is.null(state$pi_all)) check_prob(sample)
  warn_weight_edited(sample)
  # The rows the design can reach: those with a probability above 0 at draw
  # time, or none recorded because the sample is conditional on clusters
  # that left them out. A sample saved by drawn 0.2.0 has no such record, so
  # its count is read from the design's columns and the unreachable rows are
  # unknown.
  pi_all <- state$pi_all
  list(design = state$design, rows = state$rows, pop = state$pop,
       units = state$units, pi_all = pi_all, y = yv, pi = sample$.prob,
       sample = sample,
       y_name = if (is.character(y) && length(y) == 1L) y,
       y_logical = y_logical,
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
#' Returns the interval, its `scale`, `edge` (the beta interval was used for
#' an estimate of 0 or 1) and `exact`.
#'
#' @noRd
ci_for <- function(est, v, level, df, alternative = "two.sided",
                   proportion = FALSE, ci_method = "logit", pi_i = NULL) {
  if (!proportion) {
    return(list(ci = interval(est, v, level, df, alternative = alternative),
                scale = "linear", edge = FALSE, exact = FALSE))
  }
  n <- sum(!is_certain(pi_i))
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
#' design effect: see `deff_value()`.
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

  out <- list(variance = v, se = se, ci = ci$ci, level = level,
              alternative = alternative, df = df, n = n,
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
  print_estimate(x, "mean", if (identical(x$estimator, "ht"))
    "Horvitz-Thompson mean" else "Hajek mean")
}

#' @noRd
df_label <- function(df) {
  if (is.null(df) || !is.finite(df)) "normal" else
    paste0("t, ", format(signif(df, 4)), " df")
}

#' @noRd
method_label <- function(m) {
  switch(m %||% "analytic",
    deville = "Deville approximation",
    `successive difference` = "successive-difference approximation",
    `local mean` = "local-mean approximation",
    m
  )
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
