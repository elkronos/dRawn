#' Estimate population quantiles from a sample
#'
#' `ht_quantile()` estimates quantiles of `y`, such as the median or the
#' quartiles, from the Hajek estimate of its distribution function,
#' `F(t) = sum(w * (y <= t)) / sum(w)` with `w = 1 / pi`, and gives each a
#' Woodruff confidence interval and a standard error from the design's own
#' variance estimator.
#'
#' @section The quantile rule:
#' The estimate at `p` is the smallest sampled value of `y` whose estimated
#' distribution function is `p` or more, which inverts the step function
#' `F` without interpolating between sampled values. It is always a value
#' some sampled row has. This is the rule of
#' `survey::svyquantile(qrule = "math")`, survey's default, and with equal
#' weights it is `stats::quantile(type = 1)`. `p = 0` gives the smallest
#' sampled value and `p = 1` the largest.
#'
#' @section Interval and standard error:
#' The interval is Woodruff's (1952). At the estimate `q`, `F(q)` is a Hajek
#' mean of the indicator `y <= q`, and its variance is the variance of the
#' total of the linearised residuals `((y <= q) - F(q)) / N_hat`, with
#' `N_hat = sum(w)`, from the design's own estimator, as for [ht_mean()].
#' "Variance" and "Confidence intervals" in [ht_total()] describe the
#' estimators, the fallback to the jackknife, the degrees of freedom and the
#' notes. The interval `F(q)` plus or minus `t` standard errors, with `t` on
#' the design's degrees of freedom, is carried back through the quantile
#' rule, so its ends are the estimated quantiles at its two ends. The standard
#' error is the interval's width divided by `2 t`, and `variance` is its
#' square. `cdf` and `cdf_se` in the result hold `F(q)` and its standard
#' error.
#'
#' Where the interval for `F(q)` runs below 0 or above 1, the sample bounds
#' the quantile on one side only. That end of `ci` is then `NA`, and so are
#' `se` and `variance`, with a note.
#'
#' At the largest sampled value `F(q)` is 1 and its residuals are 0, so its
#' estimated variance is 0 whatever the spread of `y`, and the interval
#' would have no width. A quantile there, which a `p` close to 1 or a `y`
#' with one value gives, gets no standard error, with a note, unless every
#' frame row the design can reach was taken with certainty, when it is
#' exact. A quantile whose sampled rows lie in one primary sampling unit
#' gets none either, as a Hajek mean does, because its residuals sum to
#' zero over that unit.
#'
#' A quantile has no design effect, and `deff` is `NA`.
#'
#' @section Domains:
#' `by` gives the quantiles of each domain, each from the distribution
#' function of the domain's sampled rows. The variance of `F(q)` is that of
#' the domain's Hajek mean of `y <= q`, from the whole sample, as described
#' under "Domains" in [ht_total()], and the degrees of freedom are the
#' design's. The result is a `drawn_by` table with a row for each domain
#' and probability.
#'
#' @section Comparison with survey:
#' For the designs [as_svydesign()] maps with identical standard errors,
#' the estimate, the interval and the standard error equal those of
#' `survey::svyquantile(ci = TRUE, interval.type = "mean")` on its result,
#' whose interval is also centred on `F(q)` rather than on `p`. They differ
#' in these places:
#'
#' * At the largest sampled value survey gives an interval of no width and a
#'   standard error of 0, where `ht_quantile()` gives none.
#' * Where the interval for `F(q)` passes 0 or 1, survey gives `NaN` for
#'   that end and the standard error, where `ht_quantile()` gives `NA`.
#' * For a `p`, or an end of the interval for `F(q)`, below the share of the
#'   weight held by the smallest sampled value, survey's rule gives the
#'   second smallest value and `ht_quantile()` the smallest.
#' * With `by`, `survey::svyby()` with `survey::svyquantile()` forms each
#'   interval on the degrees of freedom of the domain's own units, and
#'   `ht_quantile()` on the design's, as for every domain estimate here.
#'   With `df = Inf` in both, the standard errors agree.
#'
#' For systematic, spread and systematic PPS samples the variance of `F(q)`
#' comes from the estimator this package uses for their totals, which
#' `survey` does not have.
#'
#' A missing `y` is refused. `survey::svyquantile(na.rm = TRUE)` on
#' [as_svydesign()] estimates the quantiles over the rows that have one.
#'
#' @inheritParams ht_total
#' @param y The variable: a numeric column name of `sample`, or a numeric
#'   vector as long as `sample`.
#' @param probs The probabilities, each from 0 to 1, such as `0.5` for the
#'   median or `c(0.25, 0.5, 0.75)` for the quartiles.
#'
#' @return For one probability without `by`, a list of class
#'   `drawn_quantile` with a `print()` method, holding `quantile`, the
#'   estimate; `prob`; `variable`, the name of `y` when it is a column name,
#'   else `NULL`; `cdf` and `cdf_se`, the estimated distribution function
#'   at the estimate and its standard error; `steps`, the estimated
#'   distribution function as the sorted sampled values `x` and the share
#'   of the weight `F` at or below each, which `confint()` reads to form the
#'   interval at another level; and `variance`, `se`, `ci`, `level`,
#'   `alternative` (always `"two.sided"`), `df`, `n`, `n_random`, `design`,
#'   `deff`, `method`, `note` and `unreachable` as for [ht_total()], with
#'   `ci_scale` `"woodruff"`.
#'
#'   For several probabilities, or with `by`, a `drawn_by` table with a row
#'   for each domain and probability: the `by` columns, `prob`, `n`,
#'   `quantile`, `variance`, `se`, `ci_lower`, `ci_upper` and `method`.
#'
#'   Both have `coef()`, `vcov()`, `confint()`, `as.data.frame()` and
#'   `tidy()` methods, as described under "Working with the result" in
#'   [ht_total()]. `confint()` at another level carries the interval for
#'   `F(q)` at that level back through the quantile rule. `vcov()` is
#'   diagonal: the covariances between quantiles are not estimated.
#'
#' @references
#' Woodruff, R. S. (1952). Confidence intervals for medians and other
#' position measures. *Journal of the American Statistical Association*,
#' 47, 635–646.
#'
#' Särndal, C.-E., Swensson, B. and Wretman, J. (1992). *Model Assisted Survey
#' Sampling*. Springer.
#'
#' @examples
#' set.seed(1)
#' pop <- data.frame(
#'   id = 1:300,
#'   site = rep(c("a", "b", "c"), each = 100),
#'   spend = round(stats::rlnorm(300, 4, 0.8))
#' )
#' s <- draw(pop, design_stratified("site", n = 60), seed = 1, weights = TRUE)
#'
#' # The median, and its Woodruff interval
#' ht_quantile(s, "spend")
#' stats::median(pop$spend)   # the truth
#'
#' # The quartiles of each site, from the whole sample
#' ht_quantile(s, "spend", probs = c(0.25, 0.5, 0.75), by = "site")
#'
#' @seealso [ht_mean()], [ht_total()], [as_svydesign()]
#' @export
ht_quantile <- function(sample, y, probs = 0.5, by = NULL,
                        variance = c("auto", "analytic", "jackknife", "none"),
                        level = 0.95, df = NULL) {
  variance <- check_choice(variance,
                           c("auto", "analytic", "jackknife", "none"),
                           "variance")
  check_probs(probs)
  parts <- prepare_sample(sample, level, df)
  parts <- add_quantile_variable(parts, y)
  parts <- finish_prepare(parts)
  if (!nrow(sample)) {
    stop("`sample` has no rows, so there is no quantile to estimate.",
         call. = FALSE)
  }
  ctx <- variance_context()
  if (is.null(by)) warn_grouped(sample, "ht_quantile")
  if (is.null(by) && length(probs) == 1L) {
    return(quantile_result(parts, quantile_cells(parts, NULL, probs, variance,
                                                 level, df, ctx)[[1]],
                           level))
  }
  quantile_table(parts, probs, by, variance, level, df, ctx)
}

#' Check the probabilities of ht_quantile()
#' @noRd
check_probs <- function(probs) {
  if (!is.numeric(probs) || !length(probs) || anyNA(probs) ||
      any(probs < 0 | probs > 1)) {
    stop("`probs` must be one or more numbers from 0 to 1, such as 0.5 for ",
         "the median, not ", format_bad(probs), ".", call. = FALSE)
  }
  twice <- unique(probs[duplicated(probs)])
  if (length(twice)) {
    stop("`probs` holds ", paste(format(twice), collapse = ", "), " more ",
         "than once. Give each probability once.", call. = FALSE)
  }
  invisible(probs)
}

#' Read the `y` of a quantile into the parts of a sample
#'
#' It must be present and finite in every sampled row.
#'
#' @noRd
add_quantile_variable <- function(parts, y) {
  v <- read_variable(parts$sample, y, "y")
  gap <- is.na(v$values) & !is.nan(v$values)
  if (any(gap)) {
    stop("`y` is missing in ", sum(gap), " of the ", length(gap), " sampled ",
         "rows. A quantile needs `y` in every sampled row. ",
         "survey::svyquantile(na.rm = TRUE) on as_svydesign() estimates it ",
         "over the rows that have one.", call. = FALSE)
  }
  check_finite_values(v$values, gap, "`y`",
                      "A quantile needs a finite `y` for every sampled row.")
  parts$y <- v$values
  parts$respond <- rep(TRUE, length(v$values))
  parts$y_name <- v$name
  parts$y_logical <- v$logical
  parts
}

#' The estimated distribution function of `y` with weights `w`
#'
#' Returns the sorted distinct values `x` and `F`, the share of the weight
#' at or below each. The last share is exactly 1.
#'
#' @noRd
cdf_steps <- function(y, w) {
  o <- order(y, method = "radix")
  y <- y[o]
  cw <- cumsum(w[o])
  last <- c(y[-1L] != y[-length(y)], TRUE)
  list(x = y[last], F = cw[last] / cw[length(cw)])
}

#' The position in `steps` of the quantile at `p`
#'
#' The first value whose share is `p` or more. A share equal to `p` can come
#' out of the cumulative sum a rounding error below it, and the tolerance
#' keeps the quantile at that value rather than the next.
#'
#' @noRd
step_index <- function(steps, p) {
  match(TRUE, steps$F >= p - 1e-9)
}

#' The quantile at `p` of an estimated distribution function
#' @noRd
step_quantile <- function(steps, p) steps$x[step_index(steps, p)]

#' The Woodruff interval of a quantile
#'
#' `q` is the estimate, `cdf` the distribution function there and `v` its
#' variance. The interval for `F(q)` at `level` on `df` degrees of freedom
#' is carried back through the quantile rule. An end past 0 or 1 is `NA`.
#' Returns `ci`, `se`, the interval's width over `2 t`, and `open`, whether
#' an end is missing for that reason.
#'
#' @noRd
woodruff <- function(steps, q, cdf, v, level, df) {
  none <- list(ci = c(NA_real_, NA_real_), se = NA_real_, open = FALSE)
  if (!is.finite(v) || v < 0) return(none)
  if (v == 0) return(list(ci = c(q, q), se = 0, open = FALSE))
  if (is.finite(df) && df < 1) return(none)
  a <- (1 - level) / 2
  t <- if (is.finite(df)) stats::qt(1 - a, df) else stats::qnorm(1 - a)
  ends <- cdf + c(-1, 1) * t * sqrt(v)
  ci <- c(if (ends[1] < 0) NA_real_ else step_quantile(steps, ends[1]),
          if (ends[2] > 1) NA_real_ else step_quantile(steps, ends[2]))
  list(ci = ci, se = if (anyNA(ci)) NA_real_ else (ci[2] - ci[1]) / (2 * t),
       open = anyNA(ci))
}

#' The quantiles of each domain at each probability, with their intervals
#'
#' `dom` is the domains of `by` from domain_keys(), or `NULL` for the whole
#' sample. For each probability the residuals of every domain go into one
#' variable, each 0 outside its domain, and the variance of each domain's
#' `F(q)` comes from domain_variances(), as for a domain's Hajek mean, or
#' for the whole sample from the design's own estimator. Returns a list with
#' an element for each domain and probability, domain by domain, holding
#' the domain's position `j`, `prob`, `n`, `n_random`, `steps`, the
#' estimate `quantile`, `cdf`, `cdf_se`, `ci`, `se`, `method`, `note` and
#' the flags `top` (a quantile at the largest value, whose variance is
#' declined), `open` and `no_df` that the notes are built from.
#'
#' @noRd
quantile_cells <- function(parts, dom, probs, variance, level, df, ctx) {
  n <- nrow(parts$sample)
  pi_i <- parts$pi
  yv <- parts$y
  groups <- if (is.null(dom)) list(seq_len(n)) else dom$idx
  k <- length(groups)
  steps <- lapply(groups, function(i) cdf_steps(yv[i], 1 / pi_i[i]))
  n_hat <- vapply(groups, function(i) sum(1 / pi_i[i]), numeric(1))
  df_use <- df %||% design_df(parts$design, parts$units, parts$pop,
                              parts$rows, pi_i)
  psu <- if (is.null(dom) && variance != "none") {
    sample_units(parts$design, parts$units, parts$pop, parts$rows)
  }
  cells <- vector("list", k * length(probs))
  for (h in seq_along(probs)) {
    fits <- lapply(steps, function(st) {
      at <- step_index(st, probs[h])
      list(q = st$x[at], cdf = st$F[at])
    })
    z <- numeric(n)
    for (j in seq_len(k)) {
      i <- groups[[j]]
      z[i] <- (as.numeric(yv[i] <= fits[[j]]$q) - fits[[j]]$cdf) / n_hat[j]
    }
    vars <- if (is.null(dom)) {
      declined <- if (variance != "none") {
        ratio_decline(parts, rep(TRUE, n), psu = psu, estimate = "quantile")
      }
      list(declined %||%
             ht_variance_dispatch(parts$design, parts$units, parts$pop,
                                  parts$rows, z, pi_i, variance, ctx = ctx))
    } else {
      domain_variances(parts, dom, list(z = z), variance, "quantile", ctx)
    }
    for (j in seq_len(k)) {
      fit <- fits[[j]]
      v <- vars[[j]]
      frame <- if (is.null(dom)) TRUE else domain_frame(dom, j)
      # At the largest value the residuals are 0 and so is the variance,
      # however `y` is spread, unless the domain is a census.
      top <- variance != "none" && fit$cdf >= 1 &&
        is.finite(v$variance) && !certain_domain(parts, frame)
      if (top) v <- list(variance = NA_real_, method = "none", note = NULL)
      w <- woodruff(steps[[j]], fit$q, fit$cdf, v$variance, level, df_use)
      i <- groups[[j]]
      cells[[(j - 1L) * length(probs) + h]] <- list(
        j = j, prob = probs[h], n = length(i),
        n_random = sum(!is_certain(pi_i[i])), steps = steps[[j]],
        quantile = fit$q, cdf = fit$cdf,
        cdf_se = if (is.finite(v$variance) && v$variance >= 0) {
          sqrt(v$variance)
        } else {
          NA_real_
        },
        ci = w$ci, se = w$se, method = v$method %||% "analytic",
        note = v$note, top = top, open = w$open,
        no_df = is.finite(v$variance) && v$variance > 0 &&
          is.finite(df_use) && df_use < 1,
        df = df_use)
    }
  }
  cells
}

#' One quantile without `by`, as a `drawn_quantile`
#' @noRd
quantile_result <- function(parts, cell, level) {
  notes <- c(cell$note,
             if (cell$top) top_note(),
             if (cell$open) open_note(),
             if (cell$no_df) df_note(cell$df))
  structure(list(
    quantile = cell$quantile, prob = cell$prob, variable = parts$y_name,
    cdf = cell$cdf, cdf_se = cell$cdf_se, variance = cell$se^2,
    se = cell$se, ci = cell$ci, level = level, alternative = "two.sided",
    df = cell$df, n = cell$n, n_random = cell$n_random,
    ci_scale = "woodruff", design = design_type(parts$design),
    deff = NA_real_, method = cell$method,
    note = if (length(notes)) paste(notes, collapse = " "),
    unreachable = parts$unreachable, steps = cell$steps
  ), class = "drawn_quantile")
}

#' Quantiles at several probabilities, or by domain, as a `drawn_by` table
#'
#' A row for each domain and probability, domain by domain. What `confint()`
#' needs to form an interval at another level, each row's distribution
#' function and the standard error of `F(q)`, is kept in the `woodruff`
#' attribute by the row's name, so that it still matches after rows are
#' selected.
#'
#' @noRd
quantile_table <- function(parts, probs, by, variance, level, df, ctx) {
  dom <- NULL
  if (!is.null(by)) {
    dom <- domain_keys(parts, by)
    check_table_names(by, c("prob", "n", "quantile", "variance", "se",
                            "ci_lower", "ci_upper", "method"))
  }
  cells <- quantile_cells(parts, dom, probs, variance, level, df, ctx)
  field <- function(name, type) vapply(cells, `[[`, type, name)
  j <- field("j", integer(1))
  out <- if (is.null(dom)) {
    data.frame(prob = field("prob", numeric(1)))
  } else {
    first <- match(dom$levels_in, dom$key)[j]
    cols <- data.frame(stats::setNames(lapply(by, function(b) {
      parts$sample[[b]][first]
    }), by), check.names = FALSE)
    cols$prob <- field("prob", numeric(1))
    cols
  }
  out$n <- field("n", integer(1))
  out$quantile <- field("quantile", numeric(1))
  se <- field("se", numeric(1))
  out$variance <- se^2
  out$se <- se
  ci <- vapply(cells, `[[`, numeric(2), "ci")
  out$ci_lower <- ci[1, ]
  out$ci_upper <- ci[2, ]
  out$method <- field("method", character(1))
  out <- structure(out, class = c("drawn_by", "data.frame"),
                   what = "quantile", level = level,
                   alternative = "two.sided", df = cells[[1]]$df,
                   design = design_type(parts$design), ci_scale = "woodruff",
                   unreachable = parts$unreachable, variable = parts$y_name,
                   by = c(by, "prob"))
  nm <- domain_names(out)
  flag <- function(name) nm[field(name, logical(1))]
  notes <- unique(unlist(lapply(cells, `[[`, "note")))
  if (length(flag("top"))) notes <- c(notes, top_note(flag("top")))
  if (length(flag("open"))) notes <- c(notes, open_note(flag("open")))
  if (length(flag("no_df"))) notes <- c(notes, df_note(cells[[1]]$df))
  attr(out, "notes") <- notes
  attr(out, "n_random") <- stats::setNames(field("n_random", integer(1)), nm)
  attr(out, "woodruff") <- stats::setNames(lapply(cells, function(cell) {
    cell[c("steps", "quantile", "cdf", "cdf_se")]
  }), nm)
  out
}

#' Why a quantile at the largest sampled value has no standard error
#'
#' `labels` names the rows of a table, as domain_names() names them.
#'
#' @noRd
top_note <- function(labels = NULL) {
  core <- paste0("where the estimated distribution function reaches 1 and ",
                 "has an estimated variance of 0 whatever the spread of ",
                 "`y`, so the Woodruff interval would have no width. A ",
                 "standard error needs a quantile below the largest sampled ",
                 "value: a smaller `probs`, or a larger sample.")
  if (is.null(labels)) {
    return(paste0("The estimate is the largest sampled value of `y`, ", core))
  }
  several <- length(labels) > 1L
  paste0("The estimate", if (several) "s", " for ", quote_labels(labels),
         if (several) " are" else " is", " the largest sampled value of `y` ",
         "in ", if (several) "their domains" else "its domain", ", ", core)
}

#' Why an end of a quantile's interval is missing
#' @noRd
open_note <- function(labels = NULL) {
  what <- if (is.null(labels)) {
    "the estimate"
  } else {
    paste0("the estimate", if (length(labels) > 1L) "s", " for ",
           quote_labels(labels))
  }
  paste0("The interval for the distribution function at ", what, ", `F(q)` ",
         "plus or minus `t` standard errors, runs below 0 or above 1, so the ",
         "sample bounds the quantile on one side only. That end of the ",
         "interval is NA, and so is the standard error, which is the ",
         "interval's width over `2 t`. A larger sample, or a probability ",
         "further from 0 and 1, narrows it.")
}

#' The Woodruff intervals of a quantile table at `level`, a row per row
#' @noRd
quantile_by_ci <- function(x, level) {
  w <- attr(x, "woodruff")
  nm <- domain_names(x)
  if (is.null(w) || !all(nm %in% names(w))) {
    stop("`object` does not hold the distribution functions its intervals ",
         "need at another level. Estimate again with `level = ",
         format(level), "`.", call. = FALSE)
  }
  out <- vapply(w[nm], function(r) {
    woodruff(r$steps, r$quantile, r$cdf, r$cdf_se^2, level,
             attr(x, "df"))$ci
  }, numeric(2))
  t(out)
}

#' @export
print.drawn_quantile <- function(x, ...) {
  print_estimate(x, "quantile",
                 paste("Quantile of", ratio_name(x$variable, "y")),
                 lines = c(prob = format(x$prob)))
}

#' @export
#' @exportS3Method survey::deff
deff.drawn_quantile <- function(object, quietly = FALSE, ...) object$deff

#' @export
coef.drawn_quantile <- function(object, ...) coef.drawn_ht(object)

#' @export
vcov.drawn_quantile <- function(object, ...) vcov.drawn_ht(object)

#' @export
confint.drawn_quantile <- function(object, parm, level = object$level, ...) {
  confint.drawn_ht(object, parm, level)
}

#' @export
as.data.frame.drawn_quantile <- function(x, row.names = NULL,
                                         optional = FALSE, ...) {
  single_frame(x, row.names)
}

#' @exportS3Method generics::tidy
tidy.drawn_quantile <- function(x, conf.int = FALSE, conf.level = x$level,
                                ...) {
  tidy.drawn_ht(x, conf.int, conf.level)
}
