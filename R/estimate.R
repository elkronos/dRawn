#' Estimate a population mean
#'
#' The Hajek estimator, `sum(y / pi) / sum(1 / pi)`, is the default. It divides
#' by the *estimated* population size rather than the known one, so a sample
#' that happens to over-represent heavy-weight rows inflates numerator and
#' denominator together and they partly cancel. `estimator = "ht"` divides by
#' the true `N` instead. For a 0/1 variable either one estimates a proportion.
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
#'   [design_systematic()] where the interval does not divide `N`;
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
#' window, outside a region, zero weight — sit inside the `N` it divides by but
#' can never enter the numerator, so the estimate is biased low by exactly their
#' share of the frame. [sample_summary()] reports how many such rows there are.
#' The Hajek mean is unaffected, because it estimates the mean of the part of
#' the frame the design can actually reach.
#'
#' @section Variance and domains:
#' The Hajek mean is a ratio, and its variance is the variance of the total of
#' the linearised residuals `(y - mean) / N_hat` (Deville 1999), computed with
#' whichever estimator suits the design — see [ht_total()], which also explains
#' the confidence interval's degrees of freedom and how `by` estimates domain
#' means from the whole sample rather than a subset of it.
#'
#' @param sample A data frame returned by [draw()] with `weights = TRUE`.
#' @param y The variable to average: a column name, or a numeric vector as long
#'   as `sample`.
#' @param estimator `"hajek"` or `"ht"`. See above.
#' @param variance,level,by,df As for [ht_total()].
#'
#' @return A list with a `print()` method, holding `mean`, `estimator`,
#'   `variance`, `se`, `ci`, `level`, `df`, `n`, `design`, `deff`, `method` and
#'   `note`. See [ht_total()] for what `method`, `df` and `note` say, and
#'   [deff()] for the design effect. With `by`, a data frame with one row per
#'   domain.
#'
#' @references
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
#' # A proportion is the mean of a 0/1 variable
#' ht_mean(s, s$spend > 250)
#'
#' # Stratified without replacement: the weights sum to N, so the two agree
#' all.equal(ht_mean(s, "spend", estimator = "ht")$mean, ht_mean(s, "spend")$mean)
#'
#' # Poisson sampling has a random size, so they part company
#' p <- draw(pop, design_weighted("spend", n = 40, method = "poisson"),
#'           seed = 1, weights = TRUE)
#' c(hajek = ht_mean(p, "spend", variance = "none")$mean,
#'   ht    = ht_mean(p, "spend", "ht", variance = "none")$mean)
#'
#' @seealso [ht_total()], [deff()], [sample_summary()]
#' @export
ht_mean <- function(sample, y, estimator = c("hajek", "ht"),
                    variance = c("auto", "analytic", "jackknife", "none"),
                    level = 0.95, by = NULL, df = NULL) {
  estimator <- match.arg(estimator)
  variance <- match.arg(variance)
  parts <- ht_prepare(sample, y, level, df)
  if (!is.null(by)) {
    return(estimate_by(parts, by, variance, level, df, what = "mean",
                       estimator = estimator))
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

  estimate_one(parts, est, z, variance, level, df, what = "mean",
               class = "drawn_mean", deff_y = yv,
               extra = list(estimator = estimator))
}

#' Estimates for each domain, from the whole sample
#'
#' A domain estimate uses the full sample with `y` zeroed outside the domain.
#' Subsetting first would treat the number of rows that happened to land in the
#' domain as fixed by design, which it is not, and would hand the variance
#' estimator a sample the design never drew.
#'
#' @noRd
estimate_by <- function(parts, by, variance, level, df, what,
                        estimator = "hajek") {
  if (!is.character(by) || !length(by) || anyNA(by) ||
      !all(by %in% names(parts$sample))) {
    stop("`by` must name one or more columns of `sample`.", call. = FALSE)
  }
  check_key_columns(parts$sample, by, "by")
  if (any(vapply(parts$sample[by], anyNA, logical(1)))) {
    stop("`by` columns must not contain missing values.", call. = FALSE)
  }
  key <- as.character(group_key(parts$sample, by))
  pop_key <- if (all(by %in% names(parts$pop))) {
    as.character(group_key(parts$pop, by))
  } else {
    NULL
  }
  levels_in <- sort(unique(key))
  pi_i <- parts$pi
  yv <- parts$y
  df_use <- df %||% design_df(parts$design, parts$sample, parts$pop,
                              parts$rows)

  rows <- lapply(levels_in, function(d) {
    ind <- as.numeric(key == d)
    if (what == "total") {
      est <- sum(ind * yv / pi_i)
      z <- ind * yv
    } else {
      denom <- if (estimator == "hajek") {
        sum(ind / pi_i)
      } else {
        if (is.null(pop_key)) {
          stop("estimator = \"ht\" needs the domain's population size, and ",
               "the `by` columns are not in the population.", call. = FALSE)
        }
        sum(pop_key == d, na.rm = TRUE)
      }
      est <- sum(ind * yv / pi_i) / denom
      z <- if (estimator == "hajek") ind * (yv - est) / denom else ind * yv / denom
    }
    v <- ht_variance_dispatch(parts$design, parts$sample, parts$pop,
                              parts$rows, z, pi_i, variance)
    ci <- interval(est, v$variance, level, df_use)
    list(n = sum(ind), est = est, se = if (is.finite(v$variance) &&
                                             v$variance >= 0)
           sqrt(v$variance) else NA_real_,
         lo = ci[1], hi = ci[2], method = v$method %||% "analytic",
         note = v$note)
  })

  first <- match(levels_in, key)
  out <- parts$sample[first, by, drop = FALSE]
  out <- as.data.frame(out)
  rownames(out) <- NULL
  out$n <- vapply(rows, `[[`, numeric(1), "n")
  out[[what]] <- vapply(rows, `[[`, numeric(1), "est")
  out$se <- vapply(rows, `[[`, numeric(1), "se")
  out$ci_lower <- vapply(rows, `[[`, numeric(1), "lo")
  out$ci_upper <- vapply(rows, `[[`, numeric(1), "hi")
  out$method <- vapply(rows, `[[`, character(1), "method")
  notes <- unique(unlist(lapply(rows, `[[`, "note")))
  structure(out, class = c("drawn_by", "data.frame"),
            what = what, estimator = if (what == "mean") estimator,
            level = level, df = df_use, design = design_type(parts$design),
            notes = notes)
}

#' @export
print.drawn_by <- function(x, ...) {
  what <- attr(x, "what")
  label <- if (what == "total") "Horvitz-Thompson total" else
    if (identical(attr(x, "estimator"), "ht")) "Horvitz-Thompson mean" else
      "Hajek mean"
  cat(label, " by domain  (", attr(x, "design"), " design, ",
      format(100 * attr(x, "level")), "% CI, ", df_label(attr(x, "df")),
      ")\n", sep = "")
  y <- x
  attributes(y) <- attributes(x)[c("names", "row.names")]
  class(y) <- "data.frame"
  print(y, digits = 4, row.names = FALSE)
  notes <- attr(x, "notes")
  notes <- notes[nzchar(notes) & !grepl("^Variance not requested", notes)]
  if (length(notes)) {
    cat("\n", strwrap(notes, prefix = "  "), sep = "\n")
    cat("\n")
  }
  invisible(x)
}

#' Design effect
#'
#' How much precision the design costs against simple random sampling of the
#' same size (Kish 1965). `deff = 1` means the design is doing as well as a
#' coin flip over the frame; above 1 it is doing worse, which is the usual price
#' of clustering; below 1 it is doing better, which is what stratification and
#' probability-proportional-to-size buy you.
#'
#' Read it as an exchange rate on sample size: at `deff = 2`, a sample of 400
#' carries about as much information as 200 drawn at random.
#'
#' @param x A result from [ht_total()] or [ht_mean()].
#'
#' @return A single number, or `NA` when the design has no variance estimate.
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

#' Shared argument handling for the estimators
#' @noRd
ht_prepare <- function(sample, y, level, df = NULL) {
  if (!is.data.frame(sample)) {
    stop("`sample` must be a data frame returned by draw().", call. = FALSE)
  }
  design <- attr(sample, "drawn_design")
  rows <- attr(sample, "drawn_rows")
  pop <- attr(sample, "drawn_population")
  if (is.null(design) || is.null(rows) || is.null(pop)) {
    stop("`sample` does not carry its design. Draw it with ",
         "draw(..., weights = TRUE), and estimate from that result before ",
         "subsetting it -- use `by` for estimates within groups.",
         call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  if (!is.null(df) && (!is.numeric(df) || length(df) != 1L || is.na(df) ||
                       df <= 0)) {
    stop("`df` must be a single positive number, Inf, or NULL.", call. = FALSE)
  }

  yv <- if (is.character(y) && length(y) == 1L) {
    if (!y %in% names(sample)) {
      stop("`sample` has no column `", y, "`.", call. = FALSE)
    }
    sample[[y]]
  } else {
    y
  }
  if (is.logical(yv)) yv <- as.numeric(yv)
  if (!is.numeric(yv) || length(yv) != nrow(sample)) {
    stop("`y` must be a numeric column name or a numeric vector as long as ",
         "`sample`.", call. = FALSE)
  }
  if (anyNA(yv)) {
    stop(sum(is.na(yv)), " value(s) of `y` are missing.", call. = FALSE)
  }
  list(design = design, rows = rows, pop = pop, y = yv, pi = sample$.prob,
       sample = sample)
}

#' A confidence interval from an estimate, a variance and degrees of freedom
#' @noRd
interval <- function(est, v, level, df) {
  if (!is.finite(v) || v < 0) return(c(NA_real_, NA_real_))
  q <- if (is.finite(df)) {
    if (df < 1) return(c(NA_real_, NA_real_))
    stats::qt(1 - (1 - level) / 2, df)
  } else {
    stats::qnorm(1 - (1 - level) / 2)
  }
  est + c(-1, 1) * q * sqrt(v)
}

#' Assemble an estimate, its interval, and its design effect
#' @noRd
finish_estimate <- function(est, var_out, level, df, n, design, y, pi_i, N,
                            what, class) {
  v <- var_out$variance
  se <- if (is.na(v) || v < 0) NA_real_ else sqrt(v)
  ci <- interval(est, v, level, df)
  note <- var_out$note
  if (!is.na(se) && is.finite(df) && df < 1) {
    note <- c(note, paste0("There are no degrees of freedom left for an ",
                           "interval: the sample has one primary sampling ",
                           "unit per stratum."))
    note <- paste(note, collapse = " ")
  }

  out <- list(variance = v, se = se, ci = ci, level = level, df = df, n = n,
              design = design_type(design),
              deff = if (is.null(y)) NA_real_ else
                deff_value(v, y, pi_i, N, n, what),
              method = var_out$method %||% "analytic", note = note)
  out[[what]] <- est
  structure(out, class = class)
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
#' @noRd
deff_value <- function(v, y, pi_i, N, n, what) {
  if (!is.finite(v) || v <= 0 || n < 2L || !is.finite(N) || N < 2L) {
    return(NA_real_)
  }
  w <- 1 / pi_i
  mu <- sum(w * y) / sum(w)
  s2 <- sum(w * (y - mu)^2) * n / (sum(w) * (n - 1))
  if (!is.finite(s2) || s2 <= 0) return(NA_real_)
  fpc <- max(0, 1 - n / N)
  v_srs <- if (what == "mean") s2 / n * fpc else N^2 * s2 / n * fpc
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

#' @noRd
print_estimate <- function(x, field, label) {
  cat(label, "  (", x$design, " design, n = ", x$n, ")\n", sep = "")
  cat("  estimate ", format(x[[field]], big.mark = ","), "\n", sep = "")
  if (is.na(x$se)) {
    cat("  se       NA\n")
  } else {
    cat("  se       ", format(x$se, big.mark = ","), "  (",
        method_label(x$method), ")\n", sep = "")
    if (all(is.finite(x$ci))) {
      cat("  ", format(100 * x$level), "% CI  ",
          format(x$ci[1], big.mark = ","), " to ",
          format(x$ci[2], big.mark = ","), "  (", df_label(x$df), ")\n",
          sep = "")
    }
    if (is.finite(x$deff)) {
      cat("  deff     ", signif(x$deff, 3), "  (",
          if (x$deff > 1.05) "worse than" else if (x$deff < 0.95) "better than"
          else "about the same as", " simple random sampling)\n", sep = "")
    }
  }
  if (!is.null(x$note) && (is.na(x$se) || !identical(x$method, "analytic"))) {
    cat("\n", strwrap(x$note, prefix = "  "), sep = "\n")
    cat("\n")
  }
  invisible(x)
}
