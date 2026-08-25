#' Estimate a population mean
#'
#' The Hajek estimator, `sum(y / pi) / sum(1 / pi)`, is the default. It divides
#' by the *estimated* population size rather than the known one, so a sample
#' that happens to over-represent heavy-weight rows inflates numerator and
#' denominator together and they partly cancel. `estimator = "ht"` divides by
#' the true `N` instead.
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
#' Hajek is usually the steadier of the two and is the default for that reason.
#' The exception is worth knowing: when `y` is close to proportional to the size
#' measure that drove selection, `y / pi` is nearly constant, the
#' Horvitz-Thompson numerator barely moves, and dividing it by the known `N`
#' beats dividing by an estimate.
#'
#' One caveat on `"ht"`. It is unbiased for the frame mean *provided every row
#' could have been selected*. Rows with inclusion probability 0 — outside a time
#' window, outside a region, zero weight — sit inside the `N` it divides by but
#' can never enter the numerator, so the estimate is biased low by exactly their
#' share of the frame. [sample_summary()] reports how many such rows there are.
#' The Hajek mean is unaffected, because it estimates the mean of the part of
#' the frame the design can actually reach.
#'
#' @param sample A data frame returned by [draw()] with `weights = TRUE`.
#' @param y The variable to average: a column name, or a numeric vector as long
#'   as `sample`.
#' @param estimator `"hajek"` or `"ht"`. See above.
#' @param variance Passed to the underlying total. See [ht_total()].
#' @param level Confidence level for the interval.
#'
#' @return A list with a `print()` method, holding `mean`, `variance`, `se`,
#'   `ci`, `level`, `n`, `design`, `deff`, `method` and `note`. See [ht_total()]
#'   for what `method` and `note` say, and [deff()] for the design effect.
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
#' # Stratified without replacement: the weights sum to N, so the two agree
#' ht_mean(s, "spend", estimator = "ht")$mean == ht_mean(s, "spend")$mean
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
                    level = 0.95) {
  estimator <- match.arg(estimator)
  variance <- match.arg(variance)
  parts <- ht_prepare(sample, y, level)

  pi_i <- parts$pi
  yv <- parts$y
  n_hat <- sum(1 / pi_i)
  N <- nrow(parts$pop)

  denom <- if (estimator == "hajek") n_hat else N
  est <- sum(yv / pi_i) / denom

  # Linearise: the variance of a ratio is the variance of a total of residuals.
  # For the HT mean there is no denominator to vary, so the residual is just y.
  z <- if (estimator == "hajek") (yv - est) / denom else yv / N

  var_out <- ht_variance_dispatch(parts$design, parts$sample, parts$pop,
                                  parts$rows, z, pi_i, variance)
  finish_estimate(est, var_out, level, nrow(sample), parts$design, yv, pi_i, N,
                  what = "mean", class = "drawn_mean")
}

#' Design effect
#'
#' How much precision the design costs against simple random sampling of the
#' same size. `deff = 1` means the design is doing as well as a coin flip over
#' the frame; above 1 it is doing worse, which is the usual price of clustering;
#' below 1 it is doing better, which is what stratification and
#' probability-proportional-to-size buy you.
#'
#' Read it as an exchange rate on sample size: at `deff = 2`, a sample of 400
#' carries about as much information as 200 drawn at random.
#'
#' @param x A result from [ht_total()] or [ht_mean()].
#'
#' @return A single number, or `NA` when the design has no variance estimate.
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
#' @seealso [ht_total()], [ht_mean()]
#' @export
deff <- function(x) {
  if (!inherits(x, c("drawn_ht", "drawn_mean"))) {
    stop("`x` must come from ht_total() or ht_mean().", call. = FALSE)
  }
  x$deff
}

#' Shared argument handling for the estimators
#' @noRd
ht_prepare <- function(sample, y, level) {
  if (!is.data.frame(sample)) {
    stop("`sample` must be a data frame returned by draw().", call. = FALSE)
  }
  design <- attr(sample, "drawn_design")
  rows <- attr(sample, "drawn_rows")
  pop <- attr(sample, "drawn_population")
  if (is.null(design) || is.null(rows) || is.null(pop)) {
    stop("`sample` does not carry its design. Draw it with ",
         "draw(..., weights = TRUE), and estimate from that result before ",
         "subsetting it.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }

  yv <- if (is.character(y) && length(y) == 1L) {
    if (!y %in% names(sample)) {
      stop("`sample` has no column `", y, "`.", call. = FALSE)
    }
    sample[[y]]
  } else {
    y
  }
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

#' Route to the requested variance estimator, falling back where allowed
#'
#' `sample` is passed through whole rather than reduced to its probabilities:
#' the jackknife reads the clustering column off it to find the primary
#' sampling units, and without that column it would delete one row at a time
#' and overstate the variance.
#'
#' @noRd
ht_variance_dispatch <- function(design, sample, pop, rows, z, pi_i, variance) {
  if (variance == "none") {
    return(list(variance = NA_real_, method = "none",
                note = "Variance not requested."))
  }
  if (variance == "jackknife") {
    return(jackknife_variance(design, sample, pop, z, pi_i))
  }

  # A Sen-Yates-Grundy estimate can come out negative, which is a failure of the
  # estimator rather than a variance -- treat it as one instead of passing a
  # number downstream that sqrt() and deff() then have to guess about.
  reason <- NULL
  got <- tryCatch({
    a <- c(ht_variance(design, pop, rows, z, pi_i), list(method = "analytic"))
    if (!is.na(a$variance) && a$variance < 0) {
      reason <- paste0("The Sen-Yates-Grundy estimator returned a negative ",
                       "variance (", signif(a$variance, 3), "), which it can ",
                       "do on an unlucky sample. There is no analytic figure ",
                       "to report.")
      NULL
    } else {
      a
    }
  }, error = function(e) {
    reason <<- conditionMessage(e)
    NULL
  })
  if (!is.null(got)) return(got)

  if (variance == "analytic") {
    return(list(variance = NA_real_, method = "analytic", note = reason))
  }

  jk <- jackknife_variance(design, sample, pop, z, pi_i)
  first <- sub("\n.*", "", reason)
  if (is.na(jk$variance)) {
    # The jackknife declined too. Say so, and keep its reason rather than
    # claiming a fallback that did not happen.
    return(list(variance = NA_real_, method = "none",
                note = paste0("No variance is available for this design. ",
                              first, " ", jk$note)))
  }
  jk$note <- paste0("No analytic variance for this design, so the jackknife ",
                    "was used instead. ", first, " ", jk$note)
  jk
}

#' Assemble an estimate, its interval, and its design effect
#' @noRd
finish_estimate <- function(est, var_out, level, n, design, y, pi_i, N,
                            what, class) {
  v <- var_out$variance
  se <- if (is.na(v) || v < 0) NA_real_ else sqrt(v)
  z <- stats::qnorm(1 - (1 - level) / 2)
  ci <- if (is.na(se)) c(NA_real_, NA_real_) else est + c(-1, 1) * z * se

  out <- list(variance = v, se = se, ci = ci, level = level, n = n,
              design = design_type(design),
              deff = deff_value(v, y, pi_i, N, n, what),
              method = var_out$method %||% "analytic", note = var_out$note)
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
print.drawn_mean <- function(x, ...) print_estimate(x, "mean", "Hajek mean")

#' @noRd
print_estimate <- function(x, field, label) {
  cat(label, "  (", x$design, " design, n = ", x$n, ")\n", sep = "")
  cat("  estimate ", format(x[[field]], big.mark = ","), "\n", sep = "")
  if (is.na(x$se)) {
    cat("  se       NA\n")
    if (!is.null(x$note)) {
      cat("\n", strwrap(x$note, prefix = "  "), sep = "\n"); cat("\n")
    }
  } else {
    cat("  se       ", format(x$se, big.mark = ","), "  (", x$method, ")\n",
        sep = "")
    cat("  ", format(100 * x$level), "% CI  ",
        format(x$ci[1], big.mark = ","), " to ",
        format(x$ci[2], big.mark = ","), "\n", sep = "")
    if (is.finite(x$deff)) {
      cat("  deff     ", signif(x$deff, 3), "  (",
          if (x$deff > 1.05) "worse than" else if (x$deff < 0.95) "better than"
          else "about the same as", " simple random sampling)\n", sep = "")
    }
  }
  invisible(x)
}
