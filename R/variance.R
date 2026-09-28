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
#'     nothing.}
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
#'     0.7 of the truth when the trend spanned 25 noise standard deviations).}
#'   \item{`"local mean"`}{[design_spread()]. A spatially balanced sample is
#'     well spread precisely because neighbours are rarely taken together, so
#'     this compares each sampled row with the local mean of its nearest
#'     sampled neighbours (Stevens and Olsen 2003; Grafström and Schelin 2014).}
#'   \item{`"jackknife"`}{The delete-one-unit jackknife, by stratum (JKn) where
#'     the design has strata, and deleting whole clusters where it has them.
#'     Used when requested, and as the fallback under `variance = "auto"`.}
#' }
#'
#' The analytic estimator declines — returning `NA` with the reason in `note`
#' — rather than understate. Two cases matter in practice:
#'
#' * **A stratum or interval with a single sampled row.** There is no
#'   variation inside it to measure, and simply leaving it out, as the
#'   Sen-Yates-Grundy sum does, reports a standard error that is too small —
#'   zero, if every stratum is like that. Draw at least two per stratum
#'   (`min_per_stratum = 2`, `per_interval = 2`).
#' * **One row per cluster in a multistage design.** The within-cluster
#'   variation is unmeasurable and the exact estimator understates by around
#'   a third. `variance = "auto"` falls back to the delete-a-cluster jackknife,
#'   which is the ultimate-cluster approximation standard in survey practice.
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
#' @section Domains:
#' `by` estimates for each level of one or more columns — sites, regions,
#' months — without subsetting the sample first. Subsetting a sample and
#' treating the piece as its own sample gets the variance wrong, because the
#' number of rows that happen to fall in a domain is itself random. Here each
#' domain's total is estimated from `y * (row in domain)` over the whole sample,
#' so the same design and the same variance machinery apply (Särndal, Swensson
#' and Wretman 1992, section 10.3).
#'
#' @param sample A data frame returned by [draw()] with `weights = TRUE`.
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
#'   (see "Confidence intervals"); `Inf` gives a normal interval.
#'
#' @return A list with a `print()` method, holding:
#'   \describe{
#'     \item{`total`}{The Horvitz-Thompson total, `sum(y / pi)`.}
#'     \item{`variance`, `se`, `ci`, `level`, `df`}{Its estimated variance,
#'       standard error, confidence interval, and the degrees of freedom the
#'       interval used. `NA` where the design supports no variance.}
#'     \item{`n`, `design`}{Rows used, and the design's type.}
#'     \item{`deff`}{The design effect — see [deff()].}
#'     \item{`method`}{Which estimator produced the variance. See "Variance".}
#'     \item{`note`}{Why a variance is missing, what approximation was used,
#'       or which fallback was taken. `NULL` when an exact estimator applied
#'       cleanly.}
#'   }
#'   With `by`, a data frame of class `drawn_by` instead: one row per domain,
#'   with the domain columns, `n`, `total`, `se`, `ci_lower`, `ci_upper` and
#'   `method`.
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
                     by = NULL, df = NULL) {
  variance <- match.arg(variance)
  parts <- ht_prepare(sample, y, level, df)
  if (!is.null(by)) {
    return(estimate_by(parts, by, variance, level, df, what = "total"))
  }
  yv <- parts$y
  pi_i <- parts$pi
  total <- sum(yv / pi_i)
  estimate_one(parts, total, z = yv, variance, level, df, what = "total",
               class = "drawn_ht", deff_y = yv)
}

#' Variance, interval and design effect for one estimate
#'
#' `z` is the variable whose Horvitz-Thompson total has the estimate's
#' variance: `y` itself for a total, the linearised residual for a ratio.
#'
#' @noRd
estimate_one <- function(parts, est, z, variance, level, df, what, class,
                         deff_y = NULL, extra = list(), logit = FALSE) {
  var_out <- ht_variance_dispatch(parts$design, parts$sample, parts$pop,
                                  parts$rows, z, parts$pi, variance)
  df <- df %||% design_df(parts$design, parts$sample, parts$pop, parts$rows)
  out <- finish_estimate(est, var_out, level, df, nrow(parts$sample),
                         parts$design, deff_y, parts$pi, nrow(parts$pop),
                         what = what, class = class, logit = logit)
  for (nm in names(extra)) out[[nm]] <- extra[[nm]]
  out
}

#' Route to the requested variance estimator, falling back where allowed
#'
#' `sample` is passed through whole rather than reduced to its probabilities:
#' the jackknife reads the clustering and stratification columns off it, and
#' without them it would delete one row at a time across strata and overstate
#' the variance many times over.
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

  # A negative Sen-Yates-Grundy estimate, or an estimator that declines with
  # NA, is a failure rather than a variance -- treat both as one instead of
  # passing a number downstream that sqrt() and deff() then have to guess about.
  reason <- NULL
  got <- tryCatch({
    a <- ht_variance(design, pop, rows, z, pi_i)
    a$method <- a$method %||% "analytic"
    if (is.na(a$variance)) {
      reason <- a$note %||% "The design's own estimator returned no figure."
      NULL
    } else if (a$variance < 0) {
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
    return(list(variance = NA_real_, method = "none", note = reason))
  }

  jk <- jackknife_variance(design, sample, pop, z, pi_i)
  first <- sub("\n.*", "", reason)
  if (is.na(jk$variance)) {
    # The jackknife declined too. Say so, and keep its reason rather than
    # claiming a fallback that did not happen.
    return(list(variance = NA_real_, method = "none",
                note = paste0("No variance is available for this sample. ",
                              first, " ", jk$note)))
  }
  jk$note <- paste0("The design's own estimator declined, so the jackknife ",
                    "was used instead. ", first, " ", jk$note)
  jk
}

#' Degrees of freedom for a design's confidence interval
#'
#' Primary sampling units minus strata, the usual rule (Korn and Graubard
#' 1999) and the one `survey::degf()` applies, so intervals agree across the two
#' packages.
#'
#' @noRd
design_df <- function(design, sample, pop, rows) {
  n <- nrow(sample)
  out <- switch(design_type(design),
    stratified = n - length(unique(group_key(sample, design$strata))),
    temporal = n - length(unique(temporal_bucket(design, sample))),
    cluster = ,
    multistage = length(unique(sample[[design$clusters]])) - 1,
    certainty = {
      sp <- certainty_split(design, pop)
      free <- !(rows %in% sp$keep[sp$take])
      if (!any(free)) {
        Inf
      } else {
        design_df(design$rest, sample[free, , drop = FALSE],
                  sp$data[sp$rest, , drop = FALSE],
                  match(rows[free], sp$keep[sp$rest]))
      }
    },
    n - 1
  )
  max(out, 0)
}

#' @noRd
ht_variance <- function(design, data, rows, y, pi_i) UseMethod("ht_variance")

#' Sen-Yates-Grundy, for fixed-size designs with closed-form joint
#' probabilities
#' @noRd
ht_variance.default <- function(design, data, rows, y, pi_i) {
  pij <- joint_inclusion(design, data, rows)
  if (length(rows) < 2L) {
    return(list(variance = NA_real_,
                note = "A variance needs at least two sampled rows."))
  }
  syg_variance(y, pi_i, pij)
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
#' @noRd
lonely_note <- function(group, rows, noun, fix) {
  # `noun` is c(singular, plural)
  keep <- !is.na(group)
  N_h <- table(group[keep])
  n_h <- table(group[rows][!is.na(group[rows])])
  lonely <- names(n_h)[n_h == 1L & N_h[names(n_h)] > 1L]
  if (!length(lonely)) return(NULL)
  shown <- group_label(utils::head(lonely, 5L))
  many <- length(lonely) > 1L
  paste0(if (many) paste(length(lonely), noun[2]) else paste("The", noun[1]),
         " ", paste0("`", shown, "`", collapse = ", "),
         if (length(lonely) > 5L) ", ..." else "",
         " ", if (many) "each have" else "has",
         " a single sampled row, so the variation within ",
         if (many) "them" else "it",
         " cannot be measured and leaving it out would understate the standard ",
         "error. ", fix)
}

#' @noRd
ht_variance.drawn_design_stratified <- function(design, data, rows, y, pi_i) {
  group <- as.character(group_key(data, design$strata))
  note <- lonely_note(group, rows, c("stratum", "strata"),
                      "Draw at least two rows per stratum (min_per_stratum = 2).")
  if (!is.null(note)) return(list(variance = NA_real_, note = note))
  NextMethod()
}

#' @noRd
ht_variance.drawn_design_temporal <- function(design, data, rows, y, pi_i) {
  bucket <- temporal_bucket(design, data)
  note <- lonely_note(bucket, rows, c("interval", "intervals"),
                      "Draw at least two rows per interval (per_interval = 2).")
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
ht_variance.drawn_design_cluster <- function(design, data, rows, y, pi_i) {
  cl <- count_clusters(design, data)
  A <- cl$total
  lab <- as.character(cl$labels[rows])
  if (anyNA(lab)) {
    return(list(variance = NA_real_,
                note = "Some sampled rows have no cluster label."))
  }
  u <- vapply(split(y / pi_i, lab), sum, numeric(1))
  m <- length(u)
  if (m < 2L) {
    return(list(variance = NA_real_, note = paste0(
      "A variance needs at least two clusters; this sample has ", m, ".")))
  }
  fpc <- if (is.finite(A) && A > m) 1 - m / A else 0
  list(variance = fpc * m * stats::var(u), note = NULL)
}

#' Two-stage variance, declining where the second stage is unmeasurable
#'
#' With one row per selected cluster, pairs inside a cluster can never be drawn
#' together, so Sen-Yates-Grundy has no way to see the within-cluster variation
#' and understates by around a third. The delete-a-cluster jackknife that
#' `variance = "auto"` falls back to is the ultimate-cluster approximation,
#' which is the standard answer.
#'
#' @noRd
ht_variance.drawn_design_multistage <- function(design, data, rows, y, pi_i) {
  if (!design$replace && design$allocation == "equal" &&
      design$n %/% design$n_clusters == 1L &&
      design$n %% design$n_clusters == 0L) {
    return(list(variance = NA_real_, note = paste0(
      "With one row per selected cluster the variation within clusters ",
      "cannot be measured, and the exact two-stage estimator would understate ",
      "the variance.")))
  }
  NextMethod()
}

#' Variance for the weighted designs
#'
#' Poisson sampling is exact: rows are independent, so no joint matrix is
#' needed. Systematic PPS has no closed-form joint probabilities; Deville's
#' approximation is the standard answer for a high-entropy fixed-size design.
#'
#' @noRd
ht_variance.drawn_design_weighted <- function(design, data, rows, y, pi_i) {
  if (design$method == "poisson") {
    return(list(variance = sum((1 - pi_i) / pi_i^2 * y^2), note = NULL))
  }
  if (design$method == "systematic") {
    return(deville_variance(y, pi_i))
  }
  stop("`design_weighted(method = \"successive\")` has no closed-form ",
       "inclusion probability.", call. = FALSE)
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
  c_k <- 1 - pi_i
  if (sum(c_k) <= 0) {
    return(list(variance = 0, method = "deville", note = paste0(
      "Every sampled row was taken with certainty, so the total is exact ",
      "rather than estimated.")))
  }
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
#' @noRd
ht_variance.drawn_design_systematic <- function(design, data, rows, y, pi_i) {
  if (!is.null(design$start)) {
    stop("A systematic design with a fixed `start` is not a probability ",
         "sample: each row is selected with probability 0 or 1.",
         call. = FALSE)
  }
  n <- length(rows)
  if (n < 2L) {
    return(list(variance = NA_real_,
                note = "A variance needs at least two sampled rows."))
  }
  pos <- systematic_positions(design, data)[rows]
  u <- (y / pi_i)[order(pos)]
  f <- 1 / design$interval
  list(variance = (1 - f) * n / (2 * (n - 1)) * sum(diff(u)^2),
       method = "successive difference",
       note = paste0("Systematic sampling has no design-unbiased variance ",
                     "estimator; this is the successive-difference ",
                     "approximation. It would be understated if the frame ",
                     "cycled with a period matching `interval`, or if a smooth ",
                     "trend along the sort order swamped all other variation."))
}

# ---- jackknife --------------------------------------------------------------

#' Delete-a-group jackknife variance
#'
#' Works from the sample alone, which is what makes it available where the
#' analytic form is not. Primary sampling units are whole clusters where the
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
#' @noRd
jackknife_variance <- function(design, sample, pop, y, pi_i) {
  refusal <- jackknife_refusal(design)
  if (!is.null(refusal)) {
    return(list(variance = NA_real_, method = "jackknife", note = refusal))
  }
  # Certainty rows are in every possible sample. Deleting one and inflating the
  # rest invents variance the design does not have -- four times too much on a
  # sample that is a third certainty -- so hand the whole problem to `rest`.
  if (inherits(design, "drawn_design_certainty")) {
    keep <- pi_i < 1
    if (!any(keep)) {
      return(list(variance = 0, method = "jackknife", note = paste0(
        "Every sampled row was taken with certainty, so the total is exact ",
        "rather than estimated.")))
    }
    sp <- certainty_split(design, pop)
    return(jackknife_variance(design$rest, sample[keep, , drop = FALSE],
                              sp$data[sp$rest, , drop = FALSE],
                              y[keep], pi_i[keep]))
  }

  psu <- jackknife_groups(design, sample)
  st <- jackknife_strata(design, sample, pop, psu)
  base <- y / pi_i
  total <- sum(base)

  v <- 0
  n_psu <- 0L
  lonely <- character(0)
  for (h in unique(st$stratum)) {
    in_h <- st$stratum == h
    g <- psu[in_h]
    ug <- unique(g)
    m <- length(ug)
    f <- st$fpc[[h]]
    n_psu <- n_psu + m
    if (m < 2L) {
      if (f > 0) lonely <- c(lonely, h)
      next
    }
    sums <- vapply(ug, function(u) sum(base[in_h][g == u]), numeric(1))
    reps <- (total - sum(sums)) + (sum(sums) - sums) * m / (m - 1)
    v <- v + f * (m - 1) / m * sum((reps - mean(reps))^2)
  }

  n_strata <- length(unique(st$stratum))
  if (n_strata == 1L && n_psu < 2L) {
    return(list(variance = NA_real_, method = "jackknife",
                note = paste0("The jackknife needs at least two primary ",
                              "sampling units; this sample has ", n_psu, ".")))
  }
  if (length(lonely)) {
    return(list(variance = NA_real_, method = "jackknife", note = paste0(
      "The stratified jackknife needs two units per stratum as well.")))
  }
  list(variance = v, method = "jackknife",
       note = paste0("Jackknife over ", n_psu, " primary sampling unit(s)",
                     if (n_strata > 1L) paste0(" in ", n_strata, " strata")
                     else "",
                     if (any(unlist(st$fpc) < 1)) ", with a finite population correction"
                     else "", "."))
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
#' `jackknife_fpc()`.
#'
#' @noRd
jackknife_strata <- function(design, sample, pop, psu) {
  grp <- switch(design_type(design),
    stratified = list(
      sample = as.character(group_key(sample, design$strata)),
      pop = as.character(group_key(pop, design$strata))
    ),
    temporal = list(
      sample = temporal_bucket(design, sample),
      pop = temporal_bucket(design, pop)
    ),
    NULL
  )
  if (is.null(grp)) {
    m <- length(unique(psu))
    return(list(stratum = rep("all", nrow(sample)),
                fpc = list(all = jackknife_fpc(design, pop, m))))
  }
  N_h <- table(grp$pop[!is.na(grp$pop)])
  n_h <- table(grp$sample)
  fpc <- lapply(stats::setNames(names(n_h), names(n_h)), function(h) {
    max(0, 1 - n_h[[h]] / N_h[[h]])
  })
  list(stratum = grp$sample, fpc = fpc)
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
#' Deliberately not applied to multi-stage designs:
#' the delete-a-cluster jackknife sees only between-cluster variation, so the
#' second stage is already missing, and correcting for the first stage on top of
#' that understates the total. Leaving it out is the ultimate-cluster
#' approximation, which errs conservative -- the usual choice in survey practice.
#'
#' @noRd
jackknife_fpc <- function(design, pop, m) {
  if (inherits(design, "drawn_design_multistage")) return(1)
  total <- tryCatch({
    if (inherits(design, "drawn_design_cluster")) {
      count_clusters(design, pop)$total
    } else if (inherits(design, "drawn_design_reservoir")) {
      reservoir_reach(design, nrow(pop))
    } else if (inherits(design, "drawn_design_spatial")) {
      sum(spatial_inside(design, pop))
    } else {
      nrow(pop)
    }
  }, error = function(e) NA_real_)
  if (!is.finite(total) || total <= m) return(1)
  1 - m / total
}

#' @noRd
jackknife_groups <- function(design, sample) {
  col <- if (inherits(design, c("drawn_design_cluster",
                                "drawn_design_multistage"))) {
    design$clusters
  } else {
    NULL   # rows are the sampling units
  }
  if (!is.null(col) && col %in% names(sample)) {
    as.character(sample[[col]])
  } else {
    as.character(seq_len(nrow(sample)))
  }
}

#' @export
print.drawn_ht <- function(x, ...) {
  print_estimate(x, "total", "Horvitz-Thompson total")
}
