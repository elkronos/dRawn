#' How large a sample do you need?
#'
#' Solves for `n` given the precision you want, rather than asking you to guess
#' it. This is the step before [draw()]: you rarely know the sample size, you
#' know the margin of error you can live with.
#'
#' @section What you need to supply:
#' A margin of error, and some idea of how much the thing you are measuring
#' varies:
#'
#' * For a **mean or total**, an `sd` — from a pilot, last year's data, or a
#'   range divided by four as a rough stand-in.
#' * For a **proportion**, a `p`. Leave it at `0.5`, the most pessimistic value,
#'   unless you have a better guess; that is the honest default because it
#'   maximises the required size. The margin is a proportion too, `0.03` for
#'   3 percentage points.
#'
#' `sd` given for a proportion, or `p` for a mean or a total, is not used,
#' and `plan_size()` warns that it is ignored.
#'
#' @section The corrections:
#' `N` applies a finite population correction: sampling 400 from 500 is very
#' different from 400 from 500,000, and past a certain point a bigger frame
#' stops mattering. `response` inflates for non-response: at `0.6` you draw
#' enough to end up with what you need.
#'
#' `deff` inflates for the design, since a clustered design of 400 may carry
#' the information of 100. Measure it with [deff()] on a comparable past
#' sample, using the estimator you will report: [ht_mean()] for a mean or a
#' proportion, and [ht_total()] for a total. The two design effects differ for
#' the same sample, and can differ many times over. In one cluster sample
#' with clusters of 2 to 40 rows, the total's design effect was 28 and the
#' mean's 0.66, and a total planned with the mean's figure came out a small
#' fraction of the size it needed. A total reported as `N` times a mean takes
#' the mean's design effect, and is planned as a mean with `margin / N`. You
#' can pass the result of `ht_mean()` or `ht_total()` itself as `deff`, and
#' `plan_size()` refuses one of the wrong kind for `target`.
#'
#' A design effect from one pilot is an estimate, and a noisy one, which
#' comes out below the truth more often than above it. Round it up. [deff()]
#' gives `NA` for a pilot without a variance estimate, such as one taken with
#' `variance = "none"`, a census, or one with a single cluster, and for one in
#' which `y` took a single value. Such a pilot cannot supply a design effect,
#' and `plan_size()` refuses an `NA`.
#'
#' @param margin Half-width of the confidence interval you want, the
#'   half-width the plan expects (see "Where the formula comes from"). For a
#'   proportion it is a proportion, and a `margin` of 1 or more is refused.
#' @param sd Standard deviation of the measure, for `target = "mean"` or
#'   `"total"`, as [stats::sd()] computes it.
#' @param p Expected proportion, for `target = "proportion"`.
#' @param N Population size, at least 1. `Inf` for an effectively unbounded
#'   frame. An estimated `N` need not be a whole number; `n` is at most
#'   `floor(N)`.
#' @param level Confidence level.
#' @param deff Design effect to inflate by: a positive number, or the result
#'   of [ht_mean()] for a mean or a proportion, or of [ht_total()] for a
#'   total, on a pilot sample. `1` assumes simple random sampling. A number
#'   should come from [deff()] on the estimator you will report, rounded up
#'   (see "The corrections").
#' @param response Expected response rate, between 0 and 1.
#' @param target `"mean"`, `"proportion"` or `"total"`. For `"total"` the
#'   margin is on the population total and `N` must be finite.
#'
#' @return A list with a `print()` method, holding:
#'   \describe{
#'     \item{`n`}{Draw this many. Never more than `N`, and at least two when
#'       the frame has two rows or more.}
#'     \item{`n_effective`}{What you expect to analyse once `response` has
#'       taken its share.}
#'     \item{`capped`}{`TRUE` when only a census meets the margin, so `n` is
#'       every row of the frame. A census meets any margin when everyone
#'       responds.}
#'     \item{`short`}{`TRUE` when a sample smaller than the frame would meet
#'       the margin if everyone responded, but not at this `response` rate,
#'       so `n` is the whole frame and the margin achieved will be wider than
#'       the one asked for. It is never `TRUE` with `capped`.}
#'     \item{`raised`}{`TRUE` when the margin alone needs one row, and
#'       `n_effective` was raised to two, the fewest that give a standard
#'       error and an interval.}
#'     \item{`margin`, `level`, `N`, `deff`, `response`, `target`}{The inputs,
#'       returned so the assumptions travel with the number. `deff` is the
#'       number used, also when an estimate was passed.}
#'     \item{`p`}{The expected proportion for a proportion, and `NA` for a
#'       mean or a total.}
#'     \item{`n_needed`}{Rows that would reach the margin if everyone
#'       responded, two at least as for `n_effective`. Equal to `n_effective`
#'       unless the whole frame is drawn at a `response` below 1.}
#'     \item{`spread`}{The standard deviation used, `sd`, or
#'       `sqrt(p * (1 - p))` for a proportion.}
#'   }
#'
#' @examples
#' # A proportion, no prior guess, 20,000 in the frame
#' plan_size(margin = 0.03, N = 20000, target = "proportion")
#'
#' # A mean, when a pilot put the spread near 40
#' plan_size(margin = 5, sd = 40, N = 20000)
#'
#' # The same, in a clustered design with 70% response
#' plan_size(margin = 5, sd = 40, N = 20000, deff = 2.5, response = 0.7)
#'
#' # A total, to within 100,000 across a 20,000-row frame
#' plan_size(margin = 1e5, sd = 40, N = 20000, target = "total")
#'
#' # The design effect of a cluster pilot, passed as the estimate itself
#' pop <- data.frame(cl = rep(paste0("c", 1:50), each = 8))
#' pop$y <- rep(seq(10, 59), each = 8) + rep(c(-3, 3), 200)
#' pilot <- draw(pop, design_cluster("cl", n_clusters = 6), seed = 1,
#'               weights = TRUE)
#' plan_size(margin = 3, sd = sd(pilot$y), N = nrow(pop),
#'           deff = ht_mean(pilot, "y"))
#'
#' @section Where the formula comes from:
#' `n0 = (z * spread / margin)^2 * deff` is the size for an unbounded frame,
#' with `z` the normal quantile for `level` and, for a total, `margin / N` as
#' the margin on the mean. For a mean or a total the finite population
#' correction gives `n0 / (1 + n0 / N)`, the smallest `n` with
#' `(1 / n - 1 / N) * deff * sd^2 <= (margin / z)^2`. For a proportion, whose
#' variance `p * (1 - p)` has the divisor `N` rather than `N - 1`, it gives
#' `n0 / (1 + (n0 - 1) / N)`. The result is rounded up, raised to two when it
#' is one, and divided by `response` (Cochran 1977, chapter 4; Valliant,
#' Dever and Kreuter 2018, chapter 3).
#'
#' `margin` is the half-width the plan expects, the usual convention for
#' sizing a sample: the estimate from a sample of size `n` lands within
#' `margin` of the truth with probability about `level`. The interval a
#' sample reports is built from that sample's own estimate of the variance,
#' which comes out above the planned variance about as often as below it, so
#' the reported half-width is wider than `margin` about half the time.
#'
#' That interval also uses t on the design's degrees of freedom rather than
#' `z`, which matters when the primary sampling units are few. With 3, 4 and
#' 6 clusters in one stratum, which give 2, 3 and 5 degrees of freedom, the
#' 95% t quantile is 2.2, 1.6 and 1.3 times `z`. `plan_size()` counts rows
#' and knows nothing of clusters, so when you turn `n` into a few clusters,
#' plan more of them.
#'
#' For a proportion, [ht_mean()] reports an interval on the logit scale, which
#' is asymmetric. When `p` is small its upper arm is longer than its lower one,
#' and exceeds `margin` in most samples even when the estimate is as precise
#' as planned.
#'
#' @references
#' Cochran, W. G. (1977). *Sampling Techniques*, 3rd ed. Wiley.
#'
#' Kish, L. (1965). *Survey Sampling*. Wiley.
#'
#' Valliant, R., Dever, J. A. and Kreuter, F. (2018). *Practical Tools for
#' Designing and Weighting Survey Samples*, 2nd ed. Springer.
#'
#' @seealso [deff()] to measure the design effect of a past sample, [draw()] to
#'   take the sample.
#' @export
plan_size <- function(margin, sd = NULL, p = 0.5, N = Inf, level = 0.95,
                      deff = 1, response = 1,
                      target = c("mean", "proportion", "total")) {
  target <- check_choice(target, c("mean", "proportion", "total"), "target")
  check_pos <- function(x, arg, allow_inf = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 ||
        (!allow_inf && !is.finite(x))) {
      stop("`", arg, "` must be a single positive number",
           if (allow_inf) " (or Inf)" else "", ".", call. = FALSE)
    }
    x
  }
  margin <- check_pos(margin, "margin")
  deff <- check_pos(plan_deff(deff, target), "deff")
  N <- check_pos(N, "N", allow_inf = TRUE)
  if (N < 1) {
    stop("`N` is ", format_bad(N), ", and a frame holds at least one row. ",
         "Give the number of rows in the frame, or Inf.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  if (!is.numeric(response) || length(response) != 1L || is.na(response) ||
      response <= 0 || response > 1) {
    stop("`response` must be a single number in (0, 1].", call. = FALSE)
  }

  if (target == "proportion") {
    if (!is.null(sd)) {
      warn_drawn(paste0(
        "`sd` is not used for target = \"proportion\", whose spread comes ",
        "from `p`. Give the expected proportion as `p`, and leave `sd` out."),
        "drawn_irrelevant_argument")
    }
    if (margin >= 1) {
      stop("`margin` is ", format_bad(margin), ", and a proportion lies ",
           "between 0 and 1, so a margin of 1 or more says nothing. Give it ",
           "as a proportion, such as ", format(margin / 100), " for ",
           format(margin), " percentage point", if (margin != 1) "s", ".",
           call. = FALSE)
    }
    if (!is.numeric(p) || length(p) != 1L || is.na(p) || p < 0 || p > 1) {
      stop("`p` must be a single number in [0, 1].", call. = FALSE)
    }
    spread <- sqrt(p * (1 - p))
    if (spread <= 0) {
      stop("`p` of exactly 0 or 1 implies no variation and no sample is ",
           "needed. Use a value strictly between 0 and 1.", call. = FALSE)
    }
    unit_margin <- margin
    p_used <- p
  } else {
    if (is.null(sd)) {
      stop("`sd` is required for target = \"", target, "\". Use a pilot ",
           "estimate, or a plausible range divided by four.",
           if (!missing(p)) {
             " To plan a proportion from `p`, use target = \"proportion\"."
           }, call. = FALSE)
    }
    if (!missing(p)) {
      warn_drawn(paste0(
        "`p` is not used for target = \"", target, "\", whose spread is ",
        "`sd`. To plan a proportion from `p`, use target = \"proportion\"; ",
        "otherwise leave `p` out."), "drawn_irrelevant_argument")
    }
    spread <- check_pos(sd, "sd")
    if (target == "total") {
      if (!is.finite(N)) {
        stop("A margin on a total needs a finite `N`.", call. = FALSE)
      }
      # A margin of M on the total is a margin of M/N on the mean.
      unit_margin <- margin / N
    } else {
      unit_margin <- margin
    }
    p_used <- NA_real_
  }

  z <- stats::qnorm(1 - (1 - level) / 2)
  n0 <- (z * spread / unit_margin)^2 * deff
  # n0 is the size for an unbounded frame. For a mean or a total, whose sd has
  # the divisor N - 1, the smallest n with (1/n - 1/N) S^2 deff <=
  # (margin/z)^2 is n0 / (1 + n0/N). For a proportion, whose p(1 - p) has the
  # divisor N, it is n0 / (1 + (n0 - 1)/N). Both are written with n0 in a
  # denominator, so an n0 too large for a double gives N, a census.
  n_eff <- if (!is.finite(N)) {
    n0
  } else if (target == "proportion") {
    N / (1 + (N - 1) / n0)
  } else {
    N / (1 + N / n0)
  }
  if (is.na(n_eff) || (!is.finite(N) && !is.finite(n_eff / response))) {
    shown <- format(level)
    if (shown == "1") shown <- format(level, digits = 17)
    stop("`margin` = ", format_bad(margin), " with a spread of ",
         format(signif(spread, 4)), " at `level` = ", shown, " needs more ",
         "rows than a number can hold. Widen `margin` or lower `level`",
         if (!is.finite(N)) ", or give a finite `N`", ".", call. = FALSE)
  }

  # The tolerance keeps a rounding error from adding a row: 21 / 0.7 is
  # 30.000000000000004 in floating point.
  n_margin <- max(1, ceiling(n_eff - 1e-9))
  if (is.finite(N)) n_margin <- min(n_margin, floor(N))
  capped <- is.finite(N) && n_margin >= floor(N)
  # One row meets a margin of about two standard deviations, but gives no
  # variance and no interval.
  raised <- n_margin < 2 && floor(N) >= 2
  n_eff <- if (raised) 2 else n_margin
  n_needed <- n_eff
  n_draw <- ceiling(n_eff / response - 1e-9)

  # Two different situations, and only one of them is a census. `capped`
  # means only a census meets the margin; `short` means a smaller sample
  # would, but not once non-response has taken its share. Either way the
  # frame is drawn whole, and fewer rows come back when response < 1.
  short <- !capped && n_draw > floor(N)
  if (n_draw > floor(N)) {
    n_draw <- floor(N)
    n_eff <- min(n_eff, floor(n_draw * response + 1e-9))
  }

  structure(
    list(n = n_draw, n_effective = n_eff, margin = margin, level = level,
         N = N, deff = deff, response = response, target = target,
         spread = spread, p = p_used, capped = capped, short = short,
         raised = raised, n_needed = n_needed),
    class = "drawn_plan"
  )
}

#' The design effect `plan_size()` inflates by
#'
#' `deff` is a number, or a result of `ht_mean()` or `ht_total()` whose own
#' design effect is used. The estimate must be of the kind `target` plans,
#' since a mean's design effect and a total's differ for the same sample. An
#' `NA`, which `deff()` gives for a pilot without a variance, is refused with
#' that reason. Any other number is left to `plan_size()`'s own check.
#'
#' @noRd
plan_deff <- function(deff, target) {
  na_fix <- paste0("Use a pilot that has a variance, or a design effect from ",
                   "a similar survey.")
  if (inherits(deff, c("drawn_mean", "drawn_ht"))) {
    if (inherits(deff, "drawn_mean") && target == "total") {
      stop("`deff` is a mean from ht_mean(), and `target` is \"total\". ",
           "Pass ht_total() on the same pilot, since the design effects of a ",
           "mean and a total can differ many times over for one sample. A ",
           "total reported as N times a mean is planned as target = \"mean\" ",
           "with `margin` divided by N.", call. = FALSE)
    }
    if (inherits(deff, "drawn_ht") && target != "total") {
      stop("`deff` is a total from ht_total(), and `target` is \"", target,
           "\". Pass ht_mean() on the same pilot, since the design effects ",
           "of a total and a mean can differ many times over for one sample.",
           call. = FALSE)
    }
    value <- deff[["deff"]]
    if (length(value) != 1L || is.na(value)) {
      stop("`deff` has a design effect of NA, which a pilot with no variance ",
           "estimate or with one value of `y` gives. ", na_fix, call. = FALSE)
    }
    return(value)
  }
  if (!is.atomic(deff)) {
    stop("`deff` is a `", class(deff)[1], "`. Give a design effect as a ",
         "single positive number, or the result of ht_mean() or ht_total() ",
         "for one variable without `by`.", call. = FALSE)
  }
  if (length(deff) == 1L && is.na(deff)) {
    stop("`deff` is NA, which deff() gives for a pilot with no variance ",
         "estimate or with one value of `y`. ", na_fix, call. = FALSE)
  }
  deff
}

#' @export
print.drawn_plan <- function(x, ...) {
  cat("Sample size for a ", x$target, "\n", sep = "")
  cat("  draw           ", fmt_n(x$n), "\n", sep = "")
  if (x$n != x$n_effective) {
    cat("  to analyse     ", fmt_n(x$n_effective), "  (after ",
        format(100 * x$response), "% response)\n", sep = "")
  }
  cat("  margin         +/- ",
      format(x$margin, big.mark = ",", scientific = FALSE), " at ",
      format(100 * x$level), "% confidence\n", sep = "")
  spread <- if (identical(x$target, "proportion") && isTRUE(x$p >= 0)) {
    paste0("p ", signif(x$p, 4))
  } else {
    paste0("sd ", signif(x$spread, 4))
  }
  frame <- if (!is.finite(x$N)) {
    "unbounded"
  } else if (x$N == floor(x$N)) {
    fmt_n(x$N)
  } else {
    format(x$N, big.mark = ",")
  }
  cat("  assuming       ", spread,
      if (x$deff != 1) paste0(", deff ", signif(x$deff, 3)) else "",
      ", N ", frame, "\n", sep = "")
  msg <- if (isTRUE(x$capped)) {
    paste0("That is a census. No sample smaller than the frame meets this ",
           "margin, so every row is needed. ",
           if (x$response < 1) paste0(
             "Even then only about ", fmt_n(x$n_effective), " will come back, ",
             "which will fall short of the margin asked for. "),
           "Widen `margin`, or accept the precision a full count gives.")
  } else if (isTRUE(x$short)) {
    paste0("The frame can reach this margin, with ", fmt_n(x$n_needed),
           " rows, but not at ", format(100 * x$response), "% response, ",
           "which would need more rows than the frame holds. Every row is ",
           "drawn, and the margin achieved will be wider than the one asked ",
           "for.")
  } else {
    NULL
  }
  if (!is.null(msg)) {
    cat_note(msg)
  }
  if (isTRUE(x$raised)) {
    cat_note(paste0("One row meets this margin, but a standard error and an ",
                    "interval need two, so the plan is for two."))
  }
  invisible(x)
}

#' Describe a drawn sample
#'
#' What you actually got, against what the design asked for. Worth a look before
#' analysing: it surfaces strata that came up short, weights that vary more than
#' you expected, and rows the design could never have reached.
#'
#' @param sample A data frame returned by [draw()] with `weights = TRUE`.
#'
#' @return A list with a `print()` method, holding:
#'   \describe{
#'     \item{`design`}{The design's type, as a string.}
#'     \item{`n`, `N`}{Rows drawn, and rows in the frame.}
#'     \item{`N_reachable`}{Rows of the frame the design gives a probability
#'       above 0, the population it samples. The printed sampling fraction is
#'       `n / N_reachable`. `NA` if the design has no closed-form inclusion
#'       probability, and the fraction is then `n / N`.}
#'     \item{`expected_n`}{The sample size the design expects, the sum of
#'       the inclusion probabilities over the frame. It is `n` for a design
#'       of fixed size, and printed beside `n` when the two differ, as they
#'       can under Poisson sampling, a [design_systematic()] `interval` that
#'       does not divide the frame, or clusters of unequal size. `NA` where
#'       `N_reachable` is, and for weights conditional on the clusters
#'       selected.}
#'     \item{`weight_range`}{The smallest and largest design weight.}
#'     \item{`weight_cv`}{Their coefficient of variation. Large values mean a
#'       few rows carry most of the estimate.}
#'     \item{`n_kish`}{Kish's effective sample size, `sum(w)^2 / sum(w^2)`
#'       over the weights `w`, which is `n / (1 + cv^2)` with the
#'       coefficient of variation computed with divisor `n`: the size of a
#'       simple random sample that would estimate a mean as precisely if
#'       unequal weights were the only cost. Clustering and stratification
#'       are not in it; [deff()] measures those.}
#'     \item{`n_certain`}{Sampled rows taken with certainty, at a probability
#'       of 1.}
#'     \item{`certain_by`}{Those rows by why they are certain, a named count:
#'       `threshold`, above the threshold of [design_certainty()]; `size`,
#'       put at probability 1 by a size measure, in [design_weighted()] or
#'       in [design_spread()] or [design_multistage()] with `size`; and
#'       `whole`, in a stratum, interval or cluster taken whole, or a
#'       census.}
#'     \item{`conditional`}{The number of clusters selected, when the
#'       weights are conditional on them, as for [design_cluster()] with
#'       `balanced = TRUE` (see "Design weights" in [draw()]); otherwise
#'       `NA`.}
#'     \item{`unreachable`}{Frame rows with inclusion probability 0 — the
#'       design could never have selected them. `NA` if the design has no
#'       closed-form inclusion probability.}
#'     \item{`by_group`}{A data frame of `group`, `drawn`, `in_frame` and
#'       `rate` per stratum, interval or cluster, or `NULL` for a design with
#'       no grouping. The strata are those of [design_stratified()], and of
#'       [design_spread()] and [design_weighted()] when they have `strata`,
#'       and every stratum of the frame is listed, one allocated no rows
#'       included. The intervals of [design_temporal()] are listed in time
#'       order and named by their start, such as `"2024-01-03 00:00"`. For
#'       [design_cluster()] and [design_multistage()], only the clusters
#'       drawn are listed. For [design_certainty()], the groups are those of
#'       its `rest` over the rows below the threshold, after a `"(certain)"`
#'       row for the rows above it.}
#'     \item{`group_col`}{The column(s) `by_group` is keyed on.}
#'   }
#'
#' @examples
#' set.seed(1)
#' pop <- data.frame(
#'   id = 1:400,
#'   site = rep(c("a", "b", "c", "d"), times = c(200, 100, 60, 40))
#' )
#' s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
#' sample_summary(s)
#'
#' @seealso [draw()], [deff()]
#' @export
sample_summary <- function(sample) {
  state <- drawn_state(sample)
  design <- state$design
  pop <- state$pop
  w <- sample$.weight
  groups <- summary_groups(design, pop, state$rows)
  certain <- is_certain(sample$.prob)

  pi_all <- state$pi_all %||%
    tryCatch(exact_inclusion(design, pop), error = function(e) NULL)

  structure(list(
    design = design_type(design),
    n = nrow(sample),
    N = nrow(pop),
    N_reachable = if (is.null(pi_all)) {
      NA_integer_
    } else {
      sum(reachable(pi_all))
    },
    # Weights conditional on the clusters selected leave the probabilities
    # of the rest of the frame unknown (NA), and their sum with them.
    expected_n = if (is.null(pi_all)) NA_real_ else sum(pi_all),
    weight_range = if (length(w)) range(w, na.rm = TRUE) else c(NA, NA),
    weight_cv = if (length(w) > 1L && mean(w, na.rm = TRUE) > 0) {
      stats::sd(w, na.rm = TRUE) / mean(w, na.rm = TRUE)
    } else NA_real_,
    n_kish = if (length(w) && sum(w^2, na.rm = TRUE) > 0) {
      sum(w, na.rm = TRUE)^2 / sum(w^2, na.rm = TRUE)
    } else NA_real_,
    n_certain = sum(certain),
    certain_by = certain_counts(design, pop, state$rows, certain),
    conditional = if (is.null(state$selection)) {
      NA_integer_
    } else {
      length(state$selection)
    },
    unreachable = if (is.null(pi_all)) {
      NA_integer_
    } else {
      sum(!reachable(pi_all))
    },
    by_group = groups$table,
    group_col = groups$col
  ), class = "drawn_summary")
}

#' Rows drawn against rows in the frame, per stratum, interval or cluster
#'
#' Keyed on the frame and counted at the sampled frame rows `rows`, so a
#' stratum or cluster column renamed or recoded in the sample after the draw
#' still counts the groups that were drawn. Strata and intervals are listed
#' whether drawn or not, so one allocated no rows shows as 0; clusters only
#' when drawn, since most clusters of a frame are not. A certainty design's
#' groups are its `rest`'s over the rows below the threshold, after a row for
#' those above it. Returns `NULL` for a design with no grouping, otherwise
#' `table` and `col`, the column(s) keyed on.
#'
#' @noRd
summary_groups <- function(design, pop, rows) {
  if (inherits(design, "drawn_design_certainty")) {
    sp <- certainty_parts(design, pop)
    top <- sp$keep[sp$take]
    pos <- match(rows, sp$keep[sp$rest])
    inner <- summary_groups(design$rest, sp$rest_data, pos[!is.na(pos)])
    if (is.null(inner)) return(NULL)
    got <- sum(rows %in% top)
    inner$table <- rbind(
      data.frame(group = "(certain)", drawn = got, in_frame = length(top),
                 rate = if (length(top)) round(got / length(top), 4) else NA,
                 row.names = NULL),
      inner$table)
    return(inner)
  }
  type <- design_type(design)
  col <- switch(type,
    stratified = ,
    spread = ,
    weighted = design$strata,
    cluster = ,
    multistage = design$clusters,
    temporal = design$time,
    NULL
  )
  if (is.null(col)) return(NULL)
  key <- if (type == "temporal") {
    temporal_bucket(design, pop)
  } else {
    group_key(pop, col)
  }
  if (walks_strata(design)) {
    # Rows the walk leaves out for a missing `order_by` are in no stratum.
    key[is.na(stratum_allocation(design, pop)$code)] <- NA
  }
  code <- as.integer(key)
  got <- tabulate(code[rows], nlevels(key))
  have <- tabulate(code, nlevels(key))
  shown <- if (type %in% c("cluster", "multistage")) {
    which(got > 0L)
  } else {
    which(have > 0L)
  }
  list(table = data.frame(
    group = group_label(levels(key)[shown]),
    drawn = got[shown],
    in_frame = have[shown],
    rate = round(got[shown] / have[shown], 4),
    row.names = NULL
  ), col = col)
}

#' Sampled rows taken with certainty, by why they are certain
#'
#' `certain` marks the sampled rows at probability 1 and `rows` gives their
#' frame rows. Rows above a certainty design's threshold count as
#' `threshold`; the others are put there by the design below it, a size
#' measure (`size`) or a group taken whole (`whole`).
#'
#' @noRd
certain_counts <- function(design, pop, rows, certain) {
  out <- c(threshold = 0L, size = 0L, whole = 0L)
  while (inherits(design, "drawn_design_certainty")) {
    sp <- certainty_parts(design, pop)
    out[["threshold"]] <- out[["threshold"]] +
      sum(certain & rows %in% sp$keep[sp$take])
    pos <- match(rows, sp$keep[sp$rest])
    certain <- certain[!is.na(pos)]
    rows <- pos[!is.na(pos)]
    pop <- sp$rest_data
    design <- design$rest
  }
  by_size <- inherits(design, "drawn_design_weighted") ||
    (inherits(design, c("drawn_design_spread", "drawn_design_multistage")) &&
       !is.null(design$size))
  src <- if (by_size) "size" else "whole"
  out[[src]] <- out[[src]] + sum(certain)
  out
}

#' @export
print.drawn_summary <- function(x, ...) {
  cat("Sample of ", fmt_n(x$n), " from ", fmt_n(x$N), "  (", x$design,
      " design)\n", sep = "")
  # The fraction of the rows the design could draw: rows at probability 0
  # are listed below as unreachable, and counting them in would understate
  # it.
  reach <- x$N_reachable
  partial <- isTRUE(reach < x$N)
  cat("  sampling fraction  ",
      signif(x$n / if (partial) reach else x$N, 3),
      if (partial) paste0("  (of the ", fmt_n(reach), " reachable rows)"),
      "\n", sep = "")
  if (isTRUE(abs(x$expected_n - x$n) > 1e-6 * max(1, x$n))) {
    cat("  expected n         ", signif(x$expected_n, 4),
        "  (the sample size is random under this design)\n", sep = "")
  }
  if (all(is.finite(x$weight_range))) {
    cat("  design weights     ", signif(x$weight_range[1], 4), " to ",
        signif(x$weight_range[2], 4), sep = "")
    if (is.finite(x$weight_cv)) {
      cat("   (cv ", signif(x$weight_cv, 3), ")", sep = "")
    }
    cat("\n")
  }
  if (isTRUE(is.finite(x$n_kish))) {
    cat("  effective n        ", signif(x$n_kish, 4),
        "  (Kish, from the weights alone)\n", sep = "")
  }
  if (isTRUE(x$n_certain > 0)) {
    why <- c(threshold = "above the threshold", size = "by size",
             whole = "taken whole")
    k <- x$certain_by[x$certain_by > 0]
    cat("  certain rows       ", fmt_n(x$n_certain), sep = "")
    if (length(k)) {
      cat("  (", paste(fmt_n(k), why[names(k)], collapse = ", "), ")",
          sep = "")
    }
    cat("\n")
  }
  if (isTRUE(x$conditional > 0)) {
    cat("                     conditional on the ", x$conditional,
        " clusters selected\n", sep = "")
  }
  if (isTRUE(x$unreachable > 0)) {
    cat("  unreachable rows   ", fmt_n(x$unreachable),
        "  (probability 0 under this design)\n", sep = "")
  }
  if (!is.null(x$by_group)) {
    cat("\n  by ", paste(x$group_col, collapse = " x "), ":\n", sep = "")
    b <- x$by_group
    w1 <- max(nchar(c(b$group, "group")))
    cat("    ", formatC("group", width = -w1), "  drawn  in frame   rate\n",
        sep = "")
    for (i in seq_len(nrow(b))) {
      cat("    ", formatC(b$group[i], width = -w1),
          formatC(b$drawn[i], width = 7),
          formatC(b$in_frame[i], width = 9),
          formatC(formatC(b$rate[i], format = "f", digits = 3), width = 7),
          "\n", sep = "")
    }
  }
  invisible(x)
}
