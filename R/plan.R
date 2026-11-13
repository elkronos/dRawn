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
#' For a stratified sample, give the size and the spread of each stratum in
#' `strata` instead (see "Planning a stratified sample").
#'
#' @section The corrections:
#' `N` applies a finite population correction: sampling 400 from 500 is very
#' different from 400 from 500,000, and past a certain point a bigger frame
#' stops mattering. `response` inflates for non-response: at `0.6` you draw
#' enough to end up with what you need.
#'
#' `deff` inflates for the design, since a clustered design of 400 may carry
#' the information of 100. Before any pilot, [deff_cluster()] gives the
#' design effect of taking `m` rows from each cluster from their intraclass
#' correlation. Measure it with [deff()] on a comparable past
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
#'   Several margins give several plans (see "Several plans at once").
#' @param sd Standard deviation of the measure, for `target = "mean"` or
#'   `"total"`, as [stats::sd()] computes it. Not used with `strata`.
#' @param p Expected proportion, for `target = "proportion"`. Not used with
#'   `strata`.
#' @param N Population size, at least 1. `Inf` for an effectively unbounded
#'   frame. An estimated `N` need not be a whole number; `n` is at most
#'   `floor(N)`. Not used with `strata`, whose sizes add up to it.
#' @param level Confidence level.
#' @param deff Design effect to inflate by: a positive number, or the result
#'   of [ht_mean()] for a mean or a proportion, or of [ht_total()] for a
#'   total, on a pilot sample. `1` assumes simple random sampling. A number
#'   should come from [deff()] on the estimator you will report, rounded up,
#'   or from [deff_cluster()] (see "The corrections"). With `strata` it
#'   inflates the variance of the stratified sample, as for clusters sampled
#'   within strata.
#' @param response Expected response rate, between 0 and 1.
#' @param target `"mean"`, `"proportion"` or `"total"`. For `"total"` the
#'   margin is on the population total and `N` must be finite.
#' @param strata `NULL`, or a data frame with a row per stratum, to plan a
#'   stratified sample (see "Planning a stratified sample"). Its columns
#'   are those that name the strata; `N`, the rows of the frame in each
#'   stratum; `sd`, the standard deviation of the measure within each
#'   stratum, or for a proportion `p`, the expected proportion within it;
#'   and, with `allocation = "neyman"` only, an optional `cost`, the cost
#'   of one row in each stratum.
#' @param allocation How `n` is split over `strata`: `"proportional"`,
#'   `"equal"` or `"neyman"`, as in [design_stratified()]. Neyman shares are
#'   proportional to `N_h * S_h`, or `N_h * S_h / sqrt(cost_h)` with a
#'   `cost` column.
#' @param min_per_stratum The fewest rows a stratum is allocated, as in
#'   [design_stratified()]. A stratum smaller than this is taken whole. The
#'   default of 2 is the fewest that give each stratum a variance, which
#'   [ht_total()] and [ht_mean()] need from every stratum not taken whole.
#'
#' @return A plan, a list of class `drawn_plan` with a `print()` method,
#'   or for several plans a data frame (see "Several plans at once"). A plan
#'   holds:
#'   \describe{
#'     \item{`n`}{Draw this many. Never more than `N`, and without `strata`
#'       at least two when the frame has two rows or more.}
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
#'       `sqrt(p * (1 - p))` for a proportion. `NA` with `strata`.}
#'   }
#'
#'   With `strata` there are four more, and `raised` is always `FALSE`,
#'   since `min_per_stratum` sets the floor:
#'   \describe{
#'     \item{`allocation`}{The rows to draw from each stratum, a data frame
#'       of the strata columns and `n`, in the order of `strata`. Pass it to
#'       [design_stratified()] as `allocation` to draw that split.}
#'     \item{`strata`}{The strata as given, checked, with `n` added.}
#'     \item{`rule`, `min_per_stratum`}{The `allocation` rule and the floor
#'       used.}
#'   }
#'
#'   A plan can be passed as `n` to any design that takes `n`, such as
#'   `design_simple(n = plan)`, which draws `plan$n` rows. `as.data.frame()`
#'   turns it into a data frame of one row, with a column for each field
#'   above, except that a stratified plan has the number of `strata` in
#'   place of the two data frames. Several plans at once give a data frame
#'   of those rows, one per plan.
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
#' # Clusters of ten rows, whose rows correlate at 0.05, before any pilot
#' plan_size(margin = 5, sd = 40, N = 20000,
#'           deff = deff_cluster(m = 10, icc = 0.05))
#'
#' # Three margins at once, one plan per row
#' plan_size(margin = c(3, 5, 10), sd = 40, N = 20000)
#'
#' # A plan is a valid `n`
#' design_simple(n = plan_size(margin = 5, sd = 40, N = 20000))
#'
#' # Three strata, with the spread of each from an earlier survey
#' bands <- data.frame(band = c("small", "medium", "large"),
#'                     N = c(1200, 600, 200), sd = c(10, 25, 80))
#' plan <- plan_size(margin = 2, strata = bands, allocation = "neyman")
#' plan
#' frame <- data.frame(id = 1:2000,
#'                     band = rep(bands$band, times = bands$N))
#' s <- draw(frame, design_stratified("band", allocation = plan$allocation),
#'           seed = 1)
#' table(s$band)
#'
#' # The design effect of a cluster pilot, passed as the estimate itself
#' pop <- data.frame(cl = rep(paste0("c", 1:50), each = 8))
#' pop$y <- rep(seq(10, 59), each = 8) + rep(c(-3, 3), 200)
#' pilot <- draw(pop, design_cluster("cl", n_clusters = 6), seed = 1,
#'               weights = TRUE)
#' plan_size(margin = 3, sd = sd(pilot$y), N = nrow(pop),
#'           deff = ht_mean(pilot, "y"))
#'
#' @section Planning a stratified sample:
#' With `strata`, `plan_size()` plans a stratified simple random sample, as
#' [design_stratified()] draws one. For each total `n` it splits `n` over
#' the strata as that design would, by the same `allocation` rule, the
#' floor `min_per_stratum` and largest-remainder rounding, with no stratum
#' given more rows than it has, and works out the variance of the
#' stratified mean,
#' `deff * sum(W_h^2 * (1 / n_h - 1 / N_h) * S_h^2)` with `W_h = N_h / N`
#' (Cochran 1977, chapter 5). For a proportion `S_h^2` is
#' `N_h * p_h * (1 - p_h) / (N_h - 1)`. The plan is the smallest `n` whose
#' variance is at most `(margin / z)^2`, with `margin / N` for a total. A
#' split that leaves a stratum with no rows does not count, since the
#' estimate then leaves that stratum out. The rounded split can take a row
#' from one stratum as `n` grows, so the variance does not always fall with
#' `n`, and the search checks each `n` from the first that could meet the
#' margin.
#'
#' `response` divides `n` as without strata, and the rows to draw are split
#' by the same rule. The split is returned as `allocation`, which
#' [design_stratified()] takes as the count for each stratum. Passing the
#' plan as `n` instead lets the design split `n` again by its own rule and
#' `min_per_stratum`, whose default there is 0.
#'
#' @section Several plans at once:
#' `margin`, `sd`, `p`, `deff` and `response` each take one value, or
#' several for several plans. Those with one value are used for every plan,
#' and those with several must have the same number, one per plan, so
#' `plan_size(margin = c(3, 5, 10), sd = 40)` plans three margins at the
#' same spread. The result is then a data frame with a row per plan, as
#' `as.data.frame()` gives for each, without the split over strata or the
#' printed notes.
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
#' @seealso [deff()] to measure the design effect of a past sample,
#'   [deff_cluster()] for one from an intraclass correlation, [draw()] to
#'   take the sample.
#' @export
plan_size <- function(margin, sd = NULL, p = 0.5, N = Inf, level = 0.95,
                      deff = 1, response = 1,
                      target = c("mean", "proportion", "total"),
                      strata = NULL,
                      allocation = c("proportional", "equal", "neyman"),
                      min_per_stratum = 2L) {
  target <- check_choice(target, c("mean", "proportion", "total"), "target")
  check_pos <- function(x, arg, allow_inf = FALSE, several = TRUE) {
    if (!is.numeric(x) || !length(x) || (!several && length(x) != 1L) ||
        anyNA(x) || any(x <= 0) || (!allow_inf && !all(is.finite(x)))) {
      stop("`", arg, "` must be ",
           if (several) "a positive number" else "a single positive number",
           if (allow_inf) " (or Inf)" else "",
           if (several) ", or a vector of them for several plans" else "",
           ".", call. = FALSE)
    }
    x
  }
  margin <- check_pos(margin, "margin")
  deff <- check_pos(plan_deff(deff, target), "deff")
  if (is.null(strata)) {
    given <- c("allocation", "min_per_stratum")[
      c(!missing(allocation), !missing(min_per_stratum))]
    if (length(given)) {
      warn_drawn(paste0(
        paste0("`", given, "`", collapse = " and "), " ",
        if (length(given) == 1L) "is" else "are", " not used without ",
        "`strata`, which describes the strata to split `n` over. Give a ",
        "row per stratum in `strata`, or leave ",
        if (length(given) == 1L) "it" else "them", " out."),
        "drawn_irrelevant_argument")
    }
    if (is.numeric(N) && length(N) > 1L) {
      stop("`N` must be a single positive number (or Inf), not ", length(N),
           " numbers. For a frame in strata, give each stratum's `N` and ",
           "`sd` in `strata`, a data frame with a row per stratum.",
           call. = FALSE)
    }
    N <- check_pos(N, "N", allow_inf = TRUE, several = FALSE)
    if (N < 1) {
      stop("`N` is ", format_bad(N), ", and a frame holds at least one row. ",
           "Give the number of rows in the frame, or Inf.", call. = FALSE)
    }
    st <- NULL
  } else {
    given <- c("sd", "p", "N")[c(!is.null(sd), !missing(p), !missing(N))]
    if (length(given)) {
      warn_drawn(paste0(
        paste0("`", given, "`", collapse = " and "), " ",
        if (length(given) == 1L) "is" else "are", " not used with ",
        "`strata`, whose columns give each stratum's `N` and its `",
        if (target == "proportion") "p" else "sd", "`. Leave ",
        if (length(given) == 1L) "it" else "them", " out."),
        "drawn_irrelevant_argument")
    }
    allocation <- check_choice(allocation,
                               c("proportional", "equal", "neyman"),
                               "allocation")
    min_per_stratum <- check_count(min_per_stratum, "min_per_stratum")
    st <- plan_strata(strata, target, allocation, min_per_stratum)
    N <- st$N
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  if (!is.numeric(response) || !length(response) || anyNA(response) ||
      any(response <= 0 | response > 1)) {
    stop("`response` must be a number in (0, 1], or a vector of them for ",
         "several plans.", call. = FALSE)
  }

  if (target == "proportion") {
    if (is.null(st) && !is.null(sd)) {
      warn_drawn(paste0(
        "`sd` is not used for target = \"proportion\", whose spread comes ",
        "from `p`. Give the expected proportion as `p`, and leave `sd` out."),
        "drawn_irrelevant_argument")
    }
    if (any(margin >= 1)) {
      m <- margin[margin >= 1][1]
      stop("`margin` is ", format_bad(m), ", and a proportion lies ",
           "between 0 and 1, so a margin of 1 or more says nothing. Give it ",
           "as a proportion, such as ", format(m / 100), " for ",
           format(m), " percentage point", if (m != 1) "s", ".",
           call. = FALSE)
    }
    if (is.null(st)) {
      if (!is.numeric(p) || !length(p) || anyNA(p) || any(p < 0 | p > 1)) {
        stop("`p` must be a number in [0, 1], or a vector of them for ",
             "several plans.", call. = FALSE)
      }
      spread <- sqrt(p * (1 - p))
      if (any(spread <= 0)) {
        stop("`p` of exactly 0 or 1 implies no variation and no sample is ",
             "needed. Use a value strictly between 0 and 1.", call. = FALSE)
      }
      p_used <- p
    } else {
      spread <- NA_real_
      p_used <- NA_real_
    }
  } else {
    if (is.null(st)) {
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
    } else {
      spread <- NA_real_
    }
    if (target == "total" && !is.finite(N)) {
      stop("A margin on a total needs a finite `N`.", call. = FALSE)
    }
    p_used <- NA_real_
  }

  several <- list(margin = margin, deff = deff, response = response)
  if (is.null(st)) {
    several[[if (target == "proportion") "p" else "sd"]] <- spread
  }
  k <- plan_count(several)
  one <- function(i) {
    at <- function(x) x[[if (length(x) == 1L) 1L else i]]
    plan_one(at(margin), at(spread), at(p_used), N, level, at(deff),
             at(response), target, st)
  }
  if (k == 1L) return(one(1L))
  do.call(rbind, lapply(seq_len(k), function(i) as.data.frame(one(i))))
}

#' How many plans the vector arguments of `plan_size()` ask for
#'
#' Each of `args` holds one value, or one per plan; one value is used for
#' every plan, as [base::Map()] recycles it. Lengths that differ otherwise
#' are refused, since recycling them would pair values by accident.
#'
#' @noRd
plan_count <- function(args) {
  lens <- lengths(args)
  k <- max(lens)
  bad <- lens != 1L & lens != k
  if (any(bad)) {
    long <- names(lens)[lens == k][1]
    other <- names(lens)[bad][1]
    stop("`", long, "` has ", k, " values and `", other, "` has ",
         lens[[other]], ". Give each of ",
         paste0("`", utils::head(names(args), -1L), "`", collapse = ", "),
         " and `", utils::tail(names(args), 1L), "` one value, or one for ",
         "each plan.", call. = FALSE)
  }
  k
}

#' One plan from single values of each argument
#'
#' `plan_size()` has checked the arguments. `st` is `NULL`, or the strata
#' from `plan_strata()`.
#'
#' @noRd
plan_one <- function(margin, spread, p, N, level, deff, response, target,
                     st) {
  # A margin of M on the total is a margin of M/N on the mean.
  unit_margin <- if (target == "total") margin / N else margin
  z <- stats::qnorm(1 - (1 - level) / 2)
  if (is.null(st)) {
    n0 <- (z * spread / unit_margin)^2 * deff
    # n0 is the size for an unbounded frame. For a mean or a total, whose sd
    # has the divisor N - 1, the smallest n with (1/n - 1/N) S^2 deff <=
    # (margin/z)^2 is n0 / (1 + n0/N). For a proportion, whose p(1 - p) has
    # the divisor N, it is n0 / (1 + (n0 - 1)/N). Both are written with n0
    # in a denominator, so an n0 too large for a double gives N, a census.
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
  } else {
    n_margin <- strata_size(st, (unit_margin / z)^2 / deff)
    capped <- n_margin >= N
    # The floor `min_per_stratum` takes the place of the two rows.
    raised <- FALSE
    n_eff <- n_margin
  }
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

  out <- list(n = n_draw, n_effective = n_eff, margin = margin,
              level = level, N = N, deff = deff, response = response,
              target = target, spread = spread, p = p, capped = capped,
              short = short, raised = raised, n_needed = n_needed)
  if (!is.null(st)) {
    tab <- st$table
    tab$n <- integer(nrow(tab))
    tab$n[st$order] <- strata_split(st, n_draw)
    out$rule <- st$rule
    out$min_per_stratum <- st$floor
    out$strata <- tab
    out$allocation <- tab[c(st$keys, "n")]
  }
  structure(out, class = "drawn_plan")
}

#' Check the strata given to `plan_size()`, and put them in key order
#'
#' `strata` is a data frame of the columns that name the strata, `N`, `sd`
#' (or `p` for a proportion) and, for Neyman allocation, an optional
#' `cost`. Returns the checked table in the order given, and in key order,
#' the order `design_stratified()` allocates in, the sizes, the variances
#' `S_h^2`, the spreads and costs for Neyman allocation, and `w2`, each
#' stratum's squared share of the frame.
#'
#' @noRd
plan_strata <- function(strata, target, allocation, min_per_stratum) {
  spread_col <- if (target == "proportion") "p" else "sd"
  if (!is.data.frame(strata)) {
    stop("`strata` must be a data frame with a row per stratum, holding the ",
         "columns that name the strata, `N` and `", spread_col, "`, not ",
         format_bad(strata), ".", call. = FALSE)
  }
  given <- names(strata)
  if ("n" %in% given) {
    stop("`strata` has a column `n`, which is what plan_size() works out. ",
         "Give the number of rows of the frame in each stratum as `N`.",
         call. = FALSE)
  }
  keys <- setdiff(given, c("N", "sd", "p", "cost"))
  if (!length(keys)) {
    stop("`strata` has no column naming the strata. Add one, such as the ",
         "column `design_stratified()` will take as `strata`, beside `N` ",
         "and `", spread_col, "`.", call. = FALSE)
  }
  need <- setdiff(c("N", spread_col), given)
  if (length(need)) {
    what <- c(N = "the number of rows of the frame in each stratum",
              sd = paste0("the standard deviation of the measure in each ",
                          "stratum, from a pilot or an earlier survey"),
              p = "the expected proportion in each stratum")
    stop("`strata` needs a column `", need[1], "`, ", what[[need[1]]], ".",
         call. = FALSE)
  }
  unused <- setdiff(intersect(c("sd", "p"), given), spread_col)
  if (length(unused)) {
    warn_drawn(paste0(
      "`strata` has a `", unused, "` column, which target = \"", target,
      "\" does not use, since its spread is the `", spread_col, "` column. ",
      "Leave it out", if (unused == "p") {
        ", or use target = \"proportion\""
      }, "."), "drawn_irrelevant_argument")
  }
  has_cost <- "cost" %in% given
  if (has_cost && allocation != "neyman") {
    stop("`strata` has a `cost` column, which applies only to allocation = ",
         "\"neyman\", and `allocation` is \"", allocation, "\". Use ",
         "allocation = \"neyman\", or drop the column.", call. = FALSE)
  }
  value <- c("N", spread_col, if (has_cost) "cost")
  spec <- strata_frame(strata, keys, "strata", value)

  sizes <- spec$N
  if (!is.numeric(sizes)) {
    stop("The stratum sizes in `strata$N` must be numbers, not ",
         class(sizes)[1], ".", call. = FALSE)
  }
  bad <- !is.finite(sizes) | sizes < 1 | sizes != trunc(sizes) |
    sizes > .Machine$integer.max
  if (any(bad)) {
    h <- which(bad)[1]
    stop("The stratum sizes in `strata$N` must be whole numbers of 1 or ",
         "more, the rows of the frame in each stratum, and stratum ",
         spec_label(spec, keys, h), " has ", format(sizes[h]), ".",
         call. = FALSE)
  }
  spec$N <- as.integer(sizes)
  sizes <- as.numeric(sizes)
  if (spread_col == "sd") {
    spread <- check_stratum_values(spec[c(keys, "sd")], keys, "strata",
                                   "sd")$sd
    spec$sd <- spread
    s2 <- spread^2
  } else {
    p <- spec$p
    if (!is.numeric(p)) {
      stop("The proportions in `strata$p` must be numbers, not ",
           class(p)[1], ".", call. = FALSE)
    }
    bad <- is.na(p) | p < 0 | p > 1
    if (any(bad)) {
      h <- which(bad)[1]
      stop("The proportions in `strata$p` must lie between 0 and 1, and ",
           "stratum ", spec_label(spec, keys, h), " has ", format(p[h]), ".",
           call. = FALSE)
    }
    # p(1 - p) has the divisor N_h; S_h^2 has N_h - 1. A stratum of one row
    # is taken whole and adds nothing.
    s2 <- ifelse(sizes > 1, sizes * p * (1 - p) / pmax(sizes - 1, 1), 0)
    spread <- sqrt(s2)
  }
  cost <- if (has_cost) {
    spec$cost <- check_stratum_values(spec[c(keys, "cost")], keys, "strata",
                                      "cost", positive = TRUE)$cost
    spec$cost
  }

  o <- order(as.integer(group_key(spec, keys)))
  # Neyman shares do not change when every spread is scaled alike, and a
  # spread near the largest double would overflow `size * spread`.
  if (any(spread > 0)) spread <- spread / max(spread)
  N <- sum(sizes)
  list(table = spec, keys = keys, order = o, sizes = sizes[o],
       s2 = s2[o], w2 = (sizes[o] / N)^2,
       spread = if (allocation == "neyman") spread[o],
       cost = if (!is.null(cost)) cost[o],
       rule = allocation, floor = min_per_stratum, N = N)
}

#' The split of `n` over the strata of a plan, in key order
#'
#' As `design_stratified()` splits it: the same rule, floor and
#' largest-remainder rounding, with no stratum given more rows than it has.
#'
#' @noRd
strata_split <- function(st, n) {
  allocate(n, st$sizes, st$rule, st$floor, cap = TRUE, spread = st$spread,
           cost = st$cost)
}

#' Variance of a stratified mean for a split `n_h` of a plan's strata
#'
#' `sum(W_h^2 * (1 / n_h - 1 / N_h) * S_h^2)`. A stratum taken whole adds
#' nothing, even when its `S_h^2` is too large for a double. A split that
#' leaves a stratum with no rows gives `Inf`, since the estimate then leaves
#' the stratum out and no margin holds for it.
#'
#' @noRd
strata_variance <- function(st, n_h) {
  if (any(n_h < 1)) return(Inf)
  sum(st$w2 * ifelse(n_h >= st$sizes, 0, (1 / n_h - 1 / st$sizes) * st$s2))
}

#' The smallest total whose split meets a variance bound
#'
#' The variance of the rounded split need not fall as `n` grows, since
#' largest-remainder rounding can give a stratum one row fewer at `n + 1`
#' than at `n`. So the search scans upward from the first `n` that a lower
#' bound on the variance allows. Rounding gives each stratum at most one row
#' more than its share before rounding, and those shares never fall as `n`
#' grows, so the variance at the shares plus one row is a lower bound that
#' falls with `n`, and bisection finds the first `n` it allows. Every
#' smaller `n` has a variance above `bound`, and the scan ends at the census
#' at the latest, whose variance is 0.
#'
#' @noRd
strata_size <- function(st, bound) {
  limit <- bound * (1 + 1e-9)
  share <- allocation_share(st$rule, st$sizes, st$spread, st$cost)
  lower <- group_floor(st$floor, st$sizes, cap = TRUE)
  floor_at <- function(n) {
    raw <- capped_shares(n, share, st$sizes, cap = TRUE, lower = lower)
    strata_variance(st, pmin(raw + 1, st$sizes))
  }
  lo <- sum(lower)
  hi <- st$N
  if (floor_at(lo) > limit) {
    while (hi - lo > 1) {
      mid <- floor((lo + hi) / 2)
      if (floor_at(mid) <= limit) hi <- mid else lo <- mid
    }
    lo <- hi
  }
  n <- lo
  while (n < st$N && strata_variance(st, strata_split(st, n)) > limit) {
    n <- n + 1
  }
  n
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
         "positive number, a vector of them for several plans, or the result ",
         "of ht_mean() or ht_total() for one variable without `by`.",
         call. = FALSE)
  }
  if (length(deff) && anyNA(deff)) {
    stop("`deff` ", if (length(deff) == 1L) "is NA" else "holds an NA",
         ", which deff() gives for a pilot with no variance estimate or with ",
         "one value of `y`. ", na_fix, call. = FALSE)
  }
  deff
}

#' @export
print.drawn_plan <- function(x, ...) {
  by_strata <- !is.null(x$allocation)
  cat("Sample size for a ", x$target,
      if (by_strata) paste0(", over ", nrow(x$allocation), " strata"), "\n",
      sep = "")
  cat("  draw           ", fmt_n(x$n), "\n", sep = "")
  if (x$n != x$n_effective) {
    cat("  to analyse     ", fmt_n(x$n_effective), "  (after ",
        format(100 * x$response), "% response)\n", sep = "")
  }
  cat("  margin         +/- ",
      format(x$margin, big.mark = ",", scientific = FALSE), " at ",
      format(100 * x$level), "% confidence\n", sep = "")
  spread <- if (by_strata) {
    paste0(if (x$rule == "neyman") "Neyman" else x$rule, " allocation",
           if (x$min_per_stratum > 0L) {
             paste0(", at least ", x$min_per_stratum, " row",
                    if (x$min_per_stratum != 1L) "s", " per stratum")
           })
  } else if (identical(x$target, "proportion") && isTRUE(x$p >= 0)) {
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
  if (by_strata) print_strata(x$strata)
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

#' The strata of a plan, ten at most, with the rows to draw from each
#' @noRd
print_strata <- function(tab, max = 10L) {
  shown <- utils::capture.output(
    print(utils::head(tab, max), row.names = FALSE))
  cat("\n", paste0("    ", shown, "\n"), sep = "")
  if (nrow(tab) > max) {
    cat("    ... and ", nrow(tab) - max, " more strata, all in ",
        "`$allocation`\n", sep = "")
  }
}

#' @export
as.data.frame.drawn_plan <- function(x, row.names = NULL, optional = FALSE,
                                     ...) {
  fields <- c("n", "n_effective", "margin", "level", "N", "deff", "response",
              "target", "spread", "p", "capped", "short", "raised",
              "n_needed")
  out <- data.frame(unclass(x)[fields], stringsAsFactors = FALSE)
  if (!is.null(x$allocation)) {
    out$strata <- nrow(x$allocation)
    out$rule <- x$rule
    out$min_per_stratum <- x$min_per_stratum
  }
  if (!is.null(row.names)) row.names(out) <- row.names
  out
}

#' Design effect of a cluster sample, from the intraclass correlation
#'
#' The design effect of taking `m` rows from each cluster, for
#' [plan_size()] before a pilot has measured one with [deff()]:
#' `1 + (m - 1) * icc`, with `icc` the intraclass correlation of the
#' measure, the correlation between two rows of one cluster (Kish 1965).
#' Rows of one cluster tend to be alike, so each row taken from a cluster
#' already sampled adds less than a row from a new cluster would, and the
#' more alike they are and the more taken from each, the larger the design
#' effect.
#'
#' The formula assumes every cluster gives the same number of rows, drawn
#' by simple random sampling of clusters. `plan_size()` counts rows, so a
#' plan of `n` rows in clusters of `m` takes `ceiling(n / m)` clusters, and
#' the "Where the formula comes from" section of [plan_size()] explains why
#' a few clusters call for more than the plan.
#'
#' @param m Rows taken from each cluster, 1 or more. With one row a
#'   cluster, the design effect is 1 whatever `icc` is.
#' @param icc Intraclass correlation, from a pilot or a published study of
#'   the same measure, between `-1 / (m - 1)` and 1. A negative one, rows
#'   of a cluster less alike than rows taken at random, gives a design
#'   effect below 1.
#'
#' @return The design effect, a number. `m` and `icc` may be vectors,
#'   recycled to a common length, and give one design effect per pair.
#'
#' @references
#' Kish, L. (1965). *Survey Sampling*. Wiley.
#'
#' @examples
#' deff_cluster(m = 10, icc = 0.05)
#'
#' # Rows to draw, and the clusters of ten rows that hold them
#' p <- plan_size(margin = 5, sd = 40, N = 20000,
#'                deff = deff_cluster(m = 10, icc = 0.05))
#' p$n
#' ceiling(p$n / 10)
#'
#' # Takes of 5, 10 and 20 rows a cluster, side by side
#' plan_size(margin = 5, sd = 40, N = 20000,
#'           deff = deff_cluster(m = c(5, 10, 20), icc = 0.05))
#'
#' @seealso [plan_size()], [deff()] to measure the design effect of a past
#'   sample.
#' @export
deff_cluster <- function(m, icc) {
  if (!is.numeric(m) || !length(m) || anyNA(m) || !all(is.finite(m)) ||
      any(m < 1)) {
    stop("`m` must be the rows taken from each cluster, a finite number of ",
         "1 or more, not ", format_bad(m, value = TRUE), ".",
         call. = FALSE)
  }
  if (!is.numeric(icc) || !length(icc) || anyNA(icc) ||
      !all(is.finite(icc))) {
    stop("`icc` must be a finite number, the intraclass correlation, not ",
         format_bad(icc, value = TRUE), ".", call. = FALSE)
  }
  k <- max(length(m), length(icc))
  if (!all(c(length(m), length(icc)) %in% c(1L, k))) {
    stop("`m` has ", length(m), " values and `icc` has ", length(icc),
         ". Give each one value, or the same number of values.",
         call. = FALSE)
  }
  m <- rep_len(m, k)
  icc <- rep_len(icc, k)
  low <- ifelse(m > 1, -1 / (m - 1), -Inf)
  # A rounding error below the bound still means the bound.
  bad <- icc > 1 | icc < low - 1e-12
  if (any(bad)) {
    i <- which(bad)[1]
    stop("`icc` must lie between -1 / (m - 1) and 1, and is ",
         format(icc[i]), " with `m` = ", format(m[i]),
         if (icc[i] < 1) paste0(", below ", format(signif(low[i], 3))),
         ". A correlation above 1 is impossible, and one below the bound ",
         "would make the design effect negative.", call. = FALSE)
  }
  pmax(1 + (m - 1) * icc, 0)
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
#'     \item{`drawn`}{How the sample was drawn, from its `drawn_meta`
#'       attribute: `seed`, `rng_kind`, `version`, `time`, `frame_rows` and
#'       `made_by` (see "Seeding" in [draw()]). `NULL` for a sample from an
#'       earlier version of drawn, which recorded none. The print shows the
#'       seed, the version and the date.}
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
    group_col = groups$col,
    drawn = state$meta
  ), class = "drawn_summary")
}

#' The line of a summary that says how the sample was drawn
#'
#' The seed, or that there was none, with the generator when it is not R's
#' default; then the version of drawn and the date. A sample rebuilt by
#' as_drawn() says so, since its seed is not known.
#'
#' @noRd
meta_line <- function(meta) {
  stamp <- paste0("drawn ", meta$version, " on ",
                  format(meta$time, "%Y-%m-%d"))
  if (!identical(meta$made_by, "draw()")) {
    return(paste0("  rebuilt            by ", meta$made_by, " of ", stamp,
                  "; the seed is not recorded"))
  }
  kind <- if (!identical(meta$rng_kind, seed_kind())) {
    paste0(" (", paste(meta$rng_kind, collapse = ", "), ")")
  }
  paste0("  drawn              ",
         if (is.null(meta$seed)) {
           "without a seed"
         } else {
           paste0("with seed ", format(meta$seed, scientific = FALSE))
         }, kind, ", by ", stamp)
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
  if (!is.null(x$drawn)) cat(meta_line(x$drawn), "\n", sep = "")
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
