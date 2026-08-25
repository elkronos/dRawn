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
#'   maximises the required size.
#'
#' @section The corrections:
#' `N` applies a finite population correction: sampling 400 from 500 is very
#' different from 400 from 500,000, and past a certain point a bigger frame
#' stops mattering. `deff` inflates for the design — pass the [deff()] from a
#' comparable past sample, since a clustered design of 400 may carry the
#' information of 100. `response` inflates for non-response: at `0.6` you draw
#' enough to end up with what you need.
#'
#' @param margin Half-width of the confidence interval you want.
#' @param sd Standard deviation of the measure, for `target = "mean"` or
#'   `"total"`.
#' @param p Expected proportion, for `target = "proportion"`.
#' @param N Population size. `Inf` for an effectively unbounded frame.
#' @param level Confidence level.
#' @param deff Design effect to inflate by. `1` assumes simple random sampling.
#' @param response Expected response rate, between 0 and 1.
#' @param target `"mean"`, `"proportion"` or `"total"`. For `"total"` the
#'   margin is on the population total and `N` must be finite.
#'
#' @return A list with a `print()` method, holding:
#'   \describe{
#'     \item{`n`}{Draw this many. Never more than `N`.}
#'     \item{`n_effective`}{What you expect to analyse once `response` has
#'       taken its share.}
#'     \item{`capped`}{`TRUE` when the frame is not large enough to reach the
#'       margin at all, however you sample it.}
#'     \item{`short`}{`TRUE` when the frame could reach the margin but not at
#'       this `response` rate, so `n` is the whole frame and the margin
#'       achieved will be wider than the one asked for.}
#'     \item{`margin`, `level`, `N`, `deff`, `response`, `target`}{The inputs,
#'       returned so the assumptions travel with the number.}
#'     \item{`n_needed`}{Rows that would reach the margin if everyone
#'       responded. Equal to `n_effective` unless `short` is `TRUE`.}
#'     \item{`spread`}{The standard deviation used — `sd`, or
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
#' @seealso [deff()] to measure the design effect of a past sample, [draw()] to
#'   take the sample.
#' @export
plan_size <- function(margin, sd = NULL, p = 0.5, N = Inf, level = 0.95,
                      deff = 1, response = 1,
                      target = c("mean", "proportion", "total")) {
  target <- match.arg(target)
  check_pos <- function(x, arg, allow_inf = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 ||
        (!allow_inf && !is.finite(x))) {
      stop("`", arg, "` must be a single positive number",
           if (allow_inf) " (or Inf)" else "", ".", call. = FALSE)
    }
    x
  }
  margin <- check_pos(margin, "margin")
  deff <- check_pos(deff, "deff")
  N <- check_pos(N, "N", allow_inf = TRUE)
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
    if (!is.numeric(p) || length(p) != 1L || is.na(p) || p < 0 || p > 1) {
      stop("`p` must be a single number in [0, 1].", call. = FALSE)
    }
    spread <- sqrt(p * (1 - p))
    if (spread <= 0) {
      stop("`p` of exactly 0 or 1 implies no variation and no sample is ",
           "needed. Use a value strictly between 0 and 1.", call. = FALSE)
    }
    unit_margin <- margin
  } else {
    if (is.null(sd)) {
      stop("`sd` is required for target = \"", target, "\". Use a pilot ",
           "estimate, or a plausible range divided by four.", call. = FALSE)
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
  }

  z <- stats::qnorm(1 - (1 - level) / 2)
  n0 <- (z * spread / unit_margin)^2 * deff
  n_eff <- if (is.finite(N)) n0 / (1 + (n0 - 1) / N) else n0
  n_eff <- ceiling(n_eff)
  if (is.finite(N)) n_eff <- min(n_eff, N)
  n_draw <- ceiling(n_eff / response)

  # Two different situations, and only one of them is a census. `capped` means
  # the frame cannot reach the margin however hard you sample; `short` means it
  # could, but not once non-response has taken its share.
  capped <- is.finite(N) && n_eff >= N
  short <- is.finite(N) && n_draw > N
  n_needed <- n_eff
  if (short) {
    n_draw <- N
    n_eff <- min(n_eff, floor(n_draw * response))
  }

  structure(
    list(n = n_draw, n_effective = n_eff, margin = margin, level = level,
         N = N, deff = deff, response = response, target = target,
         spread = spread, capped = capped, short = short,
         n_needed = n_needed),
    class = "drawn_plan"
  )
}

#' @export
print.drawn_plan <- function(x, ...) {
  cat("Sample size for a ", x$target, "\n", sep = "")
  cat("  draw           ", fmt_n(x$n), "\n", sep = "")
  if (x$n != x$n_effective) {
    cat("  to analyse     ", fmt_n(x$n_effective), "  (after ",
        format(100 * x$response), "% response)\n", sep = "")
  }
  cat("  margin         +/- ", format(x$margin), " at ",
      format(100 * x$level), "% confidence\n", sep = "")
  cat("  assuming       sd ", signif(x$spread, 4),
      if (x$deff != 1) paste0(", deff ", x$deff) else "",
      ", N ", if (is.finite(x$N)) fmt_n(x$N) else "unbounded", "\n", sep = "")
  msg <- if (isTRUE(x$capped)) {
    paste0("That is a census: the frame is not large enough to reach this ",
           "margin by sampling, so every row is needed. ",
           if (x$response < 1) paste0(
             "Even then only about ", fmt_n(x$n_effective), " will come back, ",
             "which will fall short of the margin asked for. "),
           "Widen `margin`, or accept the precision a full count gives.")
  } else if (isTRUE(x$short)) {
    paste0("The frame can reach this margin -- ", fmt_n(x$n_needed),
           " rows would do it -- but not at ", format(100 * x$response),
           "% response, which would need more rows than the frame holds. ",
           "Every row is drawn, and the margin achieved will be wider than ",
           "the one asked for.")
  } else {
    NULL
  }
  if (!is.null(msg)) {
    cat("\n", strwrap(msg, prefix = "  "), sep = "\n")
    cat("\n")
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
#'     \item{`weight_range`}{The smallest and largest design weight.}
#'     \item{`weight_cv`}{Their coefficient of variation. Large values mean a
#'       few rows carry most of the estimate.}
#'     \item{`unreachable`}{Frame rows with inclusion probability 0 — the
#'       design could never have selected them. `NA` if the design has no
#'       closed-form inclusion probability.}
#'     \item{`by_group`}{A data frame of `group`, `drawn`, `in_frame` and
#'       `rate` per stratum or cluster, or `NULL` for a design with no
#'       grouping.}
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
  if (!is.data.frame(sample) || is.null(attr(sample, "drawn_design"))) {
    stop("`sample` must be a data frame from draw(..., weights = TRUE).",
         call. = FALSE)
  }
  design <- attr(sample, "drawn_design")
  pop <- attr(sample, "drawn_population")
  w <- sample$.weight

  grp_col <- switch(design_type(design),
    stratified = design$strata,
    cluster = ,
    multistage = design$clusters,
    NULL
  )
  by_group <- NULL
  if (!is.null(grp_col) && all(grp_col %in% names(sample))) {
    got <- table(group_key(sample, grp_col))
    have <- table(group_key(pop, grp_col))
    by_group <- data.frame(
      group = group_label(names(got)),
      drawn = as.integer(got),
      in_frame = as.integer(have[names(got)]),
      rate = round(as.integer(got) / as.integer(have[names(got)]), 4),
      row.names = NULL
    )
  }

  pi_all <- tryCatch(exact_inclusion(design, pop), error = function(e) NULL)

  structure(list(
    design = design_type(design),
    n = nrow(sample),
    N = nrow(pop),
    weight_range = if (length(w)) range(w, na.rm = TRUE) else c(NA, NA),
    weight_cv = if (length(w) > 1L && mean(w, na.rm = TRUE) > 0) {
      stats::sd(w, na.rm = TRUE) / mean(w, na.rm = TRUE)
    } else NA_real_,
    unreachable = if (is.null(pi_all)) NA_integer_ else sum(pi_all == 0),
    by_group = by_group,
    group_col = grp_col
  ), class = "drawn_summary")
}

#' @export
print.drawn_summary <- function(x, ...) {
  cat("Sample of ", fmt_n(x$n), " from ", fmt_n(x$N), "  (", x$design,
      " design)\n", sep = "")
  cat("  sampling fraction  ", signif(x$n / x$N, 3), "\n", sep = "")
  if (all(is.finite(x$weight_range))) {
    cat("  design weights     ", signif(x$weight_range[1], 4), " to ",
        signif(x$weight_range[2], 4), sep = "")
    if (is.finite(x$weight_cv)) {
      cat("   (cv ", signif(x$weight_cv, 3), ")", sep = "")
    }
    cat("\n")
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
          formatC(format(b$rate[i], nsmall = 3), width = 7), "\n", sep = "")
    }
  }
  invisible(x)
}
