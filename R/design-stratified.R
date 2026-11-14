#' Stratified sampling
#'
#' Draws a sample from each stratum. Allocation uses the largest-remainder
#' method, so the returned row count matches `n` exactly rather than drifting
#' with per-stratum rounding.
#'
#' @param strata One or more column names defining the strata. Several columns
#'   are cross-classified.
#' @param n Total rows to draw across all strata, under any allocation, or a
#'   plan from [plan_size()]. With counts in `allocation` it can be left out
#'   and is their sum; when given, it must equal that sum.
#' @param allocation How `n` is split across strata: a rule, or the count for
#'   each stratum. Of the rules, `"proportional"` gives each stratum a share
#'   of `n` in proportion to its size; `"equal"` splits `n` evenly;
#'   `"neyman"` gives shares proportional to `size * sd`, with the standard
#'   deviation of each stratum from `allocation_by` (Neyman 1934). Neyman
#'   minimises the variance of a total for a fixed `n` by putting more rows
#'   where the values vary most, and is the right choice when you have a
#'   frame variable correlated with what you are measuring, or estimates of
#'   each stratum's standard deviation from a pilot or an earlier survey.
#'   Under any rule, a stratum that would be
#'   allocated more rows than it holds is taken whole, one that would get
#'   fewer than `min_per_stratum` gets that many, and the rest of `n` is
#'   re-split over the others by the same rule (Cochran 1977, section 5.9).
#'
#'   Counts are a named vector, such as `c(north = 10, south = 5, east = 5)`,
#'   when `strata` is one column, its names matched to that column's values
#'   as text; or a data frame holding the `strata` columns and a column `n`,
#'   matched on the values themselves. Every stratum of the frame needs a
#'   count, 0 included, and a count for a stratum the frame does not have is
#'   refused, so a misspelt name cannot quietly leave a stratum out. Without
#'   replacement no count may exceed the rows in its stratum. The design then
#'   draws exactly those counts, and the inclusion probabilities, variances
#'   and [as_svydesign()] follow them.
#' @param allocation_by The standard deviation of each stratum, for
#'   `allocation = "neyman"`. Ignored otherwise. Either a column of the
#'   frame, whose standard deviation within each stratum is used, and which
#'   needs a finite value on every row allocated over; or the standard
#'   deviations themselves, so the frame needs no such column. Those are a
#'   named vector, such as `c(north = 40, south = 15)`, when `strata` is one
#'   column, its names matched to that column's values as text; or a data
#'   frame holding the `strata` columns and a column `sd`, matched on the
#'   values themselves. Every stratum of the frame needs one, and each must
#'   be finite and 0 or more. A value for a stratum the frame does not have
#'   is ignored, so the same values serve a part of the frame, such as the
#'   `rest` of [design_certainty()] when every row of a stratum is above the
#'   threshold. A stratum with a standard deviation of 0 gets no share.
#' @param min_per_stratum Minimum rows from each stratum, under an allocation
#'   rule. With the default of
#'   `0`, a stratum smaller than about `N / n` can be allocated no rows at all;
#'   its rows then have inclusion probability 0 and a total estimated from the
#'   sample silently leaves them out, so [draw()] warns when that happens. `1`
#'   makes every stratum reachable, which is what an unbiased
#'   Horvitz-Thompson total needs, and `2` also lets every stratum contribute
#'   to the variance estimate. Over-sampling a small stratum this way does
#'   *not* bias an estimate made with the design weights, which correct for it;
#'   it only moves precision around. Without replacement a stratum smaller
#'   than the floor is taken whole, so `n` must cover `min_per_stratum` rows
#'   from each stratum, or all of a smaller one. With counts in `allocation`,
#'   give each stratum its count there instead.
#' @param replace Sample with replacement within each stratum?
#' @param na_rm Drop rows whose stratum key is `NA` instead of raising an error.
#' @param cost The cost of sampling one row in each stratum, for
#'   `allocation = "neyman"` only. In the same forms as standard deviations
#'   given in `allocation_by`, with the column of a data frame named `cost`,
#'   so one data frame of the `strata` columns, `sd` and `cost` can serve
#'   both. Each must be finite and above 0. The shares become proportional to
#'   `size * sd / sqrt(cost)`, which gives the smallest variance among the
#'   splits of the same total cost, `sum(cost * n_h)` over strata (Cochran
#'   1977, section 5.5). Before any stratum is held at its floor or its
#'   size, the `n` that spends a budget `C`, net of any fixed cost, is
#'   `C * sum(N_h * S_h / sqrt(c_h)) / sum(N_h * S_h * sqrt(c_h))`, rounded
#'   down. When no stratum has any spread, the shares are
#'   proportional to `size / sqrt(cost)`, as for equal spreads, and without
#'   `cost` to `size`.
#' @param within How rows are drawn inside each stratum: `NULL`, the
#'   default, for a simple random sample, or [design_systematic()] without
#'   `interval` or `n`, to walk each stratum in the order of its `order_by`
#'   column, or in frame order without one. Each stratum gets the rows the
#'   allocation gives it, at an even step from a random start of its own;
#'   see "Systematic sampling within strata". Needs `replace = FALSE`.
#'
#' @section Rounding, ties and the order of strata:
#' Each stratum first gets the whole part of its share of `n`. The rows left
#' over go one at a time to the strata with the largest fractional parts. When
#' two strata tie, the one with more rows not yet allocated goes first, and
#' after that the one that comes first in key order. `min_per_stratum` and
#' the size of each stratum bound the shares before they are rounded, so the
#' rounding never takes a stratum below its floor or past its size.
#'
#' Key order is the order of the levels for a factor column and the sorted
#' values for any other column, with text sorted byte by byte as in the C
#' locale (`"B"` before `"a"`). With several `strata` columns the first varies
#' slowest. Strata are also drawn from in this order, so a seed selects the
#' same rows in every locale. Renaming a stratum, or reordering a factor's
#' levels, can move a tied row from one stratum to another; the inclusion
#' probabilities always follow the allocation used.
#'
#' An empty string is a stratum like any other, shown as `(blank)` in
#' messages.
#'
#' @section Systematic sampling within strata:
#' With `within = design_systematic(order_by = "x")`, the strata are
#' allocated as without it, each stratum is sorted on `x`, and the `n_h`
#' rows allocated to a stratum of `N_h` are taken as
#' `design_systematic(n = n_h)` takes them from a frame of `N_h` rows: a
#' start drawn from (0, 1) and a step of `N_h / n_h` rows, which need not be
#' whole. Strata draw their starts in key order. Every row has probability
#' `n_h / N_h`, as under simple random sampling within strata, so the
#' allocation, [allocation_table()] and the weights are the same; what
#' changes is which rows can be drawn together. Two rows of one stratum have
#' the joint probability of the systematic design there, which is 0 for
#' rows less than a step apart, and two rows of different strata the product
#' of their probabilities. Rows whose `x` is missing are refused, or with
#' the `na_rm = TRUE` of the systematic design left out of their stratum
#' before it is allocated.
#'
#' Each stratum has one random start, so no design-unbiased variance exists,
#' as for [design_systematic()] on its own. [ht_total()] and [ht_mean()]
#' add up the strata's successive-difference approximations,
#' `(1 - n_h / N_h) * n_h / (2 * (n_h - 1)) * sum(diff(y / pi)^2)` over each
#' stratum's sampled rows in the order it was walked (Wolter 2007), with
#' the degrees of freedom of a stratified sample, the sampled rows less the
#' strata. A stratum with one sampled row from several is refused, as
#' without `within`. So is the jackknife: a stratum's rows all come from one
#' start, and deleting them one at a time does not mirror the design.
#' [as_svydesign()] hands `survey` the strata, whose simple-random variance
#' differs.
#'
#' Only [design_systematic()] is taken as `within`, with one start and no
#' `interval` or `n`: the allocation sets the size of each stratum. To walk
#' a stratum by several columns, give `order_by` a column of ranks, as in
#' the examples of [design_systematic()].
#'
#' @return A design object, for use with [draw()]. [allocation_table()] shows
#'   the rows it takes from each stratum of a frame before anything is drawn.
#'
#' @references
#' Neyman, J. (1934). On the two different aspects of the representative
#' method. *Journal of the Royal Statistical Society*, 97, 558–625.
#'
#' Cochran, W. G. (1977). *Sampling Techniques*, 3rd ed. Wiley.
#'
#' Wolter, K. M. (2007). *Introduction to Variance Estimation*, 2nd ed.
#' Springer.
#'
#' @examples
#' df <- data.frame(id = 1:100, site = rep(letters[1:4], each = 25))
#' table(draw(df, design_stratified("site", n = 20), seed = 1)$site)
#'
#' # Rare strata are covered only if you ask -- without the floor, the one
#' # "rare" row gets no allocation, and draw() warns that it is unreachable
#' skewed <- data.frame(id = 1:1000, g = c(rep("common", 999), "rare"))
#' draw(skewed, design_stratified("g", n = 10, min_per_stratum = 1), seed = 1)$g
#'
#' # A count for each stratum, drawn as given; `n` is their sum
#' fixed <- design_stratified("site", allocation = c(a = 8, b = 4, c = 4, d = 4))
#' table(draw(df, fixed, seed = 1)$site)
#'
#' # Neyman from standard deviations known from an earlier survey, with the
#' # cost of reaching a row in each stratum
#' pilot <- design_stratified("site", n = 20, allocation = "neyman",
#'                            allocation_by = c(a = 4, b = 1, c = 2, d = 8),
#'                            cost = c(a = 1, b = 1, c = 4, d = 4))
#' allocation_table(df, pilot)
#'
#' # Each stratum walked in order of `size`, from a start of its own
#' df$size <- (df$id * 37) %% 101
#' walked <- design_stratified("site", n = 20,
#'                             within = design_systematic(order_by = "size"))
#' draw(df, walked, seed = 1)
#'
#' @family designs
#' @seealso [draw()], [allocation_table()]
#' @export
design_stratified <- function(strata, n,
                              allocation = c("proportional", "equal", "neyman"),
                              allocation_by = NULL,
                              min_per_stratum = 0L, replace = FALSE,
                              na_rm = FALSE, cost = NULL, within = NULL) {
  if (!missing(n)) n <- plan_n(n)
  if (is.numeric(allocation) || is.data.frame(allocation)) {
    strata <- check_columns(strata, "strata")
    allocation <- check_stratum_counts(allocation, strata)
    n <- counts_n(if (!missing(n)) n, allocation)
    min_per_stratum <- check_count(min_per_stratum, "min_per_stratum")
    if (min_per_stratum > 0L) {
      stop("`min_per_stratum` (", min_per_stratum, ") is a floor for an ",
           "allocation rule, and `allocation` gives the counts. Give each ",
           "stratum a count of at least ", min_per_stratum, " there ",
           "instead.", call. = FALSE)
    }
  } else {
    if (!is.null(allocation) && !is.character(allocation)) {
      stop("`allocation` must be \"proportional\", \"equal\", \"neyman\" ",
           "or the count for each stratum, not ", format_bad(allocation), ".",
           call. = FALSE)
    }
    allocation <- check_choice(allocation,
                               c("proportional", "equal", "neyman"),
                               "allocation")
    if (allocation == "neyman" && is.null(allocation_by)) {
      stop("`allocation = \"neyman\"` needs `allocation_by`, the column ",
           "whose within-stratum spread should drive the split.",
           call. = FALSE)
    }
    strata <- check_columns(strata, "strata")
    if (missing(n) || is.null(n)) {
      stop("`n` is missing. Give the total number of rows to draw, or the ",
           "count for each stratum in `allocation`.", call. = FALSE)
    }
    if (is.numeric(n) && length(n) > 1L) {
      stop("`n` must be a single whole number, the total across strata, not ",
           format_bad(n), ". Give the count for each stratum as ",
           "`allocation`, as in allocation = c(north = 10, south = 5).",
           call. = FALSE)
    }
    n <- check_count(n, "n")
  }
  if (!is.null(cost) && !identical(allocation, "neyman")) {
    stop("`cost` applies only to allocation = \"neyman\", and ",
         if (is.character(allocation)) {
           paste0("`allocation` is \"", allocation, "\". Use allocation = ",
                  "\"neyman\" with `allocation_by`, or leave `cost` out.")
         } else {
           "`allocation` gives the count for each stratum. Leave `cost` out."
         }, call. = FALSE)
  }
  new_design("stratified", list(
    strata = strata,
    n = n,
    allocation = allocation,
    allocation_by = check_allocation_by(allocation_by, strata),
    cost = if (!is.null(cost)) {
      check_stratum_values(cost, strata, "cost", "cost", positive = TRUE)
    },
    min_per_stratum = check_count(min_per_stratum, "min_per_stratum"),
    replace = check_flag(replace, "replace"),
    na_rm = check_flag(na_rm, "na_rm"),
    within = check_within(within, replace)
  ))
}

#' Check `within`: `NULL`, or a systematic design the allocation can size
#' @noRd
check_within <- function(within, replace) {
  if (is.null(within)) return(NULL)
  if (!inherits(within, "drawn_design_systematic")) {
    stop("`within` must be NULL, for a simple random sample in each ",
         "stratum, or design_systematic(), to walk each stratum in order, ",
         "not ", if (is_design(within)) design_label(within) else
           format_bad(within), ".", call. = FALSE)
  }
  given <- c("interval", "n")[c(!is.null(within$interval),
                                 !is.null(within$n))]
  if (length(given)) {
    stop("`within` takes design_systematic() without `interval` or `n`, ",
         "and has `", given[1], "` (", within[[given[1]]], "). The ",
         "allocation sets the size of each stratum. Leave it out, as in ",
         "within = design_systematic(order_by = \"x\").", call. = FALSE)
  }
  if (isTRUE(replace)) {
    stop("`within` cannot be combined with replace = TRUE: a systematic ",
         "walk takes each row of a stratum at most once. Use replace = ",
         "FALSE.", call. = FALSE)
  }
  within
}

#' Does the design walk each stratum systematically?
#' @noRd
walks_strata <- function(design) {
  inherits(design, "drawn_design_stratified") && !is.null(design$within)
}

#' Check `allocation_by`: a column name, or a standard deviation per stratum
#'
#' Returns the column name, or the standard deviations in the form the
#' design stores per-stratum values.
#'
#' @noRd
check_allocation_by <- function(allocation_by, strata) {
  if (is.null(allocation_by)) return(NULL)
  if (is.character(allocation_by)) {
    return(check_columns(allocation_by, "allocation_by", 1L))
  }
  if (!is.numeric(allocation_by) && !is.data.frame(allocation_by)) {
    stop("`allocation_by` must be a column name, or the standard deviation ",
         "of each stratum as a named vector or a data frame of the `strata` ",
         "columns and `sd`, not ", format_bad(allocation_by), ".",
         call. = FALSE)
  }
  check_stratum_values(allocation_by, strata, "allocation_by", "sd")
}

#' Check numbers given per stratum, such as standard deviations or costs
#'
#' Each must be finite and 0 or more, or above 0 with `positive`. Returns a
#' named double vector, or a plain data frame of the `strata` columns and the
#' `value` column. Whether every stratum has one is a question about the
#' frame, answered by `match_strata()` when the design meets one.
#'
#' @noRd
check_stratum_values <- function(spec, strata, arg, value,
                                 positive = FALSE) {
  what <- c(sd = "standard deviations", cost = "costs")[[value]]
  out <- strata_frame(spec, strata, arg, value)
  v <- out[[value]]
  if (!is.numeric(v)) {
    stop("The ", what, " in `", arg, "` must be numbers, not ",
         class(v)[1], ".", call. = FALSE)
  }
  bad <- !is.finite(v) | if (positive) v <= 0 else v < 0
  if (any(bad)) {
    h <- which(bad)[1]
    stop("The ", what, " in `", arg, "` must be finite numbers ",
         if (positive) "above 0" else "of 0 or more", ", and stratum ",
         spec_label(out, strata, h), " has ", format(v[h]), ". Give every ",
         "stratum ", if (positive) {
           "the cost of sampling one row there."
         } else {
           "an estimate, such as one from a pilot or an earlier survey."
         }, call. = FALSE)
  }
  if (is.data.frame(spec)) {
    out[[value]] <- as.numeric(v)
    return(out)
  }
  stats::setNames(as.numeric(v), names(spec))
}

#' One row of per-stratum values given by the user, named for a message
#' @noRd
spec_label <- function(spec, strata, h) {
  quote_labels(do.call(paste, c(
    lapply(spec[strata], function(v) as.character(v[h])), sep = "\r")))
}

#' Check the per-stratum counts given as `allocation`
#'
#' Returns them in the form the design stores: a named integer vector, or a
#' plain data frame of the `strata` columns and an integer `n`. Whether every
#' stratum has a count is a question about the frame, answered by
#' `match_strata()` when the design meets one.
#'
#' @noRd
check_stratum_counts <- function(allocation, strata) {
  spec <- strata_frame(allocation, strata, "allocation", "n")
  n <- spec$n
  if (!is.numeric(n)) {
    stop("The counts in `allocation` must be numbers, not ",
         class(n)[1], ".", call. = FALSE)
  }
  bad <- is.na(n) | !is.finite(n) | n < 0 | n != trunc(n)
  if (any(bad)) {
    h <- which(bad)[1]
    stop("The counts in `allocation` must be whole numbers of 0 or more; ",
         "stratum ", spec_label(spec, strata, h), " has ", format(n[h]), ".",
         call. = FALSE)
  }
  big <- n > .Machine$integer.max
  if (any(big)) {
    h <- which(big)[1]
    stop("The counts in `allocation` must be at most ", .Machine$integer.max,
         ", the largest integer R holds; stratum ",
         spec_label(spec, strata, h), " has ", format(n[h]), ". Give ",
         "that stratum a smaller count.", call. = FALSE)
  }
  if (is.data.frame(allocation)) {
    spec$n <- as.integer(n)
    return(spec)
  }
  stats::setNames(as.integer(n), names(allocation))
}

#' `n` for a design given per-stratum counts
#'
#' Their sum when `n` is `NULL`, and otherwise `n`, which must equal it.
#'
#' @noRd
counts_n <- function(n, allocation) {
  total <- sum(as.numeric(stratum_counts(allocation)))
  n <- check_count(n %||% total, "n")
  if (n != total) {
    stop("`n` (", n, ") is not the sum of the counts in `allocation` (",
         total, "). Leave `n` out to draw the counts as given.",
         call. = FALSE)
  }
  n
}

#' The counts stored by a design built with counts, in the order given
#' @noRd
stratum_counts <- function(allocation) {
  if (is.data.frame(allocation)) allocation$n else unname(allocation)
}

#' Was the design given its per-stratum counts rather than a rule?
#' @noRd
has_counts <- function(design) {
  !is.null(design$allocation) && !is.character(design$allocation)
}

#' Per-stratum values for a design's printout and title
#'
#' Counts, standard deviations or costs, as the design stores them: a named
#' vector, or a data frame whose last column holds the values.
#' `c(north = 10, south = 5)` for up to six strata, with columns joined by
#' ` . ` and names that are not syntactic in backticks, and a count of the
#' strata beyond that.
#'
#' @noRd
fmt_counts <- function(x) {
  v <- if (is.data.frame(x)) x[[ncol(x)]] else unname(x)
  what <- if (is.integer(v)) "counts" else "values"
  if (length(v) > 6L) {
    return(paste0("<", what, " for ", length(v), " strata>"))
  }
  lab <- if (is.data.frame(x)) {
    keys <- names(x)[-ncol(x)]
    group_label(do.call(paste, c(lapply(x[keys], as.character),
                                 sep = "\r")))
  } else {
    group_label(names(x))
  }
  plain <- lab == make.names(lab)
  lab[!plain] <- paste0("`", lab[!plain], "`")
  shown <- if (is.integer(v)) {
    as.character(v)
  } else {
    trimws(formatC(v, format = "fg", digits = 6))
  }
  paste0("c(", paste0(lab, " = ", shown, collapse = ", "), ")")
}

#' Within-stratum spread for Neyman allocation
#'
#' One standard deviation per stratum, in the order of `idx_by_stratum` and
#' the levels of `key`: given in the design, or worked out from the column
#' `allocation_by`. A spread worked out on a smaller scale, because the
#' column's own overflows, carries that scale as the attribute `scale`.
#'
#' @noRd
stratum_spread <- function(design, data, idx_by_stratum,
                           key = group_key(data, design$strata)) {
  if (!identical(design$allocation, "neyman")) return(NULL)
  col <- design$allocation_by
  if (!is.character(col)) {
    return(match_strata(col, data, design$strata, "allocation_by",
                        value = "sd", key = key, ignore_extra = TRUE)$sd)
  }
  validate_data(data, required_columns = col)
  v <- data[[col]]
  if (!is.numeric(v)) {
    stop("`allocation_by` names `", col, "`, which must be numeric for Neyman ",
         "allocation, not ", class(v)[1], ".", call. = FALSE)
  }
  # Only the rows that survive the design's own NA handling matter. Checking
  # the whole frame made a design that draws fine fail on inclusion_prob().
  # sd() of a stratum holding an infinite value is NaN, which would read as
  # no spread and allocate that stratum nothing.
  used <- unlist(idx_by_stratum, use.names = FALSE)
  bad <- !is.finite(v[used])
  if (any(bad)) {
    stop(sum(bad), " value(s) of `", col, "` are missing or non-finite; ",
         "Neyman allocation needs a finite value for every row it allocates ",
         "over. Drop those rows, or set them to a finite value.",
         call. = FALSE)
  }
  spread <- function(scale) {
    vapply(idx_by_stratum, function(i) {
      if (length(i) < 2L) 0 else stats::sd(v[i] / scale)
    }, numeric(1), USE.NAMES = FALSE)
  }
  sp <- spread(1)
  # Values near the largest double overflow inside sd(), which squares
  # them. In the stratum's own largest value as the unit they do not, and
  # the other strata keep the spreads sd() gives them.
  big <- !is.finite(sp)
  if (any(big)) {
    sp[big] <- vapply(idx_by_stratum[big], function(i) {
      unit <- max(abs(v[i]))
      unit * stats::sd(v[i] / unit)
    }, numeric(1), USE.NAMES = FALSE)
  }
  # A spread beyond the largest double needs a smaller scale for all of
  # them. The shares depend only on the ratios of the spreads, so the same
  # spreads on a smaller scale give the same allocation.
  if (!all(is.finite(sp))) {
    scale <- max(abs(v[used]))
    sp <- structure(spread(scale), scale = scale)
  }
  sp
}

#' The cost of one row in each stratum, in the order of the levels of `key`
#' @noRd
stratum_cost <- function(design, data, key = group_key(data, design$strata)) {
  if (is.null(design$cost)) return(NULL)
  match_strata(design$cost, data, design$strata, "cost", value = "cost",
               key = key, ignore_extra = TRUE)$cost
}

#' Strata, their sizes and their allocation
#'
#' The one place a stratified design's strata and per-stratum sample sizes are
#' worked out, so the draw, the inclusion probabilities and the joint
#' probabilities cannot disagree. Everything is by position: `code` is each
#' row's stratum as a position in `sizes`, `n_alloc` and `idx`, and `NA` for a
#' row with a missing key.
#'
#' It is also where an allocation the frame cannot supply is refused, on every
#' path, so [inclusion_prob()] does not report a census for an `n` larger than
#' the frame that [draw()] refuses. Counts given in `allocation` are matched
#' to the strata here and used as they are. `spread` is each stratum's
#' standard deviation under Neyman allocation, and `NULL` otherwise; `cost`
#' is each stratum's cost per row when the design has costs, and `NULL`
#' otherwise.
#'
#' [design_spread()] with `strata` splits its `n` here too. It has no
#' `replace`, `min_per_stratum`, `allocation_by` or `cost`, and is treated as
#' a design without replacement, with no floor.
#'
#' Under `within`, a row whose `order_by` is missing is refused, or with the
#' systematic design's `na_rm` left out like a row with a missing stratum
#' key.
#'
#' @noRd
stratum_allocation <- function(design, data) {
  validate_data(data, required_columns = design$strata)
  check_key_columns(data, design$strata, "strata")
  key <- group_key(data, design$strata)
  bad <- check_na_policy(is.na(key), design$na_rm, "a missing stratum key")
  if (all(bad)) stop("Every row has a missing stratum key.", call. = FALSE)
  if (walks_strata(design)) {
    # A row with no place in its stratum's walk cannot be drawn, so it is
    # left out before the strata are counted and allocated.
    gone <- walk_missing(design$within, data, bad, key)
    if (any(gone)) {
      bad <- bad | gone
      key[gone] <- NA
    }
  }
  code <- as.integer(key)
  idx <- split(seq_len(nrow(data))[!bad], key[!bad])
  sizes <- lengths(idx, use.names = FALSE)
  spread <- NULL
  cost <- NULL

  replace <- isTRUE(design$replace)
  if (has_counts(design)) {
    n_alloc <- match_strata(design$allocation, data, design$strata,
                            "allocation", key = key)$n
    over <- which(n_alloc > sizes)
    if (!replace && length(over)) {
      h <- over[1]
      stop("Stratum `", group_label(levels(key)[h]), "` has ", sizes[h],
           " row(s) but `allocation` asks for ", n_alloc[h], ". Lower its ",
           "count", if (!is.null(design$replace)) ", or use replace = TRUE",
           ".", call. = FALSE)
    }
  } else {
    check_draw_size(design$n, sum(!bad), replace)
    spread <- stratum_spread(design, data, idx, key)
    cost <- stratum_cost(design, data, key)
    n_alloc <- unname(allocate(design$n, sizes, design$allocation,
                               design$min_per_stratum %||% 0L, cap = !replace,
                               spread = spread, cost = cost))
    if (!replace) {
      over <- which(n_alloc > sizes)
      if (length(over)) {
        h <- over[1]
        stop("Stratum `", group_label(levels(key)[h]), "` has ", sizes[h],
             " row(s) but ", n_alloc[h], " were requested. ",
             "Use replace = TRUE or lower `n`.", call. = FALSE)
      }
      if (sum(n_alloc) != design$n) {
        stop("`n` (", design$n, ") cannot be allocated ",
             if (length(sizes) == 1L) {
               "in the one stratum"
             } else {
               paste0("across these ", length(sizes), " strata")
             }, " with min_per_stratum = ",
             design$min_per_stratum, "; the allocation totals ",
             sum(n_alloc), " row(s).", call. = FALSE)
      }
    }
  }

  list(code = code, labels = levels(key), idx = unname(idx), sizes = sizes,
       n_alloc = n_alloc, spread = spread, cost = cost)
}

#' How a design splits its sample over the strata of a frame
#'
#' The rows [draw()] takes from each stratum, the inclusion probability that
#' gives each row there, and the reason for any stratum held at the floor,
#' capped at its size or left with no rows, all before anything is drawn.
#' Under Neyman allocation it also shows the standard deviation of each
#' stratum the split used, and the cost per row when the design has costs.
#' The numbers come from the same computation as [draw()] and
#' [inclusion_prob()], so they always agree with both.
#'
#' @param data A data frame.
#' @param design [design_stratified()], [design_spread()] or
#'   [design_weighted()] with `strata`, or [design_certainty()] whose `rest`
#'   is one of these. For the last the
#'   table describes `rest` over the rows below the threshold, and the rows
#'   at or above it, which are taken with certainty, are not counted in any
#'   stratum. Other designs are refused.
#'
#' @return A data frame with one row per stratum of `data`, in key order (see
#'   [design_stratified()]): the `strata` columns, then
#'   \describe{
#'     \item{`N`}{Rows in the stratum.}
#'     \item{`n`}{Rows drawn from it.}
#'     \item{`prob`}{The inclusion probability of each of its rows, as
#'       [inclusion_prob()] gives it. `NA` for [design_spread()] with `size`
#'       and for [design_weighted()], whose rows in one stratum have
#'       different probabilities; there [inclusion_prob()] gives them row by
#'       row.}
#'     \item{`sd`}{Under Neyman allocation, the stratum's standard deviation
#'       that the split used: of the `allocation_by` column within the
#'       stratum, or as given. `NA` under other allocations.}
#'     \item{`cost`}{The cost of one row, when the design has `cost`, and
#'       `NA` otherwise.}
#'     \item{`status`}{Why the stratum got `n`: `"rule"` for its share under
#'       the allocation rule, rounded; `"count"` for a count given in
#'       `allocation`; `"floor"` for a stratum raised to `min_per_stratum`
#'       from a smaller share; `"capped"` for one the rule asked more rows of
#'       than it holds, so it is taken whole and the rows it could not give
#'       go to the other strata; `"take-all"` for one taken whole for another
#'       reason, such as a size below `min_per_stratum` or a count equal to
#'       its size; and `"empty"` for one allocated no rows, whose rows can
#'       never be drawn. With `replace = TRUE` no stratum is capped or taken
#'       whole.}
#'   }
#'
#' @examples
#' df <- data.frame(id = 1:300,
#'                  site = rep(c("a", "b", "c", "d"), c(150, 100, 45, 5)),
#'                  v = c(1:150 %% 7, 1:100 %% 19, 1:45 %% 3, 1:5 * 40))
#' d <- design_stratified("site", n = 30, allocation = "neyman",
#'                        allocation_by = "v", min_per_stratum = 2)
#' allocation_table(df, d)
#'
#' # The same counts as a draw, and the same probabilities as inclusion_prob()
#' table(draw(df, d, seed = 1)$site)
#' tapply(inclusion_prob(df, d), df$site, unique)
#'
#' @seealso [design_stratified()], [design_spread()], [inclusion_prob()]
#' @export
allocation_table <- function(data, design) {
  if (!is_design(design)) {
    stop("`design` must come from one of the design_*() constructors, not ",
         a_class(design), ". See ?designs.", call. = FALSE)
  }
  validate_data(data)
  if (inherits(design, "drawn_design_certainty") &&
      has_strata(design$rest)) {
    sp <- certainty_split(design, data)
    validate_data(data, required_columns = design$rest$strata)
    below <- sp$data[sp$rest, , drop = FALSE]
    return(in_rest(design, length(sp$rest),
                   stratum_table(design$rest, below)))
  }
  if (!has_strata(design)) {
    stop("`design` is ", design_label(design), ", which has no strata to ",
         "allocate over. allocation_table() takes design_stratified(), ",
         "design_spread() or design_weighted() with `strata`, or ",
         "design_certainty() with one of those as `rest`.", call. = FALSE)
  }
  stratum_table(design, data)
}

#' Does the design split `n` over strata?
#' @noRd
has_strata <- function(design) {
  inherits(design, "drawn_design_stratified") ||
    (inherits(design, c("drawn_design_spread", "drawn_design_weighted")) &&
       !is.null(design$strata))
}

#' The table `allocation_table()` returns, for a design with strata
#'
#' A frame with no rows gives a table with no strata. Only a certainty
#' design with no rows below its threshold, where `rest` draws nothing,
#' passes one. The strata of a spread or weighted design hold the rows it
#' can draw, and have no floor.
#'
#' @noRd
stratum_table <- function(design, data) {
  cols <- c("N", "n", "prob", "sd", "cost", "status")
  clash <- intersect(design$strata, cols)
  if (length(clash)) {
    stop("`strata` column ", format_bad(clash[1]), " has the name of a ",
         "column allocation_table() adds. Rename it in `data` and in the ",
         "design.", call. = FALSE)
  }
  if (!nrow(data)) {
    keys <- lapply(design$strata, function(col) data[[col]])
    return(table_frame(keys, design$strata, list(
      N = integer(0), n = integer(0), prob = numeric(0), sd = numeric(0),
      cost = numeric(0), status = character(0))))
  }
  if (inherits(design, "drawn_design_spread")) {
    fr <- spread_frame(design, data)
    al <- list(idx = lapply(fr$idx, function(i) fr$keep[i]),
               sizes = fr$sizes, n_alloc = fr$n_alloc)
    prob <- if (is.null(design$size)) {
      fr$n_alloc / fr$sizes
    } else {
      rep(NA_real_, length(fr$sizes))
    }
  } else if (inherits(design, "drawn_design_weighted")) {
    st <- weighted_setup(design, data)
    al <- list(idx = lapply(st$idx, function(i) st$keep[i]),
               sizes = st$sizes, n_alloc = st$n_alloc)
    # Rows of one stratum have probabilities proportional to their size.
    prob <- rep(NA_real_, length(st$sizes))
  } else {
    al <- stratum_allocation(design, data)
    prob <- stratum_prob(design, al)
  }
  k <- length(al$sizes)
  first <- vapply(al$idx, function(i) i[1L], integer(1))
  keys <- lapply(design$strata, function(col) data[[col]][first])
  cap <- !isTRUE(design$replace)
  status <- if (has_counts(design)) {
    out <- rep("count", k)
    if (cap) out[al$n_alloc == al$sizes] <- "take-all"
    out[al$n_alloc == 0L] <- "empty"
    out
  } else {
    allocation_status(
      design$n, al$n_alloc, al$sizes,
      share = allocation_share(design$allocation, al$sizes, al$spread,
                               al$cost),
      lower = group_floor(design$min_per_stratum %||% 0L, al$sizes, cap),
      cap = cap)
  }
  sd <- if (is.null(al$spread)) {
    rep(NA_real_, k)
  } else {
    as.numeric(al$spread) * (attr(al$spread, "scale") %||% 1)
  }
  table_frame(keys, design$strata, list(
    N = as.integer(al$sizes),
    n = as.integer(al$n_alloc),
    prob = prob,
    sd = sd,
    cost = if (is.null(al$cost)) rep(NA_real_, k) else as.numeric(al$cost),
    status = status))
}

#' A plain data frame of key columns followed by value columns
#' @noRd
table_frame <- function(keys, strata, values) {
  out <- c(stats::setNames(keys, strata), values)
  structure(out, class = "data.frame", row.names = seq_along(out[[1L]]))
}

#' @export
draw_design.drawn_design_stratified <- function(design, data) {
  al <- stratum_allocation(design, data)
  # Neyman allocation gives no share to a stratum whose `allocation_by` does
  # not vary, unless no stratum varies and it falls back to proportional.
  flat <- if (!is.null(al$spread) && sum(al$sizes * al$spread) > 0) {
    al$spread == 0
  }
  by <- design$allocation_by
  warn_empty_groups(
    al$n_alloc, al$sizes, al$labels, flat = flat,
    by = if (is.character(by)) by,
    fix = if (has_counts(design)) {
      paste0("Give each stratum a count of at least 1 in `allocation` (or ",
             "2, so each stratum also supports a variance).")
    })
  if (walks_strata(design)) {
    return(reindex(data, take_walks(design, data, al), sort = TRUE))
  }
  reindex(data, take_within(al$idx, al$n_alloc, design$replace), sort = TRUE)
}

# ---- systematic sampling within strata --------------------------------------

#' Rows `within` cannot place on its walk
#'
#' A missing `order_by` among the rows that have a stratum, `bad` marking the
#' others and `key` giving each row's stratum. Refused unless the systematic
#' design has `na_rm = TRUE`, and refused even then for a stratum that would
#' be left with no rows to walk.
#'
#' @noRd
walk_missing <- function(within, data, bad, key) {
  col <- within$order_by
  if (is.null(col)) return(rep(FALSE, nrow(data)))
  if (!col %in% names(data)) {
    stop("`data` is missing column: `", col, "`.", call. = FALSE)
  }
  check_key_columns(data, col, "order_by")
  what <- paste0("a missing `", col, "`")
  gone <- check_na_policy(is.na(data[[col]]) & !bad, within$na_rm, what)
  code <- as.integer(key)
  k <- nlevels(key)
  emptied <- tabulate(code[!bad], k) > 0L &
    tabulate(code[!bad & !gone], k) == 0L
  if (any(emptied)) {
    one <- sum(emptied) == 1L
    stop(if (one) "Stratum " else "Strata ",
         quote_labels(levels(key)[emptied]), if (one) " has " else " have ",
         what, " in every row, so `within` has nothing to walk there. Give ",
         "those rows a `", col, "`, or leave ", if (one) "the stratum" else
           "those strata", " out of `data`.", call. = FALSE)
  }
  gone
}

#' The frame rows of every stratum, in the order the design walks them
#'
#' Stratum by stratum in key order, each sorted on the `order_by` of
#' `within`, or left in frame order without one. Text is sorted byte by byte
#' and ties keep frame order, as in [design_systematic()]. `code` is each
#' row's stratum from `stratum_allocation()`, `NA` for a row in none.
#'
#' @noRd
stratum_walk <- function(design, data, code) {
  kept <- which(!is.na(code))
  col <- design$within$order_by
  if (!is.null(col)) kept <- kept[key_order(data[[col]][kept])]
  # Radix ordering is stable, so each stratum keeps the order of its key.
  kept[order(code[kept], method = "radix")]
}

#' Each frame row's position in its stratum's walk, `NA` for a row in none
#'
#' `k` is the number of strata that `code` indexes.
#'
#' @noRd
walk_positions <- function(design, data, code, k) {
  walk <- stratum_walk(design, data, code)
  before <- c(0L, cumsum(tabulate(code, k)))
  pos <- rep(NA_integer_, length(code))
  pos[walk] <- seq_along(walk) - before[code[walk]]
  pos
}

#' The rows a stratified design with `within` draws
#'
#' Stratum `h`, in key order, takes the walk positions `fractional_walk()`
#' gives for `n_h` of `N_h` rows, from one `stats::runif()` of its own. A
#' stratum allocated no rows draws no start.
#'
#' @noRd
take_walks <- function(design, data, al) {
  walk <- stratum_walk(design, data, al$code)
  before <- c(0L, cumsum(al$sizes))
  unlist(lapply(which(al$n_alloc > 0L), function(h) {
    walk[before[h] + fractional_walk(stats::runif(1L), al$n_alloc[h],
                                     al$sizes[h])]
  }))
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_stratified <- function(design, data) {
  stratum_pi(design, stratum_allocation(design, data))
}

#' Each row's inclusion probability, from `stratum_allocation()`
#'
#' Linear in the rows, whatever the number of strata. A row with no stratum
#' gets 0.
#'
#' @noRd
stratum_pi <- function(design, al) {
  p <- stratum_prob(design, al)
  out <- numeric(length(al$code))
  keep <- !is.na(al$code)
  out[keep] <- p[al$code[keep]]
  out
}

#' The inclusion probability of a row in each stratum of `stratum_allocation()`
#' @noRd
stratum_prob <- function(design, al) {
  if (design$replace) {
    1 - (1 - 1 / al$sizes)^al$n_alloc
  } else {
    al$n_alloc / al$sizes
  }
}
