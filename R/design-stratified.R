#' Stratified sampling
#'
#' Draws a sample from each stratum. Allocation uses the largest-remainder
#' method, so the returned row count matches `n` exactly rather than drifting
#' with per-stratum rounding.
#'
#' @param strata One or more column names defining the strata. Several columns
#'   are cross-classified.
#' @param n Total rows to draw across all strata, under any allocation. With
#'   counts in `allocation` it can be left out and is their sum; when given,
#'   it must equal that sum.
#' @param allocation How `n` is split across strata: a rule, or the count for
#'   each stratum. Of the rules, `"proportional"` gives each stratum a share
#'   of `n` in proportion to its size; `"equal"` splits `n` evenly;
#'   `"neyman"` gives shares proportional to `size * sd`, using the column
#'   named by `allocation_by` (Neyman 1934). Neyman minimises the variance of
#'   a total for a fixed `n` by putting more rows where the values vary most,
#'   and is the right choice when you have a frame variable correlated with
#'   what you are measuring. Under any rule, a stratum that would be
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
#' @param allocation_by Column whose within-stratum standard deviation drives
#'   `allocation = "neyman"`. Ignored otherwise. Every row allocated over needs
#'   a finite value.
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
#' @return A design object, for use with [draw()].
#'
#' @references
#' Neyman, J. (1934). On the two different aspects of the representative
#' method. *Journal of the Royal Statistical Society*, 97, 558–625.
#'
#' Cochran, W. G. (1977). *Sampling Techniques*, 3rd ed. Wiley.
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
#' @family designs
#' @seealso [draw()]
#' @export
design_stratified <- function(strata, n,
                              allocation = c("proportional", "equal", "neyman"),
                              allocation_by = NULL,
                              min_per_stratum = 0L, replace = FALSE,
                              na_rm = FALSE) {
  if (is.numeric(allocation) || is.data.frame(allocation)) {
    strata <- check_columns(strata, "strata")
    allocation <- check_stratum_counts(allocation, strata)
    total <- sum(as.numeric(stratum_counts(allocation)))
    if (missing(n) || is.null(n)) n <- total
    n <- check_count(n, "n")
    if (n != total) {
      stop("`n` (", n, ") is not the sum of the counts in `allocation` (",
           total, "). Leave `n` out to draw the counts as given.",
           call. = FALSE)
    }
    min_per_stratum <- check_count(min_per_stratum, "min_per_stratum")
    if (min_per_stratum > 0L) {
      stop("`min_per_stratum` (", min_per_stratum, ") is a floor for an ",
           "allocation rule, and `allocation` gives the counts. Give each ",
           "stratum a count of at least ", min_per_stratum, " there ",
           "instead.", call. = FALSE)
    }
  } else {
    if (!is.character(allocation)) {
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
  new_design("stratified", list(
    strata = strata,
    n = n,
    allocation = allocation,
    allocation_by = if (is.null(allocation_by)) NULL else
      check_columns(allocation_by, "allocation_by", 1L),
    min_per_stratum = check_count(min_per_stratum, "min_per_stratum"),
    replace = check_flag(replace, "replace"),
    na_rm = check_flag(na_rm, "na_rm")
  ))
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
  bad <- is.na(n) | !is.finite(n) | n < 0 | n != trunc(n) |
    n > .Machine$integer.max
  if (any(bad)) {
    h <- which(bad)[1]
    stop("The counts in `allocation` must be whole numbers of 0 or more; ",
         "stratum ", quote_labels(do.call(paste, c(
           lapply(spec[strata], function(v) as.character(v[h])),
           sep = "\r"))), " has ", format(n[h]), ".", call. = FALSE)
  }
  if (is.data.frame(allocation)) {
    spec$n <- as.integer(n)
    return(spec)
  }
  stats::setNames(as.integer(n), names(allocation))
}

#' The counts stored by a design built with counts, in the order given
#' @noRd
stratum_counts <- function(allocation) {
  if (is.data.frame(allocation)) allocation$n else unname(allocation)
}

#' Was the design given its per-stratum counts rather than a rule?
#' @noRd
has_counts <- function(design) !is.character(design$allocation)

#' Counts for a design's printout and title
#'
#' `c(north = 10, south = 5)` for up to six strata, with columns joined by
#' ` . ` and names that are not syntactic in backticks, and a count of the
#' strata beyond that.
#'
#' @noRd
fmt_counts <- function(allocation) {
  n <- stratum_counts(allocation)
  if (length(n) > 6L) return(paste0("<counts for ", length(n), " strata>"))
  lab <- if (is.data.frame(allocation)) {
    keys <- setdiff(names(allocation), "n")
    group_label(do.call(paste, c(lapply(allocation[keys], as.character),
                                 sep = "\r")))
  } else {
    group_label(names(allocation))
  }
  plain <- lab == make.names(lab)
  lab[!plain] <- paste0("`", lab[!plain], "`")
  paste0("c(", paste0(lab, " = ", n, collapse = ", "), ")")
}

#' Within-stratum spread for Neyman allocation
#'
#' One standard deviation per stratum, in the order of `idx_by_stratum`.
#'
#' @noRd
stratum_spread <- function(design, data, idx_by_stratum) {
  if (!identical(design$allocation, "neyman")) return(NULL)
  col <- design$allocation_by
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
  # Values near the largest double overflow inside sd(). The shares depend
  # only on the ratios of the spreads, so the same spreads on a smaller scale
  # give the same allocation.
  if (!all(is.finite(sp))) sp <- spread(max(abs(v[used])))
  sp
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
#' standard deviation under Neyman allocation, and `NULL` otherwise.
#'
#' @noRd
stratum_allocation <- function(design, data) {
  validate_data(data, required_columns = design$strata)
  check_key_columns(data, design$strata, "strata")
  key <- group_key(data, design$strata)
  bad <- check_na_policy(is.na(key), design$na_rm, "a missing stratum key")
  if (all(bad)) stop("Every row has a missing stratum key.", call. = FALSE)
  code <- as.integer(key)
  idx <- split(seq_len(nrow(data))[!bad], key[!bad])
  sizes <- lengths(idx, use.names = FALSE)
  spread <- NULL

  if (has_counts(design)) {
    n_alloc <- match_strata(design$allocation, data, design$strata,
                            "allocation", key = key)$n
    over <- which(n_alloc > sizes)
    if (!design$replace && length(over)) {
      h <- over[1]
      stop("Stratum `", group_label(levels(key)[h]), "` has ", sizes[h],
           " row(s) but `allocation` asks for ", n_alloc[h], ". Lower its ",
           "count, or use replace = TRUE.", call. = FALSE)
    }
  } else {
    check_draw_size(design$n, sum(!bad), design$replace)
    spread <- stratum_spread(design, data, idx)
    n_alloc <- unname(allocate(design$n, sizes, design$allocation,
                               design$min_per_stratum, cap = !design$replace,
                               spread = spread))
    if (!design$replace) {
      over <- which(n_alloc > sizes)
      if (length(over)) {
        h <- over[1]
        stop("Stratum `", group_label(levels(key)[h]), "` has ", sizes[h],
             " row(s) but ", n_alloc[h], " were requested. ",
             "Use replace = TRUE or lower `n`.", call. = FALSE)
      }
      if (sum(n_alloc) != design$n) {
        stop("`n` (", design$n, ") cannot be allocated across these ",
             length(sizes), " strata with min_per_stratum = ",
             design$min_per_stratum, "; the allocation totals ",
             sum(n_alloc), " row(s).", call. = FALSE)
      }
    }
  }

  list(code = code, labels = levels(key), idx = unname(idx), sizes = sizes,
       n_alloc = n_alloc, spread = spread)
}

#' @export
draw_design.drawn_design_stratified <- function(design, data) {
  al <- stratum_allocation(design, data)
  # Neyman allocation gives no share to a stratum whose `allocation_by` does
  # not vary, unless no stratum varies and it falls back to proportional.
  flat <- if (!is.null(al$spread) && sum(al$sizes * al$spread) > 0) {
    al$spread == 0
  }
  warn_empty_groups(
    al$n_alloc, al$sizes, al$labels, flat = flat, by = design$allocation_by,
    fix = if (has_counts(design)) {
      paste0("Give each stratum a count of at least 1 in `allocation` (or ",
             "2, so each stratum also supports a variance).")
    })
  reindex(data, take_within(al$idx, al$n_alloc, design$replace), sort = TRUE)
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
  p <- if (design$replace) {
    1 - (1 - 1 / al$sizes)^al$n_alloc
  } else {
    al$n_alloc / al$sizes
  }
  out <- numeric(length(al$code))
  keep <- !is.na(al$code)
  out[keep] <- p[al$code[keep]]
  out
}
