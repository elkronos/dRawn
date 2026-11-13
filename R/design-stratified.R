#' Stratified sampling
#'
#' Draws a sample from each stratum. Allocation uses the largest-remainder
#' method, so the returned row count matches `n` exactly rather than drifting
#' with per-stratum rounding.
#'
#' @param strata One or more column names defining the strata. Several columns
#'   are cross-classified.
#' @param n Total rows to draw across all strata, under either allocation.
#' @param allocation How `n` is split across strata. `"proportional"` gives each
#'   stratum a share of `n` in proportion to its size; `"equal"` splits `n`
#'   evenly; `"neyman"` gives shares proportional to `size * sd`, using the
#'   column named by `allocation_by` (Neyman 1934). Neyman minimises the
#'   variance of a total for a fixed `n` by putting more rows where the values
#'   vary most, and is the right choice when you have a frame variable
#'   correlated with what you are measuring. Under any rule, a stratum that
#'   would be allocated more rows than it holds is taken whole and the rest of
#'   `n` re-split over the others by the same rule (Cochran 1977, section 5.9).
#' @param allocation_by Column whose within-stratum standard deviation drives
#'   `allocation = "neyman"`. Ignored otherwise.
#' @param min_per_stratum Minimum rows from each stratum. With the default of
#'   `0`, a stratum smaller than about `N / n` can be allocated no rows at all;
#'   its rows then have inclusion probability 0 and a total estimated from the
#'   sample silently leaves them out, so [draw()] warns when that happens. `1`
#'   makes every stratum reachable, which is what an unbiased
#'   Horvitz-Thompson total needs, and `2` also lets every stratum contribute
#'   to the variance estimate. Over-sampling a small stratum this way does
#'   *not* bias an estimate made with the design weights, which correct for it;
#'   it only moves precision around.
#' @param replace Sample with replacement within each stratum?
#' @param na_rm Drop rows whose stratum key is `NA` instead of raising an error.
#'
#' @section Rounding, ties and the order of strata:
#' Each stratum first gets the whole part of its share of `n`. The rows left
#' over go one at a time to the strata with the largest fractional parts. When
#' two strata tie, the one with more rows not yet allocated goes first, and
#' after that the one that comes first in key order. When `min_per_stratum`
#' pushes the total past `n`, rows are taken back from the strata with the
#' smallest fractional parts first, ties again going by key order.
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
#' @family designs
#' @seealso [draw()]
#' @export
design_stratified <- function(strata, n,
                              allocation = c("proportional", "equal", "neyman"),
                              allocation_by = NULL,
                              min_per_stratum = 0L, replace = FALSE,
                              na_rm = FALSE) {
  allocation <- check_choice(allocation,
                             c("proportional", "equal", "neyman"),
                             "allocation")
  if (allocation == "neyman" && is.null(allocation_by)) {
    stop("`allocation = \"neyman\"` needs `allocation_by`, the column whose ",
         "within-stratum spread should drive the split.", call. = FALSE)
  }
  new_design("stratified", list(
    strata = check_columns(strata, "strata"),
    n = check_count(n, "n"),
    allocation = allocation,
    allocation_by = if (is.null(allocation_by)) NULL else
      check_columns(allocation_by, "allocation_by", 1L),
    min_per_stratum = check_count(min_per_stratum, "min_per_stratum"),
    replace = check_flag(replace, "replace"),
    na_rm = check_flag(na_rm, "na_rm")
  ))
}


#' Within-stratum spread for Neyman allocation
#'
#' One standard deviation per stratum, in the order of `idx_by_stratum`.
#'
#' @noRd
stratum_spread <- function(design, data, idx_by_stratum) {
  if (design$allocation != "neyman") return(NULL)
  col <- design$allocation_by
  validate_data(data, required_columns = col)
  v <- data[[col]]
  if (!is.numeric(v)) {
    stop("`allocation_by` names `", col, "`, which must be numeric for Neyman ",
         "allocation, not ", class(v)[1], ".", call. = FALSE)
  }
  # Only the rows that survive the design's own NA handling matter. Checking
  # the whole frame made a design that draws fine fail on inclusion_prob().
  used <- unlist(idx_by_stratum, use.names = FALSE)
  if (anyNA(v[used])) {
    stop(sum(is.na(v[used])), " value(s) of `", col, "` are missing; Neyman ",
         "allocation needs it for every row it allocates over.", call. = FALSE)
  }
  vapply(idx_by_stratum, function(i) {
    if (length(i) < 2L) 0 else stats::sd(v[i])
  }, numeric(1), USE.NAMES = FALSE)
}

#' Strata, their sizes and their allocation
#'
#' The one place a stratified design's strata and per-stratum sample sizes are
#' worked out, so the draw, the inclusion probabilities and the joint
#' probabilities cannot disagree. Everything is by position: `code` is each
#' row's stratum as a position in `sizes`, `n_alloc` and `idx`, and `NA` for a
#' row with a missing key. `draw = TRUE` adds the checks only a draw needs.
#'
#' @noRd
stratum_allocation <- function(design, data, draw = FALSE) {
  validate_data(data, required_columns = design$strata)
  check_key_columns(data, design$strata, "strata")
  key <- group_key(data, design$strata)
  bad <- check_na_policy(is.na(key), design$na_rm, "a missing stratum key")
  if (draw) {
    if (all(bad)) stop("Every row has a missing stratum key.", call. = FALSE)
    check_draw_size(design$n, sum(!bad), design$replace)
  }
  code <- as.integer(key)
  idx <- split(seq_len(nrow(data))[!bad], key[!bad])
  sizes <- lengths(idx, use.names = FALSE)
  n_alloc <- allocate(design$n, sizes, design$allocation,
                      design$min_per_stratum, cap = !design$replace,
                      spread = stratum_spread(design, data, idx))
  list(code = code, labels = levels(key), idx = unname(idx), sizes = sizes,
       n_alloc = unname(n_alloc))
}

#' @export
draw_design.drawn_design_stratified <- function(design, data) {
  al <- stratum_allocation(design, data, draw = TRUE)
  sizes <- al$sizes
  n_alloc <- al$n_alloc
  warn_empty_groups(n_alloc, sizes, al$labels)

  if (!design$replace) {
    over <- which(n_alloc > sizes)
    if (length(over)) {
      h <- over[1]
      stop("Stratum `", group_label(al$labels[h]), "` has ", sizes[h],
           " row(s) but ", n_alloc[h], " were requested. ",
           "Use replace = TRUE or lower `n`.", call. = FALSE)
    }
    if (sum(n_alloc) != design$n) {
      stop("`n` (", design$n, ") cannot be allocated across these ",
           length(sizes), " strata with min_per_stratum = ",
           design$min_per_stratum, "; the allocation totals ", sum(n_alloc),
           " row(s).", call. = FALSE)
    }
  }

  reindex(data, take_within(al$idx, n_alloc, design$replace), sort = TRUE)
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
