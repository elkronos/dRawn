#' Stratified sampling
#'
#' Draws a sample from each stratum. Allocation uses the largest-remainder
#' method, so the returned row count matches `n` exactly rather than drifting
#' with per-stratum rounding.
#'
#' @param strata One or more column names defining the strata. Several columns
#'   are cross-classified.
#' @param n Total rows to draw across all strata, under either allocation.
#' @param allocation `"proportional"` splits `n` in proportion to stratum size;
#'   `"equal"` splits it evenly.
#' @param min_per_stratum Minimum rows from each stratum. The default of `0`
#'   leaves allocation unbiased; `1` guarantees coverage of rare strata at the
#'   cost of over-representing them.
#' @param replace Sample with replacement within each stratum?
#' @param na_rm Drop rows whose stratum key is `NA` instead of raising an error.
#'
#' @return A design object, for use with [draw()].
#'
#' @examples
#' df <- data.frame(id = 1:100, site = rep(letters[1:4], each = 25))
#' table(draw(df, design_stratified("site", n = 20), seed = 1)$site)
#'
#' # Rare strata are covered only if you ask
#' skewed <- data.frame(id = 1:1000, g = c(rep("common", 999), "rare"))
#' draw(skewed, design_stratified("g", n = 10, min_per_stratum = 1), seed = 1)$g
#'
#' @family designs
#' @seealso [draw()]
#' @export
design_stratified <- function(strata, n,
                              allocation = c("proportional", "equal"),
                              min_per_stratum = 0L, replace = FALSE,
                              na_rm = FALSE) {
  new_design("stratified", list(
    strata = check_columns(strata, "strata"),
    n = check_count(n, "n"),
    allocation = match.arg(allocation),
    min_per_stratum = check_count(min_per_stratum, "min_per_stratum"),
    replace = check_flag(replace, "replace"),
    na_rm = check_flag(na_rm, "na_rm")
  ))
}

#' @export
draw_design.drawn_design_stratified <- function(design, data) {
  validate_data(data, required_columns = design$strata)
  check_key_columns(data, design$strata, "strata")

  keys <- data[design$strata]
  bad <- Reduce(`|`, lapply(keys, is.na))
  if (any(bad)) {
    data <- drop_na_rows(data, bad, design$na_rm, "a missing stratum key")
    keys <- data[design$strata]
  }

  idx_by_stratum <- split(seq_len(nrow(data)),
                          interaction(keys, drop = TRUE, lex.order = TRUE))
  sizes <- lengths(idx_by_stratum)

  check_draw_size(design$n, nrow(data), design$replace)

  n_alloc <- allocate(design$n, sizes, design$allocation,
                      design$min_per_stratum, cap = !design$replace)

  if (!design$replace) {
    over <- n_alloc > sizes
    if (any(over)) {
      stop("Stratum `", names(sizes)[over][1], "` has ", sizes[over][1],
           " row(s) but ", n_alloc[over][1], " were requested. ",
           "Use replace = TRUE or lower `n`.", call. = FALSE)
    }
    if (sum(n_alloc) != design$n) {
      stop("`n` (", design$n, ") cannot be allocated across these ",
           length(sizes), " strata with min_per_stratum = ",
           design$min_per_stratum, "; the allocation totals ", sum(n_alloc),
           " row(s).", call. = FALSE)
    }
  }

  reindex(data, take_within(idx_by_stratum, n_alloc, design$replace), sort = TRUE)
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_stratified <- function(design, data) {
  validate_data(data, required_columns = design$strata)
  check_key_columns(data, design$strata, "strata")
  keys <- data[design$strata]
  bad <- Reduce(`|`, lapply(keys, is.na))
  group <- interaction(keys, drop = TRUE, lex.order = TRUE)

  idx_by_stratum <- split(seq_len(nrow(data))[!bad], group[!bad])
  sizes <- lengths(idx_by_stratum)
  n_alloc <- allocate(design$n, sizes, design$allocation,
                      design$min_per_stratum, cap = !design$replace)

  out <- numeric(nrow(data))
  for (h in names(idx_by_stratum)) {
    idx <- idx_by_stratum[[h]]
    out[idx] <- if (design$replace) {
      1 - (1 - 1 / length(idx))^n_alloc[[h]]
    } else {
      n_alloc[[h]] / length(idx)
    }
  }
  out
}
