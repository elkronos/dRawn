#' Multi-stage sampling
#'
#' Selects clusters at the first stage, then draws rows within them at the
#' second. `n` is the total across the selected clusters under both allocations,
#' matching [design_stratified()].
#'
#' @param clusters Column naming each row's cluster.
#' @param n_clusters Number of clusters to select at stage one.
#' @param n Total rows to draw across the selected clusters.
#' @param allocation `"equal"` splits `n` evenly across the selected clusters;
#'   `"proportional"` splits it in proportion to their size.
#' @param min_per_cluster Minimum rows from each selected cluster. Defaults to
#'   `0`, which leaves allocation unbiased.
#' @param replace Sample with replacement within each cluster?
#' @param na_rm Drop rows whose cluster label is `NA` instead of raising an
#'   error.
#'
#' @return A design object, for use with [draw()].
#'
#' @examples
#' df <- data.frame(id = 1:100, site = rep(paste0("s", 1:10), each = 10))
#' nrow(draw(df, design_multistage("site", n_clusters = 4, n = 12), seed = 1))
#'
#' @family designs
#' @seealso [draw()]
#' @export
design_multistage <- function(clusters, n_clusters, n,
                              allocation = c("equal", "proportional"),
                              min_per_cluster = 0L, replace = FALSE,
                              na_rm = FALSE) {
  new_design("multistage", list(
    clusters = check_columns(clusters, "clusters", 1L),
    n_clusters = check_count(n_clusters, "n_clusters", allow_zero = FALSE),
    n = check_count(n, "n"),
    allocation = match.arg(allocation),
    min_per_cluster = check_count(min_per_cluster, "min_per_cluster"),
    replace = check_flag(replace, "replace"),
    na_rm = check_flag(na_rm, "na_rm")
  ))
}

#' @export
draw_design.drawn_design_multistage <- function(design, data) {
  sel <- select_clusters(design, data)
  data <- sel$data
  sizes <- sel$sizes

  # Not check_draw_size(): the bound here is the rows across the SELECTED
  # clusters, not nrow(data), and the generic message would name the wrong
  # number.
  if (!design$replace && design$n > sum(sizes)) {
    stop("`n` (", design$n, ") exceeds the ", sum(sizes),
         " row(s) available across the ", length(sizes),
         " selected cluster(s). Use replace = TRUE, raise `n_clusters`, or ",
         "lower `n`.", call. = FALSE)
  }

  # The denominator is the rows in the SELECTED clusters, not the population.
  n_alloc <- allocate(design$n, sizes, design$allocation,
                      design$min_per_cluster, cap = !design$replace)

  if (!design$replace && sum(n_alloc) != design$n) {
    stop("`n` (", design$n, ") cannot be allocated across these ",
         length(sizes), " cluster(s) with min_per_cluster = ",
         design$min_per_cluster, "; the allocation totals ", sum(n_alloc),
         " row(s).", call. = FALSE)
  }

  reindex(data, take_within(sel$idx_by_cluster, n_alloc, design$replace),
          sort = TRUE)
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_multistage <- function(design, data) {
  if (design$allocation == "proportional") {
    no_closed_form(
      "`design_multistage(allocation = \"proportional\")`",
      paste0("The stage-two allocation depends on which clusters were ",
             "selected.\nUse allocation = \"equal\" for a closed form.")
    )
  }
  cl <- count_clusters(design, data)
  labels <- cl$labels
  keep <- cl$present
  sizes <- table(labels[keep])

  # Equal allocation: every selected cluster contributes the same n_h,
  # so the second stage does not depend on which clusters were picked.
  per_cluster <- allocate(design$n, rep(1L, design$n_clusters), "equal",
                          design$min_per_cluster, cap = FALSE)
  # allocate() spreads a remainder across clusters; with equal-sized targets the
  # per-cluster take differs by at most one, so average over the positions a
  # cluster could occupy.
  mean_take <- mean(per_cluster)

  stage1 <- design$n_clusters / cl$total
  out <- numeric(nrow(data))
  out[keep] <- stage1 * (mean_take / as.numeric(sizes[as.character(labels[keep])]))
  out
}
