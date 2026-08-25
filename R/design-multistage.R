#' Multi-stage sampling
#'
#' Selects clusters at the first stage, then draws rows within them at the
#' second. `n` is the total across the selected clusters under both allocations,
#' matching [design_stratified()].
#'
#' @param clusters Column naming each row's cluster.
#' @param n_clusters Number of clusters to select at stage one.
#' @param n Total rows to draw across the selected clusters. Make it a
#'   multiple of `n_clusters`, and no larger than `n_clusters` times the
#'   smallest cluster, if you intend to estimate: otherwise the per-cluster take
#'   depends on which clusters were selected and there is no closed-form
#'   inclusion probability. See [inclusion_prob()].
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
                      design$min_per_cluster, cap = !design$replace,
                      min_arg = "min_per_cluster")

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
  check_na_policy(!cl$present, design$na_rm, "a missing cluster key")
  labels <- cl$labels
  keep <- cl$present
  sizes <- table(labels[keep])

  # A closed form needs every selected cluster to contribute the same number of
  # rows, whichever clusters those turn out to be. Two things break that, and
  # both depend on the draw:
  #
  #   * `n` not divisible by `n_clusters` -- allocate() hands the remainder to
  #     the largest selected clusters, so a cluster's take depends on the company
  #     it keeps;
  #   * a cluster smaller than the per-cluster target -- it is capped at its own
  #     size and the shortfall is dealt to the others.
  #
  # Averaging over the possibilities is not a closed form, and the average is not
  # the answer: it was out by 13% on a four-cluster frame and could exceed 1.
  m <- design$n %/% design$n_clusters
  if (design$n %% design$n_clusters != 0L) {
    no_closed_form(
      paste0("`design_multistage()` with n = ", design$n, " over ",
             design$n_clusters, " clusters"),
      paste0("The remainder goes to the largest clusters selected, so a row's
",
             "probability depends on which other clusters were drawn. Choose an
",
             "`n` that divides by `n_clusters` (here, a multiple of ",
             design$n_clusters, "),
or pass simulate = TRUE.")
    )
  }
  small <- names(sizes)[sizes < m]
  if (length(small)) {
    no_closed_form(
      paste0("`design_multistage()` where ", length(small), " cluster(s) hold ",
             "fewer than ", m, " rows"),
      paste0("A cluster smaller than the per-cluster take is capped at its own
",
             "size and the shortfall is spread over whichever clusters came with
",
             "it, so there is no single answer. Lower `n`, drop the small
",
             "clusters, or pass simulate = TRUE.")
    )
  }
  if (design$min_per_cluster > m) {
    no_closed_form(
      paste0("`design_multistage(min_per_cluster = ", design$min_per_cluster,
             ")` above the per-cluster take of ", m),
      "Lower `min_per_cluster`, or pass simulate = TRUE."
    )
  }

  stage1 <- design$n_clusters / cl$total
  out <- numeric(nrow(data))
  out[keep] <- stage1 * (m / as.numeric(sizes[as.character(labels[keep])]))
  out
}
