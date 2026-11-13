#' Multi-stage sampling
#'
#' Selects clusters at the first stage, then draws rows within them at the
#' second. `n` is the total across the selected clusters under both allocations,
#' matching [design_stratified()].
#'
#' @param clusters Column naming each row's cluster.
#' @param n_clusters Number of clusters to select at stage one.
#' @param n Total rows to draw across the selected clusters, at least
#'   `n_clusters`. Make it a multiple of `n_clusters`, and no larger than
#'   `n_clusters` times the smallest cluster, if you intend to estimate:
#'   otherwise the per-cluster take depends on which clusters were selected
#'   and there is no closed-form inclusion probability. See [inclusion_prob()].
#' @param allocation `"equal"` splits `n` evenly across the selected clusters;
#'   `"proportional"` splits it in proportion to their size. Proportional
#'   allocation can give a small selected cluster no rows; [draw()] then warns
#'   with class `drawn_empty_cluster` and names it.
#' @param min_per_cluster Minimum rows from each selected cluster. Defaults to
#'   `0`, and `n_clusters * min_per_cluster` must not exceed `n`. A cluster
#'   smaller than the minimum is taken whole. Take at least two rows per
#'   cluster if you intend to estimate a variance: with one, the variation
#'   within clusters cannot be measured and [ht_total()] falls back to the
#'   ultimate-cluster jackknife.
#' @param replace Sample with replacement within each cluster? A selected
#'   cluster then gives its full share of `n` however small it is, and a row's
#'   inclusion probability is the chance of being drawn at least once,
#'   `(n_clusters / N_clusters) * (1 - (1 - 1 / N_h)^m)` for a take of `m`
#'   from a cluster of `N_h` rows.
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
  clusters <- check_columns(clusters, "clusters", 1L)
  n_clusters <- check_count(n_clusters, "n_clusters", allow_zero = FALSE)
  n <- check_count(n, "n")
  min_per_cluster <- check_count(min_per_cluster, "min_per_cluster")
  if (n < n_clusters) {
    stop("`n` (", n, ") is smaller than `n_clusters` (", n_clusters,
         "), so some selected clusters would give no rows. Raise `n` to at ",
         "least ", n_clusters, ", or lower `n_clusters`.", call. = FALSE)
  }
  if (min_per_cluster * n_clusters > n) {
    stop("`min_per_cluster` = ", min_per_cluster, " across ", n_clusters,
         " clusters needs at least ", min_per_cluster * n_clusters,
         " rows, but `n` is ", n, ". Raise `n` or lower `min_per_cluster`.",
         call. = FALSE)
  }
  new_design("multistage", list(
    clusters = clusters,
    n_clusters = n_clusters,
    n = n,
    allocation = check_choice(allocation, c("equal", "proportional"),
                              "allocation"),
    min_per_cluster = min_per_cluster,
    replace = check_flag(replace, "replace"),
    na_rm = check_flag(na_rm, "na_rm")
  ))
}

#' Refuse a multi-stage design the frame cannot supply
#'
#' The cluster checks, and an `n` above the rows held by the `n_clusters`
#' largest clusters, which no selection can supply. Shared by the draw and the
#' probabilities, so both refuse with the same message. An `n` that only some
#' selections can supply is the draw's to refuse, against the clusters it
#' selected; the probabilities have no closed form there. Returns
#' `count_clusters()`.
#'
#' @noRd
multistage_setup <- function(design, data) {
  cl <- cluster_setup(design, data)
  if (!design$replace) {
    sizes <- tabulate(cl$code[cl$present], cl$total)
    most <- sum(sort(sizes, decreasing = TRUE)[seq_len(design$n_clusters)])
    if (design$n > most) {
      stop("`n` (", design$n, ") exceeds the ", most, " row(s) available ",
           "across the ", design$n_clusters, " selected cluster(s), even ",
           "with the largest selected. Use replace = TRUE, raise ",
           "`n_clusters`, or lower `n`.", call. = FALSE)
    }
  }
  cl
}

#' @export
draw_design.drawn_design_multistage <- function(design, data) {
  sel <- select_clusters(design, data, multistage_setup(design, data))
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

  warn_empty_clusters(n_alloc, sel$labels, design$n)

  reindex(data, take_within(sel$idx_by_cluster, n_alloc, design$replace),
          sort = TRUE)
}

#' Warn when the allocation leaves a selected cluster with no rows
#'
#' Proportional allocation can round a small cluster's share down to nothing.
#' The cluster was selected, but this draw takes none of its rows, so its rows
#' are missing from the sample. Unlike an empty stratum, such a cluster may get
#' rows in another draw, alongside other clusters, so the warning speaks of
#' this draw only. A cluster smaller than `min_per_cluster` is taken whole and
#' is not empty.
#'
#' @noRd
warn_empty_clusters <- function(n_alloc, labels, n) {
  empty <- which(n_alloc == 0L)
  if (!length(empty)) return(invisible(NULL))
  shown <- paste0("`", group_label(as.character(labels[utils::head(empty, 5L)])),
                  "`", collapse = ", ")
  more <- if (length(empty) > 5L) paste0(" and ", length(empty) - 5L, " more")
          else ""
  plural <- length(empty) > 1L
  warn_drawn(paste0(
    "The allocation gives ",
    if (plural) paste0(length(empty), " selected clusters (", shown, more, ")")
    else paste0("selected cluster ", shown),
    " none of the ", n, " rows, so this draw takes nothing from ",
    if (plural) "them" else "it", ". Set min_per_cluster = 1 to take at least ",
    "one row from every selected cluster."
  ), "drawn_empty_cluster")
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_multistage <- function(design, data) {
  cl <- multistage_setup(design, data)
  if (design$allocation == "proportional") {
    no_closed_form(
      "`design_multistage(allocation = \"proportional\")`",
      paste0("The stage-two allocation depends on which clusters were ",
             "selected.\nUse allocation = \"equal\" for a closed form.")
    )
  }
  keep <- cl$present
  code <- cl$code
  sizes <- tabulate(code, cl$total)

  # A closed form needs every selected cluster to contribute the same number of
  # rows, whichever clusters those turn out to be. Two things break that, and
  # both depend on the draw:
  #
  #   * `n` not divisible by `n_clusters` -- allocate() hands the remainder to
  #     the largest selected clusters, so a cluster's take depends on the company
  #     it keeps;
  #   * a cluster smaller than the per-cluster target -- it is capped at its own
  #     size and the shortfall is dealt to the others. Sampling with
  #     replacement has no cap, so this one applies only without it.
  #
  # Averaging over the possibilities is not a closed form, and the average is not
  # the answer: it was out by 13% on a four-cluster frame and could exceed 1.
  m <- design$n %/% design$n_clusters
  if (design$n %% design$n_clusters != 0L) {
    no_closed_form(
      paste0("`design_multistage()` with n = ", design$n, " over ",
             design$n_clusters, " clusters"),
      paste0("The remainder goes to the largest clusters selected, so a row's ",
             "probability\ndepends on which other clusters were drawn. Choose ",
             "an `n` that divides by\n`n_clusters` (here, a multiple of ",
             design$n_clusters, ").")
    )
  }
  small <- if (design$replace) integer(0) else which(sizes < m)
  if (length(small)) {
    no_closed_form(
      paste0("`design_multistage()` where ", length(small), " cluster(s) hold ",
             "fewer than ", m, " rows"),
      paste0("A cluster smaller than the per-cluster take is capped at its own ",
             "size and the\nshortfall is spread over whichever clusters came ",
             "with it, so there is no single\nanswer. Lower `n`, or drop the ",
             "small clusters.")
    )
  }
  if (design$min_per_cluster > m) {
    no_closed_form(
      paste0("`design_multistage(min_per_cluster = ", design$min_per_cluster,
             ")` above the per-cluster take of ", m),
      "Lower `min_per_cluster`."
    )
  }

  stage1 <- design$n_clusters / cl$total
  N_h <- as.numeric(sizes[code[keep]])
  out <- numeric(nrow(data))
  # With replacement, m draws from a cluster of N_h rows reach a given row at
  # least once with probability 1 - (1 - 1/N_h)^m. m / N_h is the expected
  # number of times it is drawn, which overstates the probability.
  out[keep] <- stage1 * if (design$replace) {
    1 - (1 - 1 / N_h)^m
  } else {
    m / N_h
  }
  out
}
