# Stage-one cluster selection, shared by design_cluster() and
# design_multistage().
#
# Both designs select clusters the same way. Every later step compares
# cluster labels, from the probabilities and the variance to the degrees of
# freedom and the hand-off to `survey`. Comparing them as text merges two
# numbers that print alike, such as 0.3 and 0.1 + 0.2, and counting them with
# table() turns a factor's unused levels into clusters of no rows. Each row's
# cluster is therefore an integer code, and the codes are built in one place.

#' Integer codes for cluster labels
#'
#' Codes number the distinct labels in order of first appearance. Labels are
#' compared on their stored values: doubles exactly, dates and times on their
#' underlying numbers, factors by level, so a level no row uses has no code. A
#' missing label gets `NA`.
#'
#' @noRd
cluster_codes <- function(labels) {
  v <- cluster_values(labels)
  match(v, unique(v[!is.na(v)]))
}

#' The values cluster labels are compared on
#'
#' A factor's level numbers, the underlying numbers of a date or time, and
#' any other vector as it is.
#'
#' @noRd
cluster_values <- function(labels) {
  if (is.factor(labels)) {
    as.integer(labels)
  } else if (is.object(labels) && is.atomic(labels)) {
    unclass(labels)
  } else {
    labels
  }
}

#' Validate the cluster column and select clusters at random
#'
#' Clusters are numbered in order of first appearance and `n_clusters` of
#' those numbers are sampled, so the clusters chosen for a seed are the ones
#' `sample(unique(labels), n_clusters)` would give. Sampling the numbers also
#' keeps `sample()` from reading a single numeric cluster id as a range. The
#' rows of the chosen clusters are found in one pass over the frame.
#'
#' A multistage design with `size` selects the numbers by systematic PPS
#' over the clusters' stage-one probabilities in `cl$pi1`, with clusters at
#' probability 1 always among them.
#'
#' @param design A cluster or multistage design.
#' @param data The population.
#' @param cl `cluster_setup(design, data)`, when the caller has it already.
#' @return A list of `data` (with NA labels handled), the row indices of each
#'   selected cluster, their sizes, their labels, and their codes in
#'   `chosen`, all in the order the clusters were selected.
#' @noRd
select_clusters <- function(design, data, cl = cluster_setup(design, data)) {
  # Dropping rows with no label leaves every other row's code unchanged, since
  # codes number the labels in order of first appearance.
  if (!all(cl$present)) data <- reindex(data, which(cl$present))
  code <- cl$code[cl$present]
  labels <- cl$labels[cl$present]
  k <- cl$total

  chosen <- if (is.null(cl$pi1)) {
    sample.int(k, size = design$n_clusters)
  } else {
    systematic_pps(cl$pi1)
  }
  pos <- match(code, chosen)
  hit <- which(!is.na(pos))
  idx <- split(hit, structure(pos[hit],
                              levels = as.character(seq_along(chosen)),
                              class = "factor"))
  names(idx) <- NULL

  list(data = data, idx_by_cluster = idx, sizes = lengths(idx),
       labels = labels[match(chosen, code)], chosen = as.integer(chosen))
}

#' Each row's cluster, and how many distinct clusters the population holds
#'
#' Shared by the two designs' probabilities, variances and hand-off to
#' `survey`. `code` is each row's cluster from `cluster_codes()`, `NA` where
#' the label is missing, and `total` the number of clusters, the `A` in the
#' stage-one selection probability `n_clusters / A`. `values` holds the
#' value each code stands for, as `cluster_values()` gives it, so rows of
#' another frame can be coded alike.
#'
#' A multistage design with `size` cannot place a row whose size is missing,
#' so such a row has no code, like a row with no label, and the clusters are
#' numbered over the rows left.
#'
#' @noRd
count_clusters <- function(design, data) {
  validate_data(data, required_columns = c(design$clusters, design$size))
  check_key_columns(data, design$clusters, "clusters")
  labels <- data[[design$clusters]]
  v <- cluster_values(labels)
  if (!is.null(design$size)) {
    check_key_columns(data, design$size, "size")
    v[is.na(data[[design$size]])] <- NA
  }
  values <- unique(v[!is.na(v)])
  code <- match(v, values)
  present <- !is.na(code)
  list(labels = labels, code = code, present = present,
       total = if (any(present)) max(code[present]) else 0L,
       values = values)
}

#' Refuse a cluster selection the frame cannot supply
#'
#' The checks the draw and the probabilities of both cluster designs share:
#' the `na_rm` policy, a frame with no labelled row, and more clusters asked
#' for than it holds. Without them [inclusion_prob()] gave `n_clusters / A`
#' above 1 where [draw()] refused. Returns `count_clusters()`.
#'
#' @noRd
cluster_setup <- function(design, data) {
  cl <- count_clusters(design, data)
  unlabelled <- is.na(cl$labels)
  check_na_policy(unlabelled, design$na_rm, "a missing cluster label")
  if (!is.null(design$size)) {
    check_na_policy(!cl$present & !unlabelled, design$na_rm,
                    paste0("a missing size in `", design$size, "`"))
  }
  if (!any(cl$present)) {
    stop("Every row has a missing cluster label",
         if (!is.null(design$size)) " or size", ".", call. = FALSE)
  }
  if (design$n_clusters > cl$total) {
    stop("`n_clusters` (", design$n_clusters, ") exceeds the ", cl$total,
         " cluster(s) available.", call. = FALSE)
  }
  cl
}

# ---- weights conditional on the clusters selected ----------------------------

#' Does draw(weights = TRUE) weight a design conditional on its clusters?
#'
#' Balanced cluster sampling takes from each selected cluster as many rows as
#' the smallest of them holds, and a multistage design without `size` splits
#' `n` over the clusters it selected. When that split is not the same for
#' every selection, under proportional allocation, an `n` that `n_clusters`
#' does not divide, or a cluster smaller than `n / n_clusters`, a row's
#' probability depends on which other clusters were selected. Such a design
#' is weighted conditional on the clusters selected; see
#' `conditional_inclusion()`. A design drawn with replacement is not, since
#' draw() refuses it weights.
#'
#' @noRd
conditional_design <- function(design, data) {
  if (inherits(design, "drawn_design_cluster")) {
    return(isTRUE(design$balanced))
  }
  if (!inherits(design, "drawn_design_multistage") ||
      !is.null(design$size) || design$replace) {
    return(FALSE)
  }
  cl <- count_clusters(design, data)
  !multistage_fixed(design, tabulate(cl$code, cl$total))
}

#' Refuse a conditionally weighted draw before anything is sampled
#'
#' The checks of the design's own setup, and a floor of one row per selected
#' cluster. Proportional allocation with `min_per_cluster = 0` can give a
#' selected cluster no rows, and weights conditional on the selection would
#' then count that cluster's total as 0 while scaling the others by
#' `A / a` as if it had been measured. Returns the setup.
#'
#' @noRd
conditional_setup <- function(design, data) {
  if (inherits(design, "drawn_design_multistage") &&
      design$allocation == "proportional" && design$min_per_cluster < 1L) {
    stop("`min_per_cluster` is 0, so proportional allocation can give a ",
         "selected cluster no rows, and `weights = TRUE` weights each row ",
         "conditional on the clusters selected, which would count that ",
         "cluster's total as 0. Set min_per_cluster = 1 or more.",
         call. = FALSE)
  }
  if (inherits(design, "drawn_design_multistage")) {
    multistage_setup(design, data)
  } else {
    cluster_setup(design, data)
  }
}

#' The rows each selected cluster gives
#'
#' `sizes` holds the selected clusters' row counts, in the order they were
#' selected. A balanced cluster design takes the smallest of them from every
#' cluster, a whole-cluster design takes them all, and a multistage design
#' splits `n` with `multistage_takes()`.
#'
#' @noRd
cluster_takes <- function(design, sizes) {
  if (inherits(design, "drawn_design_multistage")) {
    return(multistage_takes(design, sizes))
  }
  if (design$balanced) rep(min(sizes), length(sizes)) else sizes
}

#' Inclusion probabilities conditional on the clusters selected
#'
#' For a design that `conditional_design()` accepts, `selection` holds the
#' codes of the clusters a draw selected, in the order it selected them,
#' which is the order the allocation breaks ties in. A row of selected
#' cluster `i`, of `N_i` rows, gets `(a / A) * m_i / N_i`: the chance that
#' its cluster is selected, times the chance that it is drawn once these
#' clusters are, with `m_i` the rows the design takes from cluster `i` given
#' this selection. Weighting by its reciprocal turns each selected cluster's
#' rows into `A / a` times an unbiased estimate of the cluster's total, so
#' the estimate of the total is unbiased given the selection, and therefore
#' over every selection.
#'
#' A row of a cluster that was not selected gets `NA`. The selection does not
#' fix its probability, and the design has no closed form for it. A row the
#' design drops, for a missing cluster label, gets 0.
#'
#' @noRd
conditional_inclusion <- function(design, data, selection) {
  cl <- conditional_setup(design, data)
  sizes <- tabulate(cl$code, cl$total)
  take <- cluster_takes(design, sizes[selection])
  if (any(take < 1L)) {
    stop("The allocation gives ", sum(take < 1L), " of the clusters ",
         "selected no rows, so weights conditional on the selection would ",
         "count their totals as 0. Set min_per_cluster = 1 or more.",
         call. = FALSE)
  }
  p <- rep(NA_real_, cl$total)
  p[selection] <- design$n_clusters / cl$total * take / sizes[selection]
  out <- numeric(nrow(data))
  out[cl$present] <- p[cl$code[cl$present]]
  out
}
