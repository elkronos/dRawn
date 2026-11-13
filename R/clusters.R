# Stage-one cluster selection, shared by design_cluster() and
# design_multistage().
#
# Both designs select clusters the same way, and every later step -- the
# probabilities, the variance, the degrees of freedom, the hand-off to
# `survey` -- compares cluster labels. Comparing them as text merges two
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
#'   selected cluster, their sizes, and their labels.
#' @noRd
select_clusters <- function(design, data, cl = cluster_setup(design, data)) {
  # Dropping rows with no label leaves every other row's code unchanged, since
  # codes number the labels in order of first appearance.
  if (!all(cl$present)) data <- data[cl$present, , drop = FALSE]
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
       labels = labels[match(chosen, code)])
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
