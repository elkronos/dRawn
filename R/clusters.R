# Stage-one cluster selection, shared by design_cluster() and
# design_multistage().
#
# Both designs select clusters the same way, and the two subtle traps here --
# sample() reinterpreting a length-one numeric cluster id as a range, and `==`
# against an NA label fabricating all-NA rows -- are easy to reintroduce. One
# copy means one place to get them right.

#' Validate the cluster column and select clusters at random
#'
#' @param design A cluster or multistage design.
#' @param data The population.
#' @return A list of `data` (with NA labels handled), the row indices of each
#'   selected cluster, and their sizes.
#' @noRd
select_clusters <- function(design, data) {
  validate_data(data, required_columns = design$clusters)
  check_key_columns(data, design$clusters, "clusters")

  col <- design$clusters
  data <- drop_na_rows(data, is.na(data[[col]]), design$na_rm,
                       "a missing cluster label")
  labels <- data[[col]]
  levels_present <- unique(labels)

  if (design$n_clusters > length(levels_present)) {
    stop("`n_clusters` (", design$n_clusters, ") exceeds the ",
         length(levels_present), " cluster(s) available.", call. = FALSE)
  }

  chosen <- sample_values(levels_present, design$n_clusters)

  # `%in%`, never `==`: NA == "a" is NA, and indexing with NA manufactures
  # all-NA rows that were never in the data.
  idx <- lapply(chosen, function(cl) which(labels %in% cl))
  names(idx) <- as.character(chosen)

  list(data = data, labels = labels, chosen = chosen,
       idx_by_cluster = idx, sizes = lengths(idx))
}

#' How many distinct clusters the population holds
#'
#' Shared by the two designs' inclusion probabilities, both of which need the
#' stage-one selection probability `n_clusters / n_clusters_total`.
#'
#' @noRd
count_clusters <- function(design, data) {
  validate_data(data, required_columns = design$clusters)
  check_key_columns(data, design$clusters, "clusters")
  labels <- data[[design$clusters]]
  list(labels = labels,
       present = !is.na(labels),
       total = length(unique(labels[!is.na(labels)])))
}
