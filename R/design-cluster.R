#' Cluster sampling
#'
#' Selects whole clusters at random and returns all of their rows. The row count
#' follows from which clusters were picked, so there is no `n`; use
#' [design_multistage()] when you need to control it.
#'
#' @param clusters Column naming each row's cluster.
#' @param n_clusters Number of clusters to select.
#' @param balanced Take an equal number of rows from each selected cluster,
#'   equal to the smallest selected cluster's size.
#' @param na_rm Drop rows whose cluster label is `NA` instead of raising an
#'   error. When `FALSE`, missing labels are never treated as a cluster of their
#'   own.
#'
#' @return A design object, for use with [draw()].
#'
#' @examples
#' df <- data.frame(id = 1:100, site = rep(paste0("s", 1:10), each = 10))
#' unique(draw(df, design_cluster("site", n_clusters = 3), seed = 1)$site)
#'
#' @family designs
#' @seealso [draw()]
#' @export
design_cluster <- function(clusters, n_clusters, balanced = FALSE,
                           na_rm = FALSE) {
  new_design("cluster", list(
    clusters = check_columns(clusters, "clusters", 1L),
    n_clusters = check_count(n_clusters, "n_clusters", allow_zero = FALSE),
    balanced = check_flag(balanced, "balanced"),
    na_rm = check_flag(na_rm, "na_rm")
  ))
}

#' @export
draw_design.drawn_design_cluster <- function(design, data) {
  sel <- select_clusters(design, data)
  data <- sel$data
  keep <- unlist(sel$idx_by_cluster, use.names = FALSE)

  if (design$balanced) {
    n_min <- min(sel$sizes)
    keep <- take_within(sel$idx_by_cluster, rep(n_min, length(sel$sizes)), FALSE)
  }

  reindex(data, keep, sort = TRUE)
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_cluster <- function(design, data) {
  if (design$balanced) {
    no_closed_form(
      "`design_cluster(balanced = TRUE)`",
      paste0("The per-cluster take is the smallest SELECTED cluster's size, ",
             "which is random.\nUse balanced = FALSE for a closed form.")
    )
  }
  cl <- count_clusters(design, data)
  out <- rep(design$n_clusters / cl$total, nrow(data))
  out[!cl$present] <- 0
  out
}
