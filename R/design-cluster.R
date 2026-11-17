#' Cluster sampling
#'
#' Selects whole clusters at random and returns all of their rows. The row count
#' follows from which clusters were picked, so there is no `n`; use
#' [design_multistage()] when you need to control it.
#'
#' @param clusters Column naming each row's cluster.
#' @param n_clusters Number of clusters to select.
#' @param balanced Take an equal number of rows from each selected cluster,
#'   equal to the smallest selected cluster's size. See "Balanced clusters".
#' @param na_rm Drop rows whose cluster label is `NA` instead of raising an
#'   error. When `FALSE`, missing labels are never treated as a cluster of their
#'   own.
#'
#' @section Balanced clusters:
#' With `balanced = TRUE`, each of the `a` clusters selected out of `A`
#' gives `m` rows by simple random sampling, `m` being the size of the
#' smallest of them. A row of a cluster of `N_i` rows has inclusion
#' probability `(a / A) * E(m) / N_i`, where `E(m)` is the average take
#' over the selections that include its cluster. [inclusion_prob()]
#' computes it exactly: given that cluster `i` is selected, the take is at
#' least `s` rows when the other `a - 1` clusters selected all hold `s` rows
#' or more, which happens with probability
#' `choose(c_s - 1, a - 1) / choose(A - 1, a - 1)`, `c_s` being the number
#' of clusters of `s` rows or more, and `E(m)` adds these up for `s` from 1
#' to `N_i`. On clusters of 2, 4, 6 and 8 rows with two selected, that gives
#' 0.5, 0.417, 0.333 and 0.25.
#'
#' `draw(weights = TRUE)` weights a balanced sample conditional on the
#' clusters it selected instead: `.prob` is `(a / A) * m / N_i` for the `m`
#' taken in this draw. The estimate of a total is then `A / a` times the sum
#' of the selected clusters' estimated totals, which is unbiased whichever
#' clusters were selected, and so over all of them. [ht_total()] and
#' [ht_mean()] use the two-stage variance of [design_multistage()], with the
#' variance between the estimated cluster totals and within each cluster,
#' and [as_svydesign()] hands the sample to `survey` in two stages. The
#' smallest selected cluster is taken whole and adds no variance within it.
#' When it holds one row, every cluster gives one, the variation within
#' clusters cannot be measured, and the variance falls back to the
#' delete-a-cluster jackknife. Joint probabilities have no closed form;
#' [joint_prob()] estimates them with `simulate = TRUE`.
#'
#' @return A design object, for use with [draw()].
#'
#' @examples
#' df <- data.frame(id = 1:100, site = rep(paste0("s", 1:10), each = 10))
#' unique(draw(df, design_cluster("site", n_clusters = 3), seed = 1)$site)
#'
#' # Balanced: the take is the smallest selected cluster's size, so a row of
#' # a large cluster is less likely to be drawn
#' uneven <- data.frame(site = rep(c("a", "b", "c", "d"), c(2, 4, 6, 8)))
#' bal <- design_cluster("site", n_clusters = 2, balanced = TRUE)
#' round(tapply(inclusion_prob(uneven, bal), uneven$site, unique), 3)
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
  keep <- if (design$balanced) {
    take_within(sel$idx_by_cluster, cluster_takes(design, sel$sizes), FALSE)
  } else {
    unlist(sel$idx_by_cluster, use.names = FALSE)
  }
  out <- reindex(data, keep, sort = TRUE)
  # draw(weights = TRUE) weights a balanced sample conditional on these
  # clusters, and removes the attribute otherwise.
  attr(out, "drawn_selection") <- sel$chosen
  out
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_cluster <- function(design, data) {
  cl <- cluster_setup(design, data)
  a <- design$n_clusters
  out <- numeric(nrow(data))
  if (design$balanced) {
    sizes <- tabulate(cl$code, cl$total)
    p <- a / cl$total * balanced_take(sizes, a) / sizes
    out[cl$present] <- p[cl$code[cl$present]]
    return(out)
  }
  out[cl$present] <- a / cl$total
  out
}

#' Each cluster's average take under balanced cluster sampling
#'
#' `a` of the `A` clusters, of `sizes` rows, are selected by simple random
#' sampling, and each gives as many rows as the smallest of them holds. For
#' cluster `i` this is that take averaged over the selections that include
#' `i`. Given `i`, the other `a - 1` are a simple random sample of the other
#' `A - 1` clusters, and the take is at least `s`, for `s` up to `N_i`, when
#' all of them hold `s` rows or more: probability
#' `choose(c_s - 1, a - 1) / choose(A - 1, a - 1)`, with `c_s` the number of
#' clusters, `i` among them, of `s` rows or more. The average take is the sum
#' of these over `s = 1, ..., N_i`. `c_s` changes only at the cluster sizes,
#' so the sum runs over the distinct sizes, in increasing order, each term
#' covering the values of `s` up to it from the size before. Computed on the
#' log scale, since the counts of selections overflow for a few hundred
#' clusters.
#'
#' @noRd
balanced_take <- function(sizes, a) {
  A <- length(sizes)
  d <- sort(unique(sizes))
  # Clusters of d[k] rows or more: those not of d[k] - 1 rows or fewer.
  c_s <- A - findInterval(d - 1, sort(sizes))
  p <- exp(lchoose(c_s - 1, a - 1) - lchoose(A - 1, a - 1))
  cumsum(diff(c(0, d)) * p)[match(sizes, d)]
}
