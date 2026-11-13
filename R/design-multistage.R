#' Multi-stage sampling
#'
#' Selects clusters at the first stage, then draws rows within them at the
#' second. `n` is the total across the selected clusters under both allocations,
#' matching [design_stratified()].
#'
#' @section Clusters in proportion to size:
#' With `size`, the clusters are selected with probability proportional to a
#' size measure, such as a count of households from a census, by
#' systematic PPS as in [design_weighted()]: cluster `i` gets
#' `pi_i = n_clusters * M_i / sum(M)`, and a cluster large enough that this
#' passes 1 is taken with certainty, the rest of `n_clusters` being spread
#' over the others. Each selected cluster then gives the same `m = n /
#' n_clusters` rows by simple random sampling, so a row of cluster `i`, of
#' `N_i` rows, has inclusion probability `pi_i * m / N_i`. When the measure
#' is the cluster's row count and no cluster is certain, every row gets
#' `n / N`: the design is self-weighting, with a fixed number of interviews
#' per cluster.
#'
#' The measure is a property of the cluster, so every row of a cluster must
#' carry the same value; a cluster with two values is refused. `n` must be a
#' multiple of `n_clusters`, and `allocation` must be `"equal"`.
#'
#' [ht_total()] and [ht_mean()] estimate the variance as Deville's
#' approximation over the estimated totals of the selected clusters below
#' certainty, each scaled by its own `1 - pi_i`, plus the variance within the
#' selected clusters, each divided by its `pi_i`. That is the two-stage
#' estimator of [design_multistage()] without `size` when every `pi_i` is
#' equal. A cluster taken with certainty is a stratum, and adds only the
#' variance of the simple random sample within it. The degrees of freedom
#' count the clusters below certainty less one, and in each certain cluster
#' the sampled rows less one. The variance is declined when one row is taken
#' per cluster, whose variation within clusters cannot be measured, or when
#' one cluster below certainty is selected. In simulations of 2,000 draws
#' each from frames of 78 or 80 clusters of 10 to 600 rows, taking 2 or 5
#' rows from each of 12 or 20 clusters, with and without clusters taken
#' with certainty, the variance averaged 0.94 to 1.02 times the true
#' variance of the total, and 95% intervals covered 93% to 95% of the time.
#'
#' [as_svydesign()] hands the design to `survey` in two stages with
#' Brewer's approximation at the first. Joint probabilities have no closed
#' form; [joint_prob()] estimates them with `simulate = TRUE`.
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
#' @param na_rm Drop rows whose cluster label, or `size`, is `NA` instead of
#'   raising an error.
#' @param size Optional column holding each cluster's size measure, the same
#'   on every row of a cluster: positive, finite and numeric. Clusters are
#'   then selected with probability proportional to it. See "Clusters in
#'   proportion to size".
#'
#' @return A design object, for use with [draw()].
#'
#' @examples
#' df <- data.frame(id = 1:100, site = rep(paste0("s", 1:10), each = 10))
#' nrow(draw(df, design_multistage("site", n_clusters = 4, n = 12), seed = 1))
#'
#' # Villages selected in proportion to their households, then 5 households
#' # from each: every household has the same chance, 20 / 260
#' hh <- data.frame(village = rep(paste0("v", 1:8),
#'                                times = c(10, 20, 30, 40, 20, 50, 60, 30)))
#' hh$households <- ave(seq_len(nrow(hh)), hh$village, FUN = length)
#' pps <- design_multistage("village", n_clusters = 4, n = 20,
#'                          size = "households")
#' unique(round(inclusion_prob(hh, pps), 4))
#'
#' @family designs
#' @seealso [draw()]
#' @export
design_multistage <- function(clusters, n_clusters, n,
                              allocation = c("equal", "proportional"),
                              min_per_cluster = 0L, replace = FALSE,
                              na_rm = FALSE, size = NULL) {
  clusters <- check_columns(clusters, "clusters", 1L)
  n_clusters <- check_count(n_clusters, "n_clusters", allow_zero = FALSE)
  n <- check_count(n, "n")
  min_per_cluster <- check_count(min_per_cluster, "min_per_cluster")
  allocation <- check_choice(allocation, c("equal", "proportional"),
                             "allocation")
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
  if (!is.null(size)) {
    size <- check_columns(size, "size", 1L)
    if (allocation != "equal") {
      stop("`allocation` is \"", allocation, "\", but with `size` every ",
           "selected cluster gives the same number of rows. Leave ",
           "`allocation` out, or set it to \"equal\".", call. = FALSE)
    }
    if (n %% n_clusters != 0L) {
      lower <- n %/% n_clusters * n_clusters
      stop("`n` (", n, ") is not a multiple of `n_clusters` (", n_clusters,
           "), and with `size` every selected cluster gives the same number ",
           "of rows. Use n = ", if (lower > 0L) paste0(lower, " or n = "),
           lower + n_clusters, ".", call. = FALSE)
    }
  }
  new_design("multistage", list(
    clusters = clusters,
    n_clusters = n_clusters,
    n = n,
    allocation = allocation,
    min_per_cluster = min_per_cluster,
    replace = check_flag(replace, "replace"),
    na_rm = check_flag(na_rm, "na_rm"),
    size = size
  ))
}

#' Refuse a multi-stage design the frame cannot supply
#'
#' The cluster checks, and an `n` above the rows held by the `n_clusters`
#' largest clusters, which no selection can supply. Shared by the draw and the
#' probabilities, so both refuse with the same message. An `n` that only some
#' selections can supply is the draw's to refuse, against the clusters it
#' selected; the probabilities have no closed form there. Returns
#' `count_clusters()`, with `pi1`, each cluster's stage-one probability by
#' code, for a design with `size`.
#'
#' @noRd
multistage_setup <- function(design, data) {
  cl <- cluster_setup(design, data)
  if (!is.null(design$size)) cl$pi1 <- cluster_pi(design, data, cl)
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

#' Each cluster's stage-one probability, in proportion to its size measure
#'
#' `pps_pi()` over one measure per cluster, in code order. The measure is
#' read from the rows of each cluster, which must agree. `cl` is
#' `cluster_setup()`.
#'
#' @noRd
cluster_pi <- function(design, data, cl) {
  col <- design$size
  v <- data[[col]]
  check_weights(v, col)
  code <- cl$code[cl$present]
  v <- v[cl$present]
  first <- match(seq_len(cl$total), code)
  measure <- v[first]
  off <- which(v != measure[code])
  if (length(off)) {
    i <- off[1]
    label <- cl$labels[cl$present][i]
    stop("`", col, "` must hold one size for each cluster, but cluster `",
         group_label(as.character(label)), "` has ", format(measure[code[i]]),
         " and ", format(v[i]), ". Give every row of a cluster the ",
         "cluster's size measure.", call. = FALSE)
  }
  pps_pi(measure, design$n_clusters)
}

#' The stage-one selection of a multistage design with `size`, from its frame
#'
#' Each frame row's cluster `code`, `NA` for a row the design drops; the
#' `values` those codes stand for; each cluster's stage-one probability `pi1`
#' and row count `sizes`, by code; and `m`, the rows taken from each selected
#' cluster. The variance, the degrees of freedom, the jackknife and
#' `as_svydesign()` read it.
#'
#' @noRd
cluster_stage <- function(design, data) {
  cl <- multistage_setup(design, data)
  list(code = cl$code, values = cl$values, pi1 = cl$pi1,
       sizes = tabulate(cl$code, cl$total),
       m = design$n %/% design$n_clusters)
}

#' The sampling units of a multistage sample with `size`
#'
#' A cluster selected below certainty is a primary sampling unit, and all of
#' them share one stratum. A cluster taken with certainty is a stratum of its
#' own whose units are its sampled rows, since only they were drawn at
#' random; one taken whole as well is a single unit, with no variance. For
#' the sampled rows `units` of the frame `pop`, returns each row's cluster
#' `code` in the frame's numbering, its unit `psu` and its `stratum`: 1 below
#' certainty, and `1 + code` for a cluster taken with certainty. `pop_stratum`
#' gives the stratum of every frame row, `NA` for a row the design drops.
#'
#' @noRd
cluster_units <- function(design, units, pop, ctx = NULL) {
  st <- from_context(ctx, "stage_one", cluster_stage(design, pop))
  code <- match(cluster_values(units[[design$clusters]]), st$values)
  A <- length(st$pi1)
  certain <- is_certain(st$pi1)
  whole <- certain & st$sizes <= st$m
  row_unit <- certain[code] & !whole[code]
  psu <- code
  psu[row_unit] <- A + which(row_unit)
  list(code = code, psu = psu,
       stratum = ifelse(certain[code], 1L + code, 1L),
       pop_stratum = ifelse(certain[st$code], 1L + st$code, 1L))
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

  stage1 <- if (is.null(cl$pi1)) {
    design$n_clusters / cl$total
  } else {
    cl$pi1[code[keep]]
  }
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
