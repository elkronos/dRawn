#' Joint inclusion probabilities
#'
#' The probability that a *pair* of rows both land in the sample. First-order
#' probabilities from [inclusion_prob()] give you an unbiased total;
#' second-order probabilities are what let you put a standard error on it.
#'
#' @section Which designs have them:
#' Closed forms exist, and are used, for the designs whose selection is either
#' independent across groups or a simple random sample within them:
#'
#' \tabular{ll}{
#'   [design_simple()]      \tab `n(n-1) / (N(N-1))` for a pair, without replacement \cr
#'   [design_stratified()]  \tab within a stratum as above; across strata, independent \cr
#'   [design_cluster()]     \tab same cluster: `a/A`; different clusters: `a(a-1)/(A(A-1))` \cr
#'   [design_multistage()]  \tab the two stages multiplied, where the per-cluster take is constant \cr
#'   [design_certainty()]   \tab `1` between certainty rows; otherwise the other row's own `pi` \cr
#'   [design_reservoir()]   \tab simple random sampling over the first `max_items` rows \cr
#'   [design_temporal()]    \tab within an interval as above; across intervals, independent \cr
#'   [design_spatial()]     \tab simple random sampling inside the region \cr
#'   [design_weighted()]    \tab `"poisson"` only, where rows are independent: `pi_i * pi_j` \cr
#'   [design_systematic()]  \tab `1/interval` for rows sharing a residue class, otherwise **zero** \cr
#' }
#'
#' Systematic sampling is the awkward one. Most pairs can never co-occur, so
#' their joint probability is genuinely 0 and no design-unbiased variance
#' estimator exists. [ht_total()] says so rather than returning a number. Its
#' residue classes follow the order the design walks, so `order_by` changes
#' which pairs can co-occur.
#'
#' `design_weighted(method = "systematic")` has joint probabilities, but they
#' depend on the order units are visited and need a dedicated algorithm. Use
#' `sampling::UPsystematicpi2()` for those.
#'
#' @param data A data frame.
#' @param design A design object.
#' @param rows Optional row indices. Supply these — usually the rows you drew —
#'   to get the submatrix for them instead of the full `nrow(data)` square,
#'   which is what makes this usable on a large population.
#' @param simulate Estimate the probabilities by repeated draws rather than in
#'   closed form. Works for every probability design, including the ones with no
#'   closed form, at the cost of Monte Carlo error. Refused for
#'   [design_bootstrap()], where every row appears in some replicate and the
#'   count converges to 1 for all of them.
#' @param R Number of simulated draws when `simulate = TRUE`.
#' @param seed Optional seed for the simulation.
#'
#' @return A square matrix with one row and column per element of `rows`
#'   (or per row of `data`). The diagonal holds first-order probabilities.
#'
#' @examples
#' df <- data.frame(id = 1:10)
#' round(joint_prob(df, design_simple(n = 4)), 3)
#'
#' # Only for the rows you drew
#' joint_prob(df, design_simple(n = 4), rows = c(2, 5, 7))
#'
#' @seealso [inclusion_prob()], [ht_total()]
#' @export
joint_prob <- function(data, design, rows = NULL, simulate = FALSE, R = 5000,
                       seed = NULL) {
  if (!is_design(design)) {
    stop("`design` must come from one of the design_*() constructors, not a ",
         class(design)[1], ". See ?designs.", call. = FALSE)
  }
  validate_data(data)
  check_flag(simulate, "simulate")
  rows <- rows %||% seq_len(nrow(data))
  if (!is.numeric(rows) || anyNA(rows) || any(rows < 1) ||
      any(rows > nrow(data))) {
    stop("`rows` must be valid row indices into `data`.", call. = FALSE)
  }
  rows <- as.integer(rows)
  if (isTRUE(simulate)) {
    return(simulate_joint(data, design, rows,
                          check_count(R, "R", allow_zero = FALSE), seed))
  }
  joint_inclusion(design, data, rows)
}

#' Monte Carlo joint inclusion probabilities
#'
#' Counts how often each pair of rows lands in the same sample. Slower and
#' noisier than a closed form, but available for every design -- which makes it
#' the general answer where no formula exists.
#'
#' @noRd
simulate_joint <- function(data, design, rows, R, seed) {
  refuse_simulation(design)
  key <- ".drawn_row_id"
  if (key %in% names(data)) {
    stop("`data` already has a column called `", key,
         "`, which the simulation needs. Rename it.", call. = FALSE)
  }
  tagged <- data
  tagged[[key]] <- seq_len(nrow(data))
  k <- length(rows)
  pos <- match(seq_len(nrow(data)), rows)   # NA for rows we are not tracking

  counts <- matrix(0L, k, k)
  with_seed(seed, quietly_empty({
    for (i in seq_len(R)) {
      s <- draw_design(design, tagged)
      ids <- unique(if (is.data.frame(s)) s[[key]] else unlist(s))
      hit <- pos[ids]
      hit <- hit[!is.na(hit)]
      if (length(hit)) counts[hit, hit] <- counts[hit, hit] + 1L
    }
  }))
  counts / R
}

#' @noRd
no_joint_form <- function(what, alternative) {
  stop(what, " has no closed-form joint inclusion probability.\n", alternative,
       call. = FALSE)
}

#' @noRd
joint_inclusion <- function(design, data, rows) UseMethod("joint_inclusion")

#' @noRd
joint_inclusion.default <- function(design, data, rows) {
  no_joint_form(
    paste0("`", design_type(design), "`"),
    "Only the designs listed in ?joint_prob have one."
  )
}

#' Pairwise matrix from a group id and a per-group (n, N)
#'
#' Rows in the same group are a simple random sample of `n_g` from `N_g`; rows
#' in different groups are selected independently.
#'
#' @noRd
joint_from_groups <- function(group, n_g, size_g, pi_i) {
  k <- length(group)
  out <- outer(pi_i, pi_i)                     # independent case
  same <- outer(group, group, "==")
  if (any(same)) {
    n <- n_g[as.character(group)]
    N <- size_g[as.character(group)]
    within <- outer(seq_len(k), seq_len(k), function(a, b) {
      nn <- n[a]; NN <- N[a]
      ifelse(NN > 1, nn * (nn - 1) / (NN * (NN - 1)), 0)
    })
    out[same] <- within[same]
  }
  diag(out) <- pi_i
  dimnames(out) <- NULL
  out
}

#' @noRd
joint_inclusion.drawn_design_simple <- function(design, data, rows) {
  if (design$replace) {
    no_joint_form("`design_simple(replace = TRUE)`",
                  "Use replace = FALSE, or estimate variance by resampling.")
  }
  N <- nrow(data)
  n <- design$n
  pi_i <- rep(n / N, length(rows))
  out <- matrix(if (N > 1) n * (n - 1) / (N * (N - 1)) else 0,
                length(rows), length(rows))
  diag(out) <- pi_i
  out
}

#' @noRd
joint_inclusion.drawn_design_reservoir <- function(design, data, rows) {
  reach <- reservoir_reach(design, nrow(data))
  out <- joint_inclusion.drawn_design_simple(
    new_design("simple", list(n = min(design$n, reach), replace = FALSE)),
    utils::head(data, reach), rows[rows <= reach]
  )
  if (reach == nrow(data)) return(out)
  # Rows past `max_items` are never read, so they co-occur with nothing.
  full <- matrix(0, length(rows), length(rows))
  seen <- rows <= reach
  full[seen, seen] <- out
  full
}

#' @noRd
joint_inclusion.drawn_design_stratified <- function(design, data, rows) {
  if (design$replace) {
    no_joint_form("`design_stratified(replace = TRUE)`", "Use replace = FALSE.")
  }
  validate_data(data, required_columns = design$strata)
  check_key_columns(data, design$strata, "strata")

  keys <- data[design$strata]
  bad <- Reduce(`|`, lapply(keys, is.na))
  group <- group_key(data, design$strata)
  idx <- split(seq_len(nrow(data))[!bad], group[!bad])
  sizes <- lengths(idx)
  n_alloc <- allocate(design$n, sizes, design$allocation,
                      design$min_per_stratum, cap = TRUE,
                      spread = stratum_spread(design, data, idx))

  pi_all <- exact_inclusion(design, data)
  # A dropped row belongs to no stratum. Give it a label of its own so it
  # co-occurs with nothing, rather than letting NA reach the subscript.
  g <- as.character(group[rows])
  g[bad[rows]] <- paste0("\r__dropped__", seq_len(sum(bad[rows])))
  joint_from_groups(g, n_alloc, sizes, pi_all[rows])
}

#' @noRd
joint_inclusion.drawn_design_temporal <- function(design, data, rows) {
  pi_all <- exact_inclusion(design, data)
  bucket <- temporal_bucket(design, data)
  sizes <- table(bucket[!is.na(bucket)])
  n_g <- pmin(sizes, design$per_interval)

  g <- bucket[rows]
  g[is.na(g)] <- paste0("__out__", seq_len(sum(is.na(g))))  # never co-occur
  joint_from_groups(g, stats::setNames(as.integer(n_g), names(n_g)),
                    stats::setNames(as.integer(sizes), names(sizes)),
                    pi_all[rows])
}

#' @noRd
joint_inclusion.drawn_design_spatial <- function(design, data, rows) {
  require_suggested("sf", "spatial sampling")
  inside <- spatial_inside(design, data)
  N <- sum(inside)
  n <- min(design$n, N)
  pi_all <- exact_inclusion(design, data)
  out <- matrix(if (N > 1) n * (n - 1) / (N * (N - 1)) else 0,
                length(rows), length(rows))
  out[!inside[rows], ] <- 0
  out[, !inside[rows]] <- 0
  diag(out) <- pi_all[rows]
  out
}

#' @noRd
joint_inclusion.drawn_design_cluster <- function(design, data, rows) {
  if (design$balanced) {
    no_joint_form("`design_cluster(balanced = TRUE)`",
                  "Use balanced = FALSE.")
  }
  cl <- count_clusters(design, data)
  A <- cl$total
  a <- design$n_clusters
  lab <- as.character(cl$labels[rows])

  same <- outer(lab, lab, "==")
  same[is.na(same)] <- FALSE
  out <- matrix(if (A > 1) a * (a - 1) / (A * (A - 1)) else 0,
                length(rows), length(rows))
  out[same] <- a / A
  # A row with no cluster -- a dropped key -- is in no sample at all, so it
  # co-occurs with nothing. Without this it inherits the between-cluster rate.
  gone <- is.na(lab)
  if (any(gone)) {
    out[gone, ] <- 0
    out[, gone] <- 0
  }
  pi_all <- exact_inclusion(design, data)
  diag(out) <- pi_all[rows]
  out
}

#' @noRd
joint_inclusion.drawn_design_multistage <- function(design, data, rows) {
  if (design$allocation == "proportional") {
    no_joint_form("`design_multistage(allocation = \"proportional\")`",
                  "Use allocation = \"equal\".")
  }
  if (design$replace) {
    no_joint_form("`design_multistage(replace = TRUE)`", "Use replace = FALSE.")
  }
  if (design$n %% design$n_clusters != 0L) {
    no_joint_form(
      paste0("`design_multistage()` with n = ", design$n, " over ",
             design$n_clusters, " clusters"),
      paste0("The per-cluster take is not constant, so pairs in different\n",
             "clusters have no single joint probability. Choose an `n` that\n",
             "divides evenly by `n_clusters` (here, a multiple of ",
             design$n_clusters, ").")
    )
  }
  cl <- count_clusters(design, data)
  A <- cl$total
  a <- design$n_clusters
  m <- design$n %/% design$n_clusters
  sizes <- table(cl$labels[cl$present])

  lab <- as.character(cl$labels[rows])
  Nh <- as.numeric(sizes[lab])
  same <- outer(lab, lab, "==")
  same[is.na(same)] <- FALSE

  # Different clusters: both clusters selected, then each row within its own.
  out <- outer(m / Nh, m / Nh) * (if (A > 1) a * (a - 1) / (A * (A - 1)) else 0)
  # Same cluster: that cluster selected, then both rows drawn from it.
  win <- outer(seq_along(rows), seq_along(rows), function(i, j) {
    NN <- Nh[i]
    ifelse(NN > 1, (a / A) * (m * (m - 1)) / (NN * (NN - 1)), 0)
  })
  out[same] <- win[same]
  pi_all <- exact_inclusion(design, data)
  diag(out) <- pi_all[rows]
  out
}

#' @noRd
joint_inclusion.drawn_design_weighted <- function(design, data, rows) {
  if (design$method == "poisson") {
    pi_all <- exact_inclusion(design, data)
    out <- outer(pi_all[rows], pi_all[rows])   # independent
    diag(out) <- pi_all[rows]
    return(out)
  }
  if (design$method == "systematic") {
    no_joint_form(
      "`design_weighted(method = \"systematic\")`",
      paste0("Its joint probabilities depend on the order units are visited ",
             "and need a\ndedicated algorithm. `sampling::UPsystematicpi2()` ",
             "computes them.")
    )
  }
  no_joint_form("`design_weighted(method = \"successive\")`",
                "It has no closed-form first-order probability either.")
}

#' @noRd
joint_inclusion.drawn_design_systematic <- function(design, data, rows) {
  if (!is.null(design$start)) {
    stop("A systematic design with a fixed `start` is not a probability ",
         "sample.", call. = FALSE)
  }
  k <- design$interval
  # Residue classes come from the order the design walks, which `order_by`
  # changes. Taking them from frame order instead reports 0 for pairs that
  # always co-occur and 1/k for pairs that never can.
  pos <- systematic_positions(design, data)[rows]
  res <- (pos - 1L) %% k
  out <- matrix(0, length(rows), length(rows))
  same <- outer(res, res, "==")
  same[is.na(same)] <- FALSE
  out[same] <- 1 / k
  diag(out) <- ifelse(is.na(pos), 0, 1 / k)
  out
}

#' @noRd
joint_inclusion.drawn_design_bootstrap <- function(design, data, rows) {
  no_joint_form("`design_bootstrap()`",
                "A bootstrap is not a probability sample of a finite population.")
}


