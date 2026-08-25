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
  with_seed(seed, {
    for (i in seq_len(R)) {
      s <- draw_design(design, tagged)
      ids <- unique(if (is.data.frame(s)) s[[key]] else unlist(s))
      hit <- pos[ids]
      hit <- hit[!is.na(hit)]
      if (length(hit)) counts[hit, hit] <- counts[hit, hit] + 1L
    }
  })
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


# ---- Horvitz-Thompson estimation -------------------------------------------

#' Estimate a population total from a sample
#'
#' Forms the Horvitz-Thompson total `sum(y / pi)` and, where the design allows
#' it, a design-unbiased variance and confidence interval.
#'
#' @section Variance:
#' The estimator is chosen to match how the design actually randomises:
#'
#' * **Fixed-size designs** use the Sen-Yates-Grundy estimator, which is
#'   non-negative more often than the general Horvitz-Thompson form and is the
#'   usual choice.
#' * **Poisson sampling** has a random size and independent rows, so the
#'   independent-units form `sum((1 - pi) / pi^2 * y^2)` is used instead.
#' * **Cluster designs** take whole clusters, so the number of *rows* is random
#'   whenever the clusters differ in size. The cluster is the sampling unit, and
#'   the estimator is applied at that level — algebraically the same thing as
#'   the delete-a-cluster jackknife.
#' * **Certainty designs** hand the problem to `rest` over the rows below the
#'   threshold, since the certainty rows are in every possible sample and
#'   contribute nothing to the variance.
#'
#' Sen-Yates-Grundy can still return a negative number on an unlucky sample.
#' That is a failure of the estimator rather than a variance, so it is reported
#' as one: `variance` is `NA`, `note` says what happened, and `variance = "auto"`
#' falls through to the jackknife.
#'
#' A variance needs joint inclusion probabilities, and not every design has
#' them — see [joint_prob()]. Where they are unavailable the estimate is still
#' returned, with `variance` as `NA` and a note saying why. Systematic sampling
#' is the notable case: most pairs of rows can never co-occur, so no
#' design-unbiased variance estimator exists at all.
#'
#' @param sample A data frame returned by [draw()] with `weights = TRUE`.
#' @param y The variable to total: a column name, or a numeric vector as long as
#'   `sample`.
#' @param variance How to compute it. `"auto"` uses the analytic estimator when
#'   the design has one and falls back to the jackknife when it does not;
#'   `"analytic"` insists on the analytic form, returning `NA` with the reason
#'   in `note` rather than falling back; `"jackknife"` always resamples;
#'   `"none"` skips it. The result reports which was used in `method` — and
#'   reports `"none"` when neither could produce a figure, rather than naming a
#'   method that declined.
#' @param level Confidence level for the interval.
#'
#' @return A list with a `print()` method, holding:
#'   \describe{
#'     \item{`total`}{The Horvitz-Thompson total, `sum(y / pi)`.}
#'     \item{`variance`, `se`, `ci`, `level`}{Its estimated variance, standard
#'       error and confidence interval. `NA` where the design supports none.}
#'     \item{`n`, `design`}{Rows used, and the design's type.}
#'     \item{`deff`}{The design effect — see [deff()].}
#'     \item{`method`}{`"analytic"`, `"jackknife"` or `"none"`.}
#'     \item{`note`}{Why a variance is missing, or which fallback was taken.
#'       `NULL` when the analytic estimator applied cleanly.}
#'   }
#'
#' @examples
#' set.seed(1)
#' pop <- data.frame(
#'   id = 1:200,
#'   site = rep(c("a", "b"), times = c(150, 50)),
#'   spend = round(stats::runif(200, 10, 500))
#' )
#'
#' s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
#' ht_total(s, "spend")
#'
#' sum(pop$spend)   # the truth
#'
#' @seealso [joint_prob()], [inclusion_prob()]
#' @export
ht_total <- function(sample, y, variance = c("auto", "analytic", "jackknife",
                                            "none"), level = 0.95) {
  variance <- match.arg(variance)
  if (!is.data.frame(sample)) {
    stop("`sample` must be a data frame returned by draw().", call. = FALSE)
  }
  parts <- ht_prepare(sample, y, level)
  design <- parts$design; rows <- parts$rows; pop <- parts$pop
  yv <- parts$y; pi_i <- parts$pi
  total <- sum(yv / pi_i)

  var_out <- ht_variance_dispatch(design, parts$sample, pop, rows, yv, pi_i,
                                  variance)

  finish_estimate(total, var_out, level, nrow(sample), design, yv, pi_i,
                  nrow(pop), what = "total", class = "drawn_ht")
}

#' Delete-a-group jackknife variance
#'
#' Works from the sample alone, which is what makes it available where the
#' analytic form is not. Groups are the primary sampling units: whole clusters
#' where the design has them, otherwise individual rows. Deleting a group,
#' inflating the surviving weights to compensate, and looking at how far the
#' estimate moves is a direct measure of how much the estimate depended on
#' which groups were drawn.
#'
#' Two designs are declined rather than approximated. A systematic sample has a
#' single primary sampling unit -- the random start -- so deleting rows does not
#' reflect its randomness at all. Poisson sampling has an exact variance and a
#' random size; deleting rows from it understates by a factor of three.
#'
#' @noRd
jackknife_variance <- function(design, sample, pop, y, pi_i) {
  refusal <- jackknife_refusal(design)
  if (!is.null(refusal)) {
    return(list(variance = NA_real_, method = "jackknife", note = refusal))
  }
  # Certainty rows are in every possible sample. Deleting one and inflating the
  # rest invents variance the design does not have -- four times too much on a
  # sample that is a third certainty -- so hand the whole problem to `rest`.
  if (inherits(design, "drawn_design_certainty")) {
    keep <- pi_i < 1
    if (!any(keep)) {
      return(list(variance = 0, method = "jackknife", note = paste0(
        "Every sampled row was taken with certainty, so the total is exact ",
        "rather than estimated.")))
    }
    sp <- certainty_split(design, pop)
    return(jackknife_variance(design$rest, sample[keep, , drop = FALSE],
                              sp$data[sp$rest, , drop = FALSE],
                              y[keep], pi_i[keep]))
  }
  grp <- jackknife_groups(design, sample)
  m <- length(unique(grp))
  if (m < 2L) {
    return(list(variance = NA_real_, method = "jackknife",
                note = paste0("The jackknife needs at least two primary ",
                              "sampling units; this sample has ", m, ".")))
  }
  base <- y / pi_i
  ug <- unique(grp)
  reps <- vapply(ug, function(g) sum(base[grp != g]) * m / (m - 1), numeric(1))

  # Without a finite population correction the jackknife treats the frame as
  # infinite and overstates the variance -- by a factor of four when a quarter
  # of the clusters were taken.
  fpc <- jackknife_fpc(design, pop, m)
  v <- fpc * ((m - 1) / m) * sum((reps - mean(reps))^2)
  list(variance = v, method = "jackknife",
       note = paste0("Jackknife over ", m, " primary sampling unit(s)",
                     if (fpc < 1) paste0(", with a finite population ",
                                         "correction of ", signif(fpc, 3))
                     else "", "."))
}

#' Designs the jackknife should decline rather than approximate
#'
#' @noRd
jackknife_refusal <- function(design) {
  if (inherits(design, "drawn_design_systematic")) {
    return(paste0("A systematic sample has one primary sampling unit -- the ",
                  "random start -- so deleting rows does not reflect the ",
                  "design's randomness and the jackknife is not valid here."))
  }
  if (inherits(design, "drawn_design_weighted") && design$method == "poisson") {
    return(paste0("Poisson sampling selects rows independently and has an ",
                  "exact variance, sum((1 - pi) / pi^2 * y^2). Deleting rows ",
                  "from a sample whose size is itself random understates it ",
                  "by a factor of three. Use variance = \"analytic\"."))
  }
  NULL
}

#' The share of primary sampling units NOT taken
#'
#' Applied for single-stage designs, where deleting a sampling unit accounts for
#' all the randomness there is. That includes probability-proportional-to-size
#' selection, which is still a fixed-size draw without replacement: measured
#' against the empirical sampling variance over four frames and two response
#' variables, the corrected estimator sits between 0.86 and 1.11 of the truth
#' where the uncorrected one runs from 0.98 to 1.77, overstating badly once the
#' sampling fraction gets large.
#'
#' Deliberately not applied to multi-stage designs:
#' the delete-a-cluster jackknife sees only between-cluster variation, so the
#' second stage is already missing, and correcting for the first stage on top of
#' that understates the total. Leaving it out is the ultimate-cluster
#' approximation, which errs conservative -- the usual choice in survey practice.
#'
#' @noRd
jackknife_fpc <- function(design, pop, m) {
  if (inherits(design, "drawn_design_multistage")) return(1)
  total <- tryCatch({
    if (inherits(design, "drawn_design_cluster")) {
      count_clusters(design, pop)$total
    } else {
      nrow(pop)
    }
  }, error = function(e) NA_real_)
  if (!is.finite(total) || total <= m) return(1)
  1 - m / total
}

#' @noRd
jackknife_groups <- function(design, sample) {
  col <- if (inherits(design, c("drawn_design_cluster",
                                "drawn_design_multistage"))) {
    design$clusters
  } else if (inherits(design, "drawn_design_stratified")) {
    NULL   # rows within strata are the sampling units
  } else {
    NULL
  }
  if (!is.null(col) && col %in% names(sample)) {
    as.character(sample[[col]])
  } else {
    as.character(seq_len(nrow(sample)))
  }
}

#' @noRd
ht_variance <- function(design, data, rows, y, pi_i) UseMethod("ht_variance")

#' Variance of a single-stage cluster total
#'
#' Whole clusters are taken, so the number of *rows* is random whenever the
#' clusters differ in size. Sen-Yates-Grundy assumes a fixed size, and applied
#' row by row here it understates the variance badly -- by a factor of five on a
#' frame whose clusters vary from 2 to 10 rows -- and can return a negative
#' number or a zero-width interval.
#'
#' The cluster is the sampling unit, so the estimator belongs at that level: the
#' clusters are a simple random sample of `a` from `A`, and the quantity summed
#' over them is each cluster's contribution to the total. That is the textbook
#' form, and it is algebraically identical to the delete-a-cluster jackknife.
#'
#' @noRd
ht_variance.drawn_design_cluster <- function(design, data, rows, y, pi_i) {
  cl <- count_clusters(design, data)
  A <- cl$total
  a <- design$n_clusters
  lab <- as.character(cl$labels[rows])
  if (anyNA(lab)) {
    return(list(variance = NA_real_,
                note = "Some sampled rows have no cluster label."))
  }
  u <- vapply(split(y / pi_i, lab), sum, numeric(1))
  m <- length(u)
  if (m < 2L) {
    return(list(variance = NA_real_, note = paste0(
      "A variance needs at least two clusters; this sample has ", m, ".")))
  }
  fpc <- if (is.finite(A) && A > m) 1 - m / A else 0
  list(variance = fpc * m * stats::var(u), note = NULL)
}

#' @noRd
ht_variance.default <- function(design, data, rows, y, pi_i) {
  if (inherits(design, "drawn_design_systematic")) {
    stop("Systematic sampling has no design-unbiased variance estimator: most ",
         "pairs of rows can never appear together, so their joint inclusion ",
         "probability is 0. Repeat the draw with different starts, or treat ",
         "the sample as simple random, which is conservative.", call. = FALSE)
  }

  # Poisson: units are independent, so no joint matrix is needed.
  if (inherits(design, "drawn_design_weighted") && design$method == "poisson") {
    return(list(variance = sum((1 - pi_i) / pi_i^2 * y^2), note = NULL))
  }

  pij <- joint_inclusion(design, data, rows)
  k <- length(rows)
  if (k < 2L) {
    return(list(variance = NA_real_,
                note = "A variance needs at least two sampled rows."))
  }

  # Sen-Yates-Grundy, for fixed-size designs.
  yk <- y / pi_i
  d <- outer(yk, yk, "-")^2
  num <- outer(pi_i, pi_i) - pij
  w <- num / pij
  w[!is.finite(w)] <- 0                 # a zero joint probability contributes nothing
  diag(w) <- 0
  list(variance = 0.5 * sum(w * d), note = NULL)
}

#' @export
print.drawn_ht <- function(x, ...) {
  print_estimate(x, "total", "Horvitz-Thompson total")
}
