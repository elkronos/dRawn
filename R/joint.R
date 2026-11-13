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
#'   [design_multistage()]  \tab the two stages multiplied, where the per-cluster take is constant and `size` is not given \cr
#'   [design_certainty()]   \tab `1` between certainty rows; otherwise the other row's own `pi` \cr
#'   [design_reservoir()]   \tab simple random sampling over the first `max_items` rows \cr
#'   [design_temporal()]    \tab within an interval as above; across intervals, independent \cr
#'   [design_spatial()]     \tab simple random sampling inside the region \cr
#'   [design_weighted()]    \tab `"poisson"` only, where rows are independent: `pi_i * pi_j` \cr
#'   [design_systematic()]  \tab `1/interval` for rows sharing a residue class, otherwise **zero**; with `n`, the overlap of the two rows' stretches, from 0 to `n / N`; with `n_starts = m`, `1/interval` for rows sharing a residue class modulo `m * interval`, otherwise `m(m-1) / (mk(mk-1))` with `k` the interval \cr
#' }
#'
#' Systematic sampling is the awkward one. Most pairs can never co-occur, so
#' their joint probability is genuinely 0 and no design-unbiased variance
#' estimator exists; [ht_total()] uses the successive-difference approximation
#' instead, and says so. Its residue classes follow the order the design walks,
#' so `order_by` changes which pairs can co-occur. Given `n`, the design walks
#' at a step of `N / n` rows from a start `u` uniform on (0, 1), each row
#' covering a stretch of length `n / N`, and a pair is drawn together when `u`
#' falls where both stretches, moved into (0, 1), overlap. With a fixed
#' `start` the design is deterministic, and it is refused as in
#' [inclusion_prob()]. With several starts (`n_starts`), rows in different
#' residue classes can be drawn together too, which is what gives a
#' replicated sample its unbiased variance. Within the strata of
#' `design_stratified(within = design_systematic())`, pairs in one stratum
#' get the overlap of their stretches, as with `n`, and pairs in two strata
#' the product of their probabilities.
#'
#' Three designs have joint probabilities that are positive but have no
#' closed form: `design_weighted(method = "systematic")`, which shuffles the
#' rows before walking them so every pair can co-occur, with or without
#' `strata`; [design_multistage()] with `size`, which does the same with the
#' clusters; and [design_spread()]. For all three, [ht_total()] uses a
#' variance approximation that needs only first-order probabilities, and
#' `simulate = TRUE` here estimates the joint probabilities themselves.
#'
#' @section Size:
#' The matrix has a row and a column for each of `rows`, by default every row
#' of `data`, and holds 8 bytes for each pair: 800 MB for 10,000 rows. Past
#' about 1 GB, from about 11,200 rows, `joint_prob()` refuses and gives the
#' size, rather than risk exhausting memory, which can end the R session.
#' [ht_total()] and [ht_mean()] do not build this matrix for any design in
#' the package.
#'
#' @param data A data frame.
#' @param design A design object.
#' @param rows Optional row indices, usually the rows you drew, which are
#'   `s$.row` for a sample `s` from `draw(weights = TRUE)`. The result is then
#'   the submatrix for these rows instead of the full `nrow(data)` square,
#'   which is what makes this usable on a large population. See "Size".
#'   Indices must be whole numbers. A row given more than once is one unit,
#'   so its entry with itself, off the diagonal too, is its own inclusion
#'   probability.
#' @param simulate Estimate the probabilities by repeated draws rather than in
#'   closed form. Works for every probability design, including the ones with no
#'   closed form, at the cost of Monte Carlo error. Refused for
#'   [design_bootstrap()], where every row appears in some replicate and the
#'   count converges to 1 for all of them.
#' @param R Number of simulated draws when `simulate = TRUE`. A simulation
#'   costs `R` draws; see "What a simulation costs" in [inclusion_prob()].
#' @param seed Optional seed for the simulation, applied as in [draw()] (see
#'   "Seeding" there). Checked whether or not `simulate` is `TRUE`.
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
  check_design(design, data, "joint_prob")
  validate_data(data)
  check_flag(simulate, "simulate")
  check_seed(seed)
  warn_region(design, data)
  given <- !is.null(rows)
  rows <- rows %||% seq_len(nrow(data))
  if (!is.numeric(rows) || anyNA(rows) || any(rows < 1) ||
      any(rows > nrow(data))) {
    stop("`rows` must be valid row indices into `data`.", call. = FALSE)
  }
  frac <- rows != trunc(rows)
  if (any(frac)) {
    stop("`rows` must be whole row numbers, not ", format_rows(rows[frac]),
         ". Pass the rows of `data` to pair, such as `s$.row` for a sample ",
         "`s`.", call. = FALSE)
  }
  rows <- as.integer(rows)
  check_pair_size(length(rows), given)
  # A row given twice is one unit, so the pair it makes with itself is in a
  # sample exactly when the row is, with the row's own `pi`. The methods work
  # on distinct rows and the result is spread back out.
  u <- unique(rows)
  out <- if (isTRUE(simulate)) {
    simulate_joint(data, design, u, check_count(R, "R", allow_zero = FALSE),
                   seed)
  } else {
    joint_inclusion(design, data, u)
  }
  if (length(u) < length(rows)) {
    i <- match(rows, u)
    out <- out[i, i, drop = FALSE]
  }
  out
}

#' Refuse a joint matrix too large to hold
#'
#' `k` rows make a `k` x `k` matrix of doubles. Past `pair_limit()` bytes the
#' request is refused before anything is allocated, with the size and the
#' way round it. `given` says whether the caller passed `rows`.
#'
#' @noRd
check_pair_size <- function(k, given) {
  bytes <- 8 * as.numeric(k)^2
  if (bytes <= pair_limit()) return(invisible(NULL))
  k_text <- format(k, big.mark = ",")
  shape <- paste0(k_text, " x ", k_text)
  if (given) {
    stop("`joint_prob()` would build a ", shape, " matrix for the ", k_text,
         " `rows` given, which needs about ", format_gb(bytes),
         " of memory. Pass fewer `rows`.", call. = FALSE)
  }
  stop("`joint_prob()` would build a ", shape, " matrix for every row of ",
       "`data`, which needs about ", format_gb(bytes), " of memory. Pass ",
       "`rows`, such as `rows = s$.row` for the rows of a sample `s` drawn ",
       "with `weights = TRUE`.", call. = FALSE)
}

#' The most memory a matrix over pairs of rows may take, in bytes
#'
#' About 1 GB. Shared by `joint_prob()` and the Sen-Yates-Grundy estimator.
#'
#' @noRd
pair_limit <- function() 1e9

#' A size in bytes, for a message
#' @noRd
format_gb <- function(bytes) {
  paste(format(signif(bytes / 1e9, 2), big.mark = ","), "GB")
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
  design <- main_design(design)
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
  simulate_draws(R, seed, function() {
    ids <- unique(drawn_ids(draw_design(design, tagged), key, design))
    hit <- pos[ids]
    hit <- hit[!is.na(hit)]
    if (length(hit)) counts[hit, hit] <<- counts[hit, hit] + 1L
  })
  counts / R
}

#' Refuse joint probabilities the design does not define
#'
#' `simulate` adds the pointer to a simulation, for the designs where one
#' estimates them: every design but the bootstrap.
#'
#' @noRd
no_joint_form <- function(what, alternative, simulate = TRUE) {
  stop_declined(paste0(
    what, " has no closed-form joint inclusion probability.\n", alternative,
    if (simulate) {
      paste0("\nOr estimate them by Monte Carlo with ",
             "joint_prob(data, design, simulate = TRUE).")
    }))
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

#' A matrix over pairs of rows, built one group's block at a time
#'
#' Pairs of rows in different groups get `between`, times
#' `product[i] * product[j]` when `product` is given. Pairs in the same group
#' get that group's entry of `within`. `group` is each row's position in
#' `within`, or `NA` for a row in no group, which co-occurs with nothing. The
#' diagonal is `on_diag`.
#'
#' The matrix is allocated once and each group's block is assigned in place,
#' so the memory used beyond the result is one block at a time. Comparing
#' every pair of codes with `outer()` held up to twelve matrices of the
#' result's size. `diag<-` is a closure and copies the matrix it is given, so
#' the diagonal is assigned by position instead.
#'
#' @noRd
pair_matrix <- function(group, within, on_diag, between = 0, product = NULL) {
  k <- length(group)
  gone <- which(is.na(group))
  out <- if (is.null(product)) {
    matrix(between, k, k)
  } else {
    product <- as.numeric(product)
    product[gone] <- 0
    tcrossprod(product) * between
  }
  if (length(gone)) {
    out[gone, ] <- 0
    out[, gone] <- 0
  }
  idx <- split(seq_len(k), structure(as.integer(group),
                                     levels = as.character(seq_along(within)),
                                     class = "factor"))
  for (g in which(lengths(idx) > 1L)) {
    out[idx[[g]], idx[[g]]] <- within[[g]]
  }
  out[diag_index(k)] <- on_diag
  out
}

#' Positions of the diagonal of a `k` x `k` matrix
#' @noRd
diag_index <- function(k) seq.int(1, by = k + 1, length.out = k)

#' Pairwise matrix from a group code and a per-group (n, N)
#'
#' Rows in the same group are a simple random sample of `n_g` from `N_g`; rows
#' in different groups are selected independently. `group` is each row's
#' position in `n_g` and `size_g`, or `NA` for a row in no group, which
#' co-occurs with nothing.
#'
#' @noRd
joint_from_groups <- function(group, n_g, size_g, pi_i) {
  n <- as.numeric(n_g)
  N <- as.numeric(size_g)
  within <- ifelse(N > 1, n * (n - 1) / (N * (N - 1)), 0)
  pair_matrix(group, within, pi_i, between = 1, product = pi_i)
}

#' @noRd
joint_inclusion.drawn_design_simple <- function(design, data, rows) {
  if (design$replace) {
    no_joint_form("`design_simple(replace = TRUE)`",
                  "Use replace = FALSE, or estimate variance by resampling.")
  }
  keep <- simple_setup(design, data)
  N <- length(keep)
  n <- design$n
  k <- length(rows)
  out <- matrix(if (N > 1) n * (n - 1) / (N * (N - 1)) else 0, k, k)
  out[diag_index(k)] <- n / N
  # Rows use_reserves() found unavailable are never drawn.
  gone <- !rows %in% keep
  if (any(gone)) {
    out[gone, ] <- 0
    out[, gone] <- 0
  }
  out
}

#' @noRd
joint_inclusion.drawn_design_reservoir <- function(design, data, rows) {
  # A simple random sample from the rows before `max_items`. Rows past it are
  # never read, so they co-occur with nothing.
  N <- reservoir_reach(design, nrow(data))
  n <- min(design$n, N)
  seen <- rows <= N
  pair_matrix(ifelse(seen, 1L, NA_integer_),
              if (N > 1) n * (n - 1) / (N * (N - 1)) else 0,
              ifelse(seen, n / N, 0))
}

#' @noRd
joint_inclusion.drawn_design_stratified <- function(design, data, rows) {
  if (design$replace) {
    no_joint_form("`design_stratified(replace = TRUE)`", "Use replace = FALSE.")
  }
  al <- stratum_allocation(design, data)
  if (walks_strata(design)) return(walk_joint(design, data, rows, al))
  # A dropped row has no stratum code, so it co-occurs with nothing.
  joint_from_groups(al$code[rows], al$n_alloc, al$sizes,
                    stratum_pi(design, al)[rows])
}

#' Joint probabilities of a stratified design that walks each stratum
#'
#' Strata are drawn independently, so a pair in two strata gets the product
#' of its probabilities. A pair in one stratum gets the overlap of its
#' stretches on that stratum's walk, from `fractional_joint()` with the
#' stratum's `n_h` and `N_h`. A row in no stratum co-occurs with nothing.
#'
#' @noRd
walk_joint <- function(design, data, rows, al) {
  pos <- walk_positions(design, data, al$code, length(al$sizes))[rows]
  code <- al$code[rows]
  p <- stratum_pi(design, al)[rows]
  out <- tcrossprod(p)
  for (h in unique(code[!is.na(code)])) {
    i <- which(code == h)
    out[i, i] <- fractional_joint(pos[i], al$n_alloc[h], al$sizes[h])
  }
  out
}

#' @noRd
joint_inclusion.drawn_design_temporal <- function(design, data, rows) {
  bucket <- temporal_setup(design, data)$bucket
  sizes <- tabulate(bucket[!is.na(bucket)])

  # A row outside every interval has no code, so it co-occurs with nothing.
  joint_from_groups(bucket[rows], pmin(sizes, design$per_interval), sizes,
                    temporal_pi(design, bucket)[rows])
}

#' @noRd
joint_inclusion.drawn_design_spatial <- function(design, data, rows) {
  inside <- spatial_setup(design, data)
  N <- sum(inside)
  n <- design$n
  pi_all <- exact_inclusion(design, data)
  # A row outside the region co-occurs with nothing.
  pair_matrix(ifelse(inside[rows], 1L, NA_integer_),
              if (N > 1) n * (n - 1) / (N * (N - 1)) else 0, pi_all[rows])
}

#' @noRd
joint_inclusion.drawn_design_cluster <- function(design, data, rows) {
  cl <- cluster_setup(design, data)
  if (design$balanced) {
    no_joint_form("`design_cluster(balanced = TRUE)`",
                  "Use balanced = FALSE.")
  }
  A <- cl$total
  a <- design$n_clusters
  pi_all <- exact_inclusion(design, data)
  # Same cluster: that cluster selected. Different clusters: both selected. A
  # row with no cluster, dropped under `na_rm = TRUE`, is in no sample, so it
  # co-occurs with nothing.
  pair_matrix(cl$code[rows], rep(a / A, A), pi_all[rows],
              between = if (A > 1) a * (a - 1) / (A * (A - 1)) else 0)
}

#' @noRd
joint_inclusion.drawn_design_multistage <- function(design, data, rows) {
  cl <- multistage_setup(design, data)
  if (!is.null(design$size)) {
    no_joint_form(
      "`design_multistage(size = )`",
      paste0("Its clusters are shuffled before the systematic PPS walk, so ",
             "every pair of clusters\ncan be selected together, but their ",
             "joint probabilities average over every ordering\nand have no ",
             "closed form. ht_total() uses Deville's approximation for the ",
             "clusters.")
    )
  }
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
  A <- cl$total
  a <- design$n_clusters
  m <- design$n %/% design$n_clusters
  sizes <- as.numeric(tabulate(cl$code, A))
  code <- cl$code[rows]
  pi_all <- exact_inclusion(design, data)

  # Different clusters: both clusters selected, then each row within its own.
  # Same cluster: that cluster selected, then both rows drawn from it. A row
  # with no cluster co-occurs with nothing.
  pair_matrix(code,
              ifelse(sizes > 1,
                     (a / A) * (m * (m - 1)) / (sizes * (sizes - 1)), 0),
              pi_all[rows],
              between = if (A > 1) a * (a - 1) / (A * (A - 1)) else 0,
              product = m / sizes[code])
}

#' @noRd
joint_inclusion.drawn_design_weighted <- function(design, data, rows) {
  weighted_setup(design, data)
  if (design$method == "poisson") {
    pi_all <- exact_inclusion(design, data)
    p <- pi_all[rows]
    out <- tcrossprod(p)                       # independent
    out[diag_index(length(rows))] <- p
    return(out)
  }
  if (design$method == "systematic") {
    no_joint_form(
      "`design_weighted(method = \"systematic\")`",
      paste0("The rows are shuffled before the systematic walk, so every pair ",
             "can co-occur,\nbut the joint probabilities average over every ",
             "ordering and have no closed form.\nht_total() uses Deville's ",
             "approximation.")
    )
  }
  no_joint_form("`design_weighted(method = \"successive\")`",
                "It has no closed-form first-order probability either.")
}

#' @noRd
joint_inclusion.drawn_design_systematic <- function(design, data, rows) {
  refuse_fixed_start(design)
  # Residue classes come from the order the design walks, which `order_by`
  # changes. Taking them from frame order instead reports 0 for pairs that
  # always co-occur and 1/k for pairs that never can.
  walk <- systematic_positions(design, data)
  pos <- walk[rows]
  if (!is.null(design$n)) {
    return(fractional_joint(pos, design$n, sum(!is.na(walk))))
  }
  # With `m` starts drawn from 1:(m * k), a pair in one residue class modulo
  # m * k is taken when its start is among them, m times in m * k, and a
  # pair in two classes when both starts are. One start is the case m = 1.
  # Only the classes of `rows` get a block.
  k <- design$interval
  m <- systematic_starts(design)
  step <- m * k
  cls <- as.integer((pos - 1L) %% step)
  live <- unique(cls[!is.na(cls)])
  pair_matrix(match(cls, live), rep(1 / k, length(live)),
              ifelse(is.na(pos), 0, 1 / k),
              between = if (m > 1L) (m - 1) / (k * (step - 1)) else 0)
}

#' Joint probabilities of a systematic sample of fixed size `n`
#'
#' Walk position `p` is taken when the start `u` lands in its stretch moved
#' into (0, 1): an arc of length `n / N` beginning at `(p - 1) * n / N`
#' modulo 1. Two positions are taken together when `u` lands where their arcs
#' overlap. For positions `d` apart, the arcs begin `gap / N` apart, where
#' `gap` is `(d * n) %% N`, and overlap by
#' `(max(n - gap, 0) + max(gap + n - N, 0)) / N`. The arithmetic is on whole
#' numbers, so when `n` divides `N` the result is `n / N` within a residue
#' class and 0 between them, as for `interval = N / n`.
#'
#' @noRd
fractional_joint <- function(pos, n, N) {
  k <- length(pos)
  out <- matrix(0, k, k)
  ok <- which(!is.na(pos))
  if (length(ok) == 0L) return(out)
  # The overlap depends only on the distance along the walk, modulo N.
  gap <- (as.double(seq.int(0L, N - 1L)) * n) %% N
  overlap <- (pmax(n - gap, 0) + pmax(gap + n - N, 0)) / N
  p <- pos[ok]
  for (j in seq_along(ok)) {
    out[ok, ok[j]] <- overlap[(p - p[j]) %% N + 1L]
  }
  out
}

#' @noRd
joint_inclusion.drawn_design_bootstrap <- function(design, data, rows) {
  no_joint_form("`design_bootstrap()`",
                "A bootstrap is not a probability sample of a finite population.",
                simulate = FALSE)
}


