#' Weighted sampling
#'
#' Draws rows with probability governed by a weights column. The weights are
#' used as supplied: `prob=` normalises by their sum, so rescaling them changes
#' nothing and the design does not offer it.
#'
#' @section Choosing a method:
#' The three methods differ in what the weights actually control, which decides
#' whether the sample can be used for estimation.
#'
#' \describe{
#'   \item{`"successive"`}{The default, and what [base::sample()] does. The
#'     weights govern each sequential draw, not the probability that a row ends
#'     up in the sample, and the realised inclusion probabilities are **not**
#'     proportional to the weights. How far off depends on the weights: with a
#'     single dominant unit the two nearly coincide (about 1% apart), but with a
#'     moderate spread the gap is large — weights `1:10` at `n = 5` give
#'     inclusion probabilities up to about 36% away from proportional, by
#'     exact enumeration of the ordered draws. Treating them as
#'     if they were proportional biases a Horvitz-Thompson total by around 1% at
#'     `n = 15` out of 30 in simulation, small in size but unmistakable in sign.
#'     There is no closed form for `pi`, so [inclusion_prob()] refuses to give
#'     one. Fine when you want a weighted selection; wrong as the basis for an
#'     estimate.}
#'   \item{`"systematic"`}{Randomised systematic probability-proportional-to-
#'     size (Madow 1949; Hartley and Rao 1962). Shuffles the rows, then walks
#'     the cumulative weights with a fixed step from a random start, giving
#'     `pi_i = n * p_i` exactly. Rows heavy enough that `n * p_i > 1` are taken
#'     with certainty and the rest rescaled, repeatedly, until every
#'     probability is valid. Fixed sample size. The shuffle means every pair of
#'     rows can appear together, but their joint probabilities have no closed
#'     form, so [ht_total()] uses Deville's approximation for the variance,
#'     which needs only the first-order probabilities.}
#'   \item{`"poisson"`}{Each row is included independently with probability
#'     `pi_i = n * p_i`. Rows heavy enough that `n * p_i > 1` are taken with
#'     certainty and the rest rescaled, repeatedly, as for `"systematic"`, so
#'     the probabilities still sum to `n`. Inclusion probabilities are exactly
#'     proportional to size among the rows below certainty, and every pair can
#'     co-occur, at the cost of a **random sample size** averaging `n`.}
#' }
#'
#' `replace = TRUE` is with-replacement PPS and applies only to
#' `method = "successive"`.
#'
#' @section Strata:
#' With `strata`, `n` is split over the strata by `allocation`, as
#' [design_stratified()] splits it over the rows in each stratum, and each
#' stratum is then sampled in proportion to size on its own: a row of stratum
#' `h` gets `n_h * w_i / sum(w)`, the sum running over the stratum, with rows
#' past certainty taken with certainty within it. Under `"systematic"` each
#' stratum is shuffled and walked from a start of its own, in key order;
#' under `"poisson"` every row is included independently with the probability
#' its stratum gives it, so the size of each stratum's sample is random.
#' `"successive"` has no closed-form inclusion probability and is refused with
#' `strata`. To split `n` in proportion to each stratum's total size rather
#' than its rows, work out the counts and give them as `allocation`.
#'
#' The strata are sampled independently. [ht_total()] and [ht_mean()] add up
#' Deville's approximation over the strata for `"systematic"`, and the
#' Poisson variance needs no change. The degrees of freedom are the sampled
#' rows below certainty less the strata that hold them. A stratum with one
#' sampled row below certainty has nothing to compare it with, and the
#' variance is refused with the stratum named; draw at least two in each.
#'
#' In simulations of 2,000 draws each from frames of 600 rows in four strata,
#' drawing 24 or 60 rows by `"systematic"`, the variance averaged 0.94 to
#' 1.04 times the true variance of the total, and 95% intervals covered 88%
#' to 96% of the time. The low coverage came where `y / pi` was skewed, so
#' that samples with a low total also had a low estimated variance; the same
#' frames without strata covered 92% to 94%.
#'
#' [as_svydesign()] hands the strata to `survey`, with Brewer's approximation
#' within each, [allocation_table()] shows the split before anything is
#' drawn, and [sample_summary()] counts the rows drawn from each stratum.
#'
#' Drawing each stratum with a design of its own and binding the samples with
#' [rbind()] does not give a sample to estimate from. The bound sample keeps
#' the first draw's frame and design, so [ht_total()] refuses it; use
#' `strata` instead.
#'
#' @param weights Column of positive, finite, numeric weights. A weight of 0
#'   is refused under every method, not given a probability of 0; drop such
#'   rows from the frame first, since a design that selects in proportion to
#'   size could never select them.
#' @param n Number of rows to draw, the total across `strata` when given, or
#'   a plan from [plan_size()].
#'   Under `method = "poisson"` this is the expected number, not a guarantee.
#'   With counts in `allocation` it can be left out and is their sum; when
#'   given, it must equal that sum.
#' @param replace Sample with replacement? Only available for
#'   `method = "successive"`.
#' @param method One of `"successive"`, `"systematic"` or `"poisson"`. See
#'   "Choosing a method".
#' @param na_rm Drop rows whose weight or stratum key is `NA` instead of
#'   raising an error.
#' @param strata Optional column name(s) defining strata, to sample in
#'   proportion to size within each one. Several columns are
#'   cross-classified. Needs `method = "systematic"` or `"poisson"`. See
#'   "Strata".
#' @param allocation How `n` is split across `strata`: `"proportional"` to
#'   the number of rows in each stratum, or `"equal"`, by the
#'   largest-remainder method of [design_stratified()], with a stratum that
#'   would get more rows than it holds taken whole; or the count for each
#'   stratum, as a named vector such as `c(north = 10, south = 5)` when
#'   `strata` is one column, or a data frame of the `strata` columns and `n`.
#'   Counts follow the rules of [design_stratified()]: every stratum of the
#'   frame needs one, a count for a stratum the frame does not have is
#'   refused, and no count may exceed the rows in its stratum. Only with
#'   `strata`.
#'
#' @return A design object, for use with [draw()].
#'
#' @references
#' Madow, W. G. (1949). On the theory of systematic sampling, II. *Annals of
#' Mathematical Statistics*, 20, 333–354.
#'
#' Hartley, H. O. and Rao, J. N. K. (1962). Sampling with unequal
#' probabilities and without replacement. *Annals of Mathematical Statistics*,
#' 33, 350–374.
#'
#' Brewer, K. R. W. and Hanif, M. (1983). *Sampling with Unequal
#' Probabilities*. Springer.
#'
#' Tillé, Y. (2006). *Sampling Algorithms*. Springer.
#'
#' @examples
#' df <- data.frame(id = 1:20, w = 1:20)
#' draw(df, design_weighted("w", n = 5), seed = 1)
#'
#' # Inclusion probabilities proportional to size, and checkable
#' d <- design_weighted("w", n = 5, method = "systematic")
#' round(inclusion_prob(df, d), 3)
#'
#' # Which is what makes an unbiased total possible
#' s <- draw(df, d, seed = 1, weights = TRUE)
#' sum(s$w * s$.weight)   # estimates sum(df$w) = 210
#'
#' # In proportion to size within each of two regions
#' df$region <- rep(c("east", "west"), each = 10)
#' by_region <- design_weighted("w", n = 8, method = "systematic",
#'                              strata = "region")
#' allocation_table(df, by_region)
#' round(inclusion_prob(df, by_region), 3)
#'
#' @family designs
#' @seealso [draw()], [inclusion_prob()], [certainty_set()] for the rows
#'   taken with certainty, [design_spread()] to select in
#'   proportion to size *and* spread the sample across a map. The `sampling`
#'   package implements further fixed-size unequal-probability algorithms
#'   (Brewer, Sampford, Tillé, maximum entropy) as vector functions.
#' @export
design_weighted <- function(weights, n, replace = FALSE,
                            method = c("successive", "systematic", "poisson"),
                            na_rm = FALSE, strata = NULL,
                            allocation = c("proportional", "equal")) {
  if (!missing(n)) n <- plan_n(n)
  method <- check_choice(method, c("successive", "systematic", "poisson"),
                         "method")
  replace <- check_flag(replace, "replace")
  if (replace && method != "successive") {
    stop("`replace = TRUE` is with-replacement PPS and cannot be combined with ",
         "method = \"", method, "\". Use method = \"successive\", or set ",
         "replace = FALSE.", call. = FALSE)
  }
  if (is.null(strata)) {
    if (!missing(allocation) && !is.null(allocation)) {
      stop("`allocation` splits `n` over strata, and `strata` is NULL. Name ",
           "the strata columns in `strata`, or leave `allocation` out.",
           call. = FALSE)
    }
    allocation <- NULL
    n <- check_count(n, "n")
  } else {
    strata <- check_columns(strata, "strata")
    if (method == "successive") {
      stop("`strata` cannot be combined with method = \"successive\", whose ",
           "inclusion probabilities have no closed form, so a total could not ",
           "be estimated from the strata. Use method = \"systematic\" or ",
           "\"poisson\".", call. = FALSE)
    }
    if (is.numeric(allocation) || is.data.frame(allocation)) {
      allocation <- check_stratum_counts(allocation, strata)
      n <- counts_n(if (!missing(n)) n, allocation)
    } else {
      if (!is.null(allocation) && !is.character(allocation)) {
        stop("`allocation` must be \"proportional\", \"equal\" or the count ",
             "for each stratum, not ", format_bad(allocation), ".",
             call. = FALSE)
      }
      allocation <- check_choice(allocation, c("proportional", "equal"),
                                 "allocation")
      if (missing(n) || is.null(n)) {
        stop("`n` is missing. Give the total number of rows to draw, or the ",
             "count for each stratum in `allocation`.", call. = FALSE)
      }
      n <- check_count(n, "n")
    }
  }
  new_design("weighted", list(
    weights = check_columns(weights, "weights", 1L),
    n = n,
    replace = replace,
    method = method,
    na_rm = check_flag(na_rm, "na_rm"),
    strata = strata,
    allocation = allocation
  ))
}

#' Weights a weighted design can use, and a size the frame can supply
#'
#' The checks the draw and the probabilities share, so [inclusion_prob()]
#' refuses what [draw()] refuses, with the same message: weights that are not
#' numeric, missing under `na_rm = FALSE`, missing in every row, not finite or
#' not positive, and an `n` above the rows left. Returns `keep`, the rows the
#' design uses, and `w`, their weights.
#'
#' With `strata`, a row with a missing stratum key is refused or dropped as a
#' missing weight is, and `n` is split over the strata of the rows kept by
#' `stratum_allocation()`. The result then also holds `code`, each kept
#' row's stratum as a position in `labels`, `sizes` and `n_alloc`, and `idx`,
#' the positions in `keep` of each stratum's rows, in key order.
#'
#' @noRd
weighted_setup <- function(design, data) {
  col <- design$weights
  strata <- design$strata
  validate_data(data, required_columns = c(col, strata))
  check_key_columns(data, col, "weights")
  if (!is.null(strata)) check_key_columns(data, strata, "strata")
  w <- data[[col]]
  check_weights(w, col)
  what <- paste0("a missing weight in `", col, "`")
  bad <- check_na_policy(is.na(w), design$na_rm, what)
  if (all(bad)) stop("Every row has ", what, ".", call. = FALSE)
  if (!is.null(strata)) return(weighted_strata(design, data, w, bad))
  keep <- which(!bad)
  n_rows <- length(keep)
  if (design$method == "successive") {
    check_draw_size(design$n, n_rows, design$replace)
  } else if (design$n > n_rows) {
    stop("`n` (", design$n, ") cannot exceed the number of rows (", n_rows,
         ") for method = \"", design$method, "\".", call. = FALSE)
  }
  list(keep = keep, w = w[keep])
}

#' The strata of a weighted design, and its split of `n` over them
#'
#' `bad` marks the rows dropped for a missing weight. The allocation runs
#' over the rows left, so a stratum's size counts only rows the design can
#' draw.
#'
#' @noRd
weighted_strata <- function(design, data, w, bad) {
  strata <- design$strata
  key <- group_key(data, strata)
  gone <- check_na_policy(is.na(key) & !bad, design$na_rm,
                          "a missing stratum key")
  keep <- which(!bad & !gone)
  if (!length(keep)) {
    stop("Every row has a missing weight in `", design$weights, "` or a ",
         "missing stratum key.", call. = FALSE)
  }
  refuse_emptied_counts(design, data, key, keep,
                        paste0("a missing weight in `", design$weights, "`"))
  # These methods have no `replace`, so a count above its stratum's rows is
  # refused without offering one.
  split_by <- design
  split_by$replace <- NULL
  al <- stratum_allocation(split_by, table_frame(
    lapply(strata, function(col) data[[col]][keep]), strata, list()))
  list(keep = keep, w = w[keep], code = al$code, labels = al$labels,
       idx = al$idx, sizes = al$sizes, n_alloc = al$n_alloc)
}

#' Each kept row's inclusion probability under a PPS method
#'
#' `st` is `weighted_setup()`. Without strata, `pps_pi()` over all the rows;
#' with them, `pps_pi()` within each stratum, for its share of `n`.
#'
#' @noRd
weighted_pi <- function(design, st) {
  if (is.null(design$strata)) return(pps_pi(st$w, design$n))
  pi <- numeric(length(st$keep))
  for (h in seq_along(st$idx)) {
    i <- st$idx[[h]]
    pi[i] <- pps_pi(st$w[i], st$n_alloc[h])
  }
  pi
}

#' @export
draw_design.drawn_design_weighted <- function(design, data) {
  st <- weighted_setup(design, data)
  if (length(st$keep) < nrow(data)) data <- reindex(data, st$keep)
  w <- st$w
  n_rows <- nrow(data)

  if (design$method == "successive") {
    idx <- sample.int(n_rows, size = design$n, replace = design$replace,
                      prob = finite_total(w))
    return(reindex(data, idx))
  }

  if (!is.null(design$strata)) {
    warn_empty_groups(st$n_alloc, st$sizes, st$labels, fix = paste0(
      if (has_counts(design)) "Give " else "Raise `n`, or give ",
      "each stratum a count of at least 1 in `allocation` (or 2, so each ",
      "stratum also supports a variance)."))
  }
  pi <- weighted_pi(design, st)

  idx <- if (design$method == "poisson") {
    which(stats::runif(n_rows) < pi)
  } else if (!is.null(design$strata)) {
    # Each stratum is shuffled and walked on its own, in key order.
    unlist(lapply(st$idx, function(i) i[systematic_pps(pi[i])]),
           use.names = FALSE)
  } else {
    systematic_pps(pi)
  }

  reindex(data, idx, sort = TRUE)
}

#' Inclusion probabilities for a fixed-size PPS design
#'
#' `n * p_i` can exceed 1 for a heavy unit, which is not a probability. The
#' standard remedy is to take such units with certainty, remove them, and
#' redistribute the remaining sample size over the rest -- repeating, because
#' removing units raises everyone else's share.
#'
#' A share that rounds to slightly below 1 is certain too, and gets exactly 1,
#' so every check on `.prob` sees the same certainty units as the draw.
#'
#' @noRd
pps_pi <- function(w, n) {
  w <- finite_total(w)
  pi <- rep(0, length(w))
  certain <- rep(FALSE, length(w))

  repeat {
    free <- !certain
    n_free <- n - sum(certain)
    if (n_free <= 0 || !any(free)) break
    # Positive weights can still divide to 0 against the largest, when they
    # span more orders of magnitude than a double holds; see the check below.
    if (!(sum(w[free]) > 0)) break
    p <- w[free] / sum(w[free])
    cand <- n_free * p
    over <- is_certain(cand)
    if (!any(over)) {
      pi[free] <- cand
      break
    }
    idx_free <- which(free)
    certain[idx_free[over]] <- TRUE
  }

  pi[certain] <- 1
  if (!isTRUE(abs(sum(pi) - n) <= 1e-6 * max(n, 1))) {
    stop("The size measure gives inclusion probabilities summing to ",
         format(sum(pi)), ", where `n` is ", n, ", because its values span ",
         "more orders of magnitude than double precision holds. Take the ",
         "largest units with design_certainty() and sample the rest.",
         call. = FALSE)
  }
  pi
}

#' Weights whose sum is finite
#'
#' A probability proportional to size does not depend on the scale of the
#' sizes. Positive weights near `1e306` sum to `Inf`, which made every
#' probability 0, so such weights are divided by the largest of them. Weights
#' with a finite sum are returned as they are, which keeps their
#' probabilities, and the rows a seed selects, unchanged.
#'
#' @noRd
finite_total <- function(w) {
  if (length(w) && !is.finite(sum(w))) w <- w / max(w)
  w
}

#' Systematic PPS selection given inclusion probabilities
#'
#' Walks the cumulative probabilities with a fixed step of 1 from a random start
#' in [0, 1). Yields exactly `sum(pi)` units, each with probability `pi_i`.
#'
#' @noRd
systematic_pps <- function(pi) {
  n <- as.integer(round(sum(pi)))
  if (n <= 0) return(integer(0))
  # Randomise the order first: systematic PPS is order-dependent, and leaving
  # the data's own order in place would make the design depend on it.
  ord <- sample.int(length(pi))
  cum <- cumsum(pi[ord])
  u <- stats::runif(1)
  picks <- findInterval(u + seq.int(0, n - 1), cum) + 1L
  picks <- picks[picks <= length(pi)]
  sort(ord[unique(picks)])
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_weighted <- function(design, data) {
  st <- weighted_setup(design, data)
  if (design$method == "successive") {
    no_closed_form(
      "`design_weighted(method = \"successive\")`",
      paste0("Its weights govern each successive draw, not each row's chance ",
             "of ending up in\nthe sample. Use method = \"systematic\" or ",
             "\"poisson\" for a design whose inclusion\nprobabilities really ",
             "are proportional to the weights.")
    )
  }
  out <- numeric(nrow(data))
  out[st$keep] <- weighted_pi(design, st)
  out
}

#' Weights a probability-proportional-to-size design can use
#'
#' Numeric, finite and positive. Missing values are left to the design's
#' `na_rm`. A negative weight sends the rescaling in `pps_pi()` to a vector of
#' 0s and 1s summing to the wrong total. A weight of 0 is refused too, on
#' every path, and the message says to drop such rows from the frame. Every
#' message names the column.
#'
#' @noRd
check_weights <- function(w, col) {
  if (!is.numeric(w)) {
    stop("`", col, "` must be numeric, not ", class(w)[1], ".", call. = FALSE)
  }
  w <- w[!is.na(w)]
  if (any(!is.finite(w))) {
    stop("`", col, "` contains ", sum(!is.finite(w)), " non-finite weight(s).",
         call. = FALSE)
  }
  if (any(w <= 0)) {
    stop("`", col, "` must be positive, but ", sum(w <= 0), " weight(s) are ",
         "0 or negative. Drop those rows from `data` first; a design that ",
         "selects in proportion to size can never select them.",
         call. = FALSE)
  }
  invisible(w)
}
