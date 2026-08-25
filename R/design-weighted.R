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
#'     inclusion probabilities up to 35% away from proportional. Treating them as
#'     if they were proportional biases a Horvitz-Thompson total by around 1% at
#'     `n = 15` out of 30 in simulation, small in size but unmistakable in sign.
#'     There is no closed form for `pi`, so [inclusion_prob()] refuses to give
#'     one. Fine when you want a weighted selection; wrong as the basis for an
#'     estimate.}
#'   \item{`"systematic"`}{Systematic probability-proportional-to-size. Walks the
#'     cumulative weights with a fixed step from a random start, giving
#'     `pi_i = n * p_i` exactly. Rows heavy enough that `n * p_i > 1` are taken
#'     with certainty and the rest rescaled, repeatedly, until every probability
#'     is valid. Fixed sample size. Some pairs of rows can never appear together,
#'     so joint inclusion probabilities are zero for them and variance estimation
#'     needs care.}
#'   \item{`"poisson"`}{Each row is included independently with probability
#'     `pi_i = n * p_i`, capped at 1. Inclusion probabilities are exactly
#'     proportional to size and every pair can co-occur, at the cost of a
#'     **random sample size** averaging `n`.}
#' }
#'
#' `replace = TRUE` is with-replacement PPS and applies only to
#' `method = "successive"`.
#'
#' @param weights Column of positive, finite weights.
#' @param n Number of rows to draw. Under `method = "poisson"` this is the
#'   expected number, not a guarantee.
#' @param replace Sample with replacement? Only available for
#'   `method = "successive"`.
#' @param method One of `"successive"`, `"systematic"` or `"poisson"`. See
#'   "Choosing a method".
#' @param na_rm Drop rows whose weight is `NA` instead of raising an error.
#'
#' @return A design object, for use with [draw()].
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
#' @family designs
#' @seealso [draw()], [inclusion_prob()]
#' @export
design_weighted <- function(weights, n, replace = FALSE,
                            method = c("successive", "systematic", "poisson"),
                            na_rm = FALSE) {
  method <- match.arg(method)
  replace <- check_flag(replace, "replace")
  if (replace && method != "successive") {
    stop("`replace = TRUE` is with-replacement PPS and cannot be combined with ",
         "method = \"", method, "\". Use method = \"successive\", or set ",
         "replace = FALSE.", call. = FALSE)
  }
  new_design("weighted", list(
    weights = check_columns(weights, "weights", 1L),
    n = check_count(n, "n"),
    replace = replace,
    method = method,
    na_rm = check_flag(na_rm, "na_rm")
  ))
}

#' @export
draw_design.drawn_design_weighted <- function(design, data) {
  validate_data(data, required_columns = design$weights)

  col <- design$weights
  check_key_columns(data, col, "weights")
  w <- data[[col]]
  if (!is.numeric(w)) {
    stop("`", col, "` must be numeric, not ", class(w)[1], ".", call. = FALSE)
  }
  if (anyNA(w)) {
    data <- drop_na_rows(data, is.na(w), design$na_rm,
                         paste0("a missing weight in `", col, "`"))
    w <- data[[col]]
  }
  check_weights(w, col)

  n_rows <- nrow(data)

  if (design$method == "successive") {
    check_draw_size(design$n, n_rows, design$replace)
    idx <- sample.int(n_rows, size = design$n, replace = design$replace,
                      prob = w)
    return(reindex(data, idx))
  }

  if (design$n > n_rows) {
    stop("`n` (", design$n, ") cannot exceed the number of rows (", n_rows,
         ") for method = \"", design$method, "\".", call. = FALSE)
  }

  pi <- pps_pi(w, design$n)

  idx <- if (design$method == "poisson") {
    which(stats::runif(n_rows) < pi)
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
#' @noRd
pps_pi <- function(w, n) {
  pi <- rep(0, length(w))
  certain <- rep(FALSE, length(w))

  repeat {
    free <- !certain
    n_free <- n - sum(certain)
    if (n_free <= 0 || !any(free)) break
    p <- w[free] / sum(w[free])
    cand <- n_free * p
    over <- cand >= 1
    if (!any(over)) {
      pi[free] <- cand
      break
    }
    idx_free <- which(free)
    certain[idx_free[over]] <- TRUE
  }

  pi[certain] <- 1
  pi
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
  if (design$method == "successive") {
    no_closed_form(
      "`design_weighted(method = \"successive\")`",
      paste0("Successive sampling has no closed-form inclusion probability.\n",
             "Use method = \"systematic\" or \"poisson\" for a design whose ",
             "inclusion\nprobabilities really are proportional to the weights.")
    )
  }
  validate_data(data, required_columns = design$weights)
  col <- design$weights
  check_key_columns(data, col, "weights")
  w <- data[[col]]
  keep <- !check_na_policy(is.na(w), design$na_rm,
                           paste0("a missing weight in `", col, "`"))
  # pps_pi()'s rescaling loop runs off the rails on a non-positive weight --
  # `n_free` goes negative and it returns a vector of 0s and 1s summing to the
  # wrong total. draw_design() rejects these; so must this.
  check_weights(w[keep], col)
  out <- numeric(nrow(data))
  out[keep] <- pps_pi(w[keep], design$n)
  out
}

#' Weights a probability-proportional-to-size design can actually use
#' @noRd
check_weights <- function(w, col) {
  if (any(!is.finite(w))) {
    stop("`", col, "` contains ", sum(!is.finite(w)), " non-finite weight(s).",
         call. = FALSE)
  }
  if (any(w <= 0)) {
    stop("All weights must be positive; ", sum(w <= 0), " are not.",
         call. = FALSE)
  }
  invisible(w)
}
