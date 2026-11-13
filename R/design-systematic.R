#' Systematic sampling
#'
#' Walks the rows at a fixed step. Give `interval` to take every
#' `interval`-th row, beginning at `start`, which is drawn at random unless
#' given; the number of rows drawn then follows from the interval and the size
#' of the data. Give `n` instead to draw exactly `n` rows from a random start:
#' the step is the number of rows walked divided by `n`, and need not be a
#' whole number.
#'
#' @param interval Sampling interval, a whole number of 1 or more: every
#'   `interval`-th row of the walk is taken. Give `interval` or `n`, not both.
#'   For a step that is not a whole number, give `n`.
#' @param start Starting row, in `1:interval`. Drawn at random from that range
#'   when `NULL`. A fixed start makes the design deterministic: each row is
#'   selected with probability 0 or 1, so [inclusion_prob()],
#'   [sampling_weight()], [joint_prob()], `draw(weights = TRUE)`,
#'   [ht_total()] and [ht_mean()] refuse it, while
#'   `inclusion_prob(simulate = TRUE)` reports those 0s and 1s. Not used with
#'   `n`, whose start is always random.
#' @param order_by Optional column to sort by before walking the data. It
#'   changes which rows can appear together, so [inclusion_prob()] and
#'   [joint_prob()] compute against the sorted order too. Text is sorted byte
#'   by byte, as in the C locale, so the walk is the same in every session;
#'   ties keep frame order.
#' @param na_rm When `order_by` is given, drop rows whose sort key is `NA`
#'   instead of raising an error.
#' @param n Number of rows to draw, in place of `interval`. With `N` rows
#'   walked, each row covers a stretch of length `n / N` on a line from 0 to
#'   `n`, a start `u` is drawn uniformly from (0, 1), and the rows under `u`,
#'   `u + 1`, ..., `u + n - 1` are taken. The step is `N / n` rows, every row
#'   has probability `n / N`, and the sample always holds `n` rows. When `n`
#'   divides `N` the probabilities are those of `interval = N / n`.
#'
#' @section Variance:
#' A systematic sample has a single random start, so no design-unbiased
#' variance estimator exists. [ht_total()] uses the successive-difference
#' approximation (Wolter 2007), which compares neighbouring sampled rows in
#' the order the design walked. Sorting on a variable related to what you
#' measure (`order_by`) is what makes systematic sampling efficient, and this
#' estimator is the one that can see it; see [ht_total()] for its limits.
#' When `interval` does not divide the number of rows walked, the sample size
#' is random, and a term for that is added to the variance. A design given
#' `n` has a fixed size and needs no such term.
#'
#' @return A design object, for use with [draw()].
#'
#' @references
#' Madow, W. G. and Madow, L. H. (1944). On the theory of systematic sampling,
#' I. *Annals of Mathematical Statistics*, 15, 1–24.
#'
#' Wolter, K. M. (2007). *Introduction to Variance Estimation*, 2nd ed.
#' Springer.
#'
#' @examples
#' df <- data.frame(id = 1:100, value = (1:100) / 10)
#' draw(df, design_systematic(interval = 10, start = 3))
#' draw(df, design_systematic(interval = 25, order_by = "value"), seed = 1)
#'
#' # 30 rows from 100, a step of 3.33
#' draw(df, design_systematic(n = 30), seed = 1)$id
#'
#' @family designs
#' @seealso [draw()]
#' @export
design_systematic <- function(interval = NULL, start = NULL, order_by = NULL,
                              na_rm = FALSE, n = NULL) {
  if (is.null(interval) == is.null(n)) {
    stop(if (is.null(n)) {
      "`design_systematic()` needs `interval` or `n`, and has neither. "
    } else {
      "`design_systematic()` takes `interval` or `n`, not both. "
    }, "Give `interval` for every k-th row, or `n` for a fixed sample size.",
    call. = FALSE)
  }
  if (!is.null(n)) {
    n <- check_count(n, "n")
    if (!is.null(start)) {
      stop("`start` cannot be fixed with `n`, whose start is drawn from ",
           "(0, 1). Leave `start` as NULL, or give `interval` instead of `n`.",
           call. = FALSE)
    }
  } else {
    if (is.numeric(interval) && length(interval) == 1L &&
        is.finite(interval) && interval > 1 && interval != trunc(interval)) {
      stop("`interval` must be a single whole number, not ", interval, ". ",
           "For a step that is not whole, give the sample size as `n`.",
           call. = FALSE)
    }
    interval <- check_count(interval, "interval", allow_zero = FALSE)
  }
  if (!is.null(start)) {
    start <- check_count(start, "start", allow_zero = FALSE)
    if (start > interval) {
      stop("`start` (", start, ") must lie in 1:interval (1:", interval, ").",
           call. = FALSE)
    }
  }
  new_design("systematic", list(
    interval = interval,
    n = n,
    start = start,
    order_by = if (is.null(order_by)) NULL else check_columns(order_by, "order_by", 1L),
    na_rm = check_flag(na_rm, "na_rm")
  ))
}

#' @export
draw_design.drawn_design_systematic <- function(design, data) {
  validate_data(data)

  # The walk visits rows in `walk` order; the result comes back in frame order
  # like every other design that selects a set of rows.
  pos <- systematic_positions(design, data)
  if (anyNA(pos)) {
    keep <- which(!is.na(pos))
    data <- data[keep, , drop = FALSE]
    pos <- pos[keep]
  }
  walk <- order(pos)

  if (!is.null(design$n)) {
    if (design$n == 0L) return(empty_like(data))
    taken <- fractional_walk(stats::runif(1L), design$n, nrow(data))
    return(reindex(data, walk[taken], sort = TRUE))
  }

  start <- design$start %||% sample.int(design$interval, 1L)

  if (start > nrow(data)) {
    warning("`start` (", start, ") is past the last row (", nrow(data),
            "); returning no rows.", call. = FALSE)
    return(empty_like(data))
  }

  picked <- walk[seq.int(from = start, to = nrow(data), by = design$interval)]
  reindex(data, picked, sort = TRUE)
}

#' Walk positions a fixed-size systematic sample takes
#'
#' Walk position `j` covers `[(j - 1) * n / N, j * n / N)` on a line from 0 to
#' `n`, and the positions under `u`, `u + 1`, ..., `u + n - 1` are taken. No
#' stretch is longer than 1, so no position is hit twice and exactly `n` are
#' taken. Each break is `j * n / N` with one rounding, rather than a running
#' sum, so the last one is `n` exactly.
#'
#' @noRd
fractional_walk <- function(u, n, N) {
  breaks <- as.double(seq_len(N)) * n / N
  # A point that rounds up to `n` lies in the last stretch.
  pmin(findInterval(u + seq.int(0L, n - 1L), breaks) + 1L, N)
}

# ---- inclusion probability ------------------------------------------------

#' Where each row of `data` sits in the order the design walks
#'
#' `NA` for a row the design drops. Systematic sampling is the one design whose
#' answer depends on *position*, so both the first- and second-order
#' probabilities have to be computed against the order `draw_design()` actually
#' uses -- the sorted, NA-dropped frame -- not against frame order.
#'
#' It is also the design's one set of checks, shared by the draw, the
#' probabilities and the successive-difference variance: an `order_by` column
#' the frame lacks, which gave every row a probability of 0, the `na_rm`
#' policy, a key missing in every row, and an `n` larger than the rows walked.
#'
#' @noRd
systematic_positions <- function(design, data) {
  n <- nrow(data)
  out <- seq_len(n)
  if (!is.null(design$order_by)) {
    col <- design$order_by
    if (!col %in% names(data)) {
      stop("`data` is missing column: `", col, "`.", call. = FALSE)
    }
    check_key_columns(data, col, "order_by")
    key <- data[[col]]
    what <- paste0("a missing `", col, "`")
    bad <- check_na_policy(is.na(key), design$na_rm, what)
    if (n > 0L && all(bad)) stop("Every row has ", what, ".", call. = FALSE)
    out <- rep(NA_integer_, n)
    kept <- which(!bad)
    # Byte order for text, so the walk is the same in every locale. Ties keep
    # frame order.
    out[kept[key_order(key[kept])]] <- seq_along(kept)
  }
  if (!is.null(design$n)) {
    walked <- sum(!is.na(out))
    if (design$n > walked) {
      stop("`n` (", design$n, ") cannot exceed the number of rows walked (",
           walked, "). Give `n` of at most ", walked, ".", call. = FALSE)
    }
  }
  out
}

#' Refuse a fixed start on every probability and variance path
#'
#' A fixed start makes the design deterministic. Its probabilities of 0 and 1
#' are true, and the simulation reports them, but as weights they would make
#' the sum over the sample a total of the frame.
#'
#' @noRd
refuse_fixed_start <- function(design) {
  if (!is.null(design$start)) {
    stop_declined(paste0(
      "A systematic design with a fixed `start` is not a probability sample: ",
      "each row is selected with probability 0 or 1. Leave `start` as NULL ",
      "for a random start."))
  }
  invisible(design)
}

#' @noRd
exact_inclusion.drawn_design_systematic <- function(design, data) {
  refuse_fixed_start(design)
  pos <- systematic_positions(design, data)
  if (!is.null(design$n)) {
    # Each walked row covers a stretch of length n / N, and the start is
    # uniform on (0, 1).
    return(ifelse(is.na(pos), 0, design$n / sum(!is.na(pos))))
  }
  # Row j is selected exactly when start == ((j - 1) %% interval) + 1, which is
  # one of the `interval` equally likely starts. Exact for any N.
  ifelse(is.na(pos), 0, 1 / design$interval)
}
