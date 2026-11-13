#' Systematic sampling
#'
#' Takes every `interval`-th row, beginning at `start`. The number of rows drawn
#' follows from the interval and the size of the data, so there is no `n`.
#'
#' @param interval Sampling interval. A whole number of 1 or more. Fractional
#'   intervals are rejected: they produce uneven gaps, which is not a systematic
#'   design.
#' @param start Starting row, in `1:interval`. Drawn at random from that range
#'   when `NULL`.
#' @param order_by Optional column to sort by before walking the data. It
#'   changes which rows can appear together, so [inclusion_prob()] and
#'   [joint_prob()] compute against the sorted order too. Text is sorted byte
#'   by byte, as in the C locale, so the walk is the same in every session;
#'   ties keep frame order.
#' @param na_rm When `order_by` is given, drop rows whose sort key is `NA`
#'   instead of raising an error.
#'
#' @section Variance:
#' A systematic sample has a single random start, so no design-unbiased
#' variance estimator exists. [ht_total()] uses the successive-difference
#' approximation (Wolter 2007), which compares neighbouring sampled rows in
#' the order the design walked. Sorting on a variable related to what you
#' measure (`order_by`) is what makes systematic sampling efficient, and this
#' estimator is the one that can see it; see [ht_total()] for its limits.
#' When `interval` does not divide the number of rows walked, the sample size
#' is random, and a term for that is added to the variance.
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
#' @family designs
#' @seealso [draw()]
#' @export
design_systematic <- function(interval, start = NULL, order_by = NULL,
                              na_rm = FALSE) {
  interval <- check_count(interval, "interval", allow_zero = FALSE)
  if (!is.null(start)) {
    start <- check_count(start, "start", allow_zero = FALSE)
    if (start > interval) {
      stop("`start` (", start, ") must lie in 1:interval (1:", interval, ").",
           call. = FALSE)
    }
  }
  new_design("systematic", list(
    interval = interval,
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

  start <- design$start %||% sample.int(design$interval, 1L)

  if (start > nrow(data)) {
    warning("`start` (", start, ") is past the last row (", nrow(data),
            "); returning no rows.", call. = FALSE)
    return(empty_like(data))
  }

  picked <- walk[seq.int(from = start, to = nrow(data), by = design$interval)]
  reindex(data, picked, sort = TRUE)
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
#' policy, and a key missing in every row.
#'
#' @noRd
systematic_positions <- function(design, data) {
  n <- nrow(data)
  if (is.null(design$order_by)) return(seq_len(n))
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
  out
}

#' @noRd
exact_inclusion.drawn_design_systematic <- function(design, data) {
  if (!is.null(design$start)) {
    stop_declined(paste0(
      "A systematic design with a fixed `start` is not a probability sample: ",
      "each row is selected with probability 0 or 1. Leave `start` as NULL ",
      "for a random start."))
  }
  # Row j is selected exactly when start == ((j - 1) %% interval) + 1, which is
  # one of the `interval` equally likely starts. Exact for any N.
  pos <- systematic_positions(design, data)
  ifelse(is.na(pos), 0, 1 / design$interval)
}
