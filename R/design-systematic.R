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
#' @param order_by Optional column to sort by before walking the data.
#' @param na_rm When `order_by` is given, drop rows whose sort key is `NA`
#'   instead of raising an error.
#'
#' @return A design object, for use with [draw()].
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
  validate_data(data, required_columns = design$order_by)

  if (!is.null(design$order_by)) {
    check_key_columns(data, design$order_by, "order_by")
    key <- data[[design$order_by]]
    if (anyNA(key)) {
      data <- drop_na_rows(data, is.na(key), design$na_rm,
                           paste0("a missing `", design$order_by, "`"))
      key <- data[[design$order_by]]
    }
    data <- data[order(key), , drop = FALSE]
  }

  start <- design$start %||% sample.int(design$interval, 1L)

  if (start > nrow(data)) {
    warning("`start` (", start, ") is past the last row (", nrow(data),
            "); returning no rows.", call. = FALSE)
    return(empty_like(data))
  }

  reindex(data, seq.int(from = start, to = nrow(data), by = design$interval))
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_systematic <- function(design, data) {
  if (!is.null(design$start)) {
    stop("A systematic design with a fixed `start` is not a probability ",
         "sample: each row is selected with probability 0 or 1. Leave `start` ",
         "as NULL for a random start.", call. = FALSE)
  }
  # Row j is selected exactly when start == ((j - 1) %% interval) + 1, which is
  # one of the `interval` equally likely starts. Exact for any N.
  rep(1 / design$interval, nrow(data))
}
