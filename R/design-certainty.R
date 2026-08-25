#' Take some rows with certainty, sample the rest
#'
#' Composes a design out of two parts: rows at or above a threshold are all
#' taken, and any other design is applied to what remains. This is the standard
#' shape of an audit or financial sample — every invoice over a hundred
#' thousand gets examined, the long tail below is sampled — and it is the
#' cheapest way to cut variance when a few units dominate the total.
#'
#' Certainty rows have an inclusion probability of exactly 1, so they contribute
#' their own value to a Horvitz-Thompson total with a weight of 1 and add
#' nothing to its variance. That is the whole point: the uncertainty in the
#' estimate comes only from the part you sampled.
#'
#' @param above Column holding the size measure.
#' @param threshold Rows with `above >= threshold` are taken with certainty.
#' @param rest A design applied to the rows below the threshold. Its `n` is the
#'   number drawn from *those* rows, not the total.
#' @param na_rm Drop rows whose size measure is `NA` instead of raising an
#'   error.
#'
#' @return A design object, for use with [draw()].
#'
#' @examples
#' set.seed(1)
#' invoices <- data.frame(
#'   id = 1:500,
#'   site = rep(c("a", "b"), times = c(300, 200)),
#'   value = round(stats::rlnorm(500, 7, 1.4))
#' )
#'
#' # Everything over 20,000 examined; 40 drawn from the rest
#' d <- design_certainty("value", threshold = 20000,
#'                       rest = design_stratified("site", n = 40))
#' s <- draw(invoices, d, seed = 1, weights = TRUE)
#'
#' table(certain = s$.prob == 1)
#' sum(invoices$value >= 20000)
#'
#' @family designs
#' @seealso [draw()], [design_weighted()] for probability-proportional-to-size,
#'   which handles dominant units by taking them with certainty automatically.
#' @export
design_certainty <- function(above, threshold, rest, na_rm = FALSE) {
  if (!is_design(rest)) {
    stop("`rest` must be a design for the rows below the threshold, not a ",
         class(rest)[1], ". See ?designs.", call. = FALSE)
  }
  if (inherits(rest, "drawn_design_certainty")) {
    stop("`rest` cannot itself be a certainty design; use one threshold.",
         call. = FALSE)
  }
  if (inherits(rest, "drawn_design_bootstrap")) {
    stop("`rest` cannot be a bootstrap design: it returns replicates rather ",
         "than a sample of rows, and composing the two would drop the ",
         "`.replicate` column and present resampled duplicates as distinct ",
         "units.", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
    stop("`threshold` must be a single number.", call. = FALSE)
  }
  new_design("certainty", list(
    above = check_columns(above, "above", 1L),
    threshold = threshold,
    rest = rest,
    na_rm = check_flag(na_rm, "na_rm")
  ))
}

#' Split a frame into the certainty part and the rest
#'
#' `take` and `rest` index `data`, the frame with any dropped rows removed;
#' `keep` maps those positions back to rows of the frame that was passed in.
#' Every caller needs both, because `draw_design()` works in the reduced frame
#' while the probability methods must answer in the original one.
#'
#' @noRd
certainty_split <- function(design, data) {
  validate_data(data, required_columns = design$above)
  check_key_columns(data, design$above, "above")
  v <- data[[design$above]]
  if (!is.numeric(v)) {
    stop("`above` names `", design$above, "`, which must be numeric, not ",
         class(v)[1], ".", call. = FALSE)
  }
  keep <- seq_len(nrow(data))
  if (anyNA(v)) {
    bad <- is.na(v)
    data <- drop_na_rows(data, bad, design$na_rm,
                         paste0("a missing `", design$above, "`"))
    keep <- keep[!bad]
    v <- data[[design$above]]
  }
  take <- which(v >= design$threshold)
  list(data = data, keep = keep, take = take,
       rest = setdiff(seq_len(nrow(data)), take))
}

#' @export
draw_design.drawn_design_certainty <- function(design, data) {
  sp <- certainty_split(design, data)
  data <- sp$data

  if (length(sp$rest) == 0L) {
    return(reindex(data, sp$take, sort = TRUE))
  }
  below <- data[sp$rest, , drop = FALSE]
  key <- ".drawn_certainty_id"
  if (key %in% names(below)) {
    stop("`data` already has a column called `", key, "`. Rename it.",
         call. = FALSE)
  }
  below[[key]] <- sp$rest
  drawn_below <- draw_design(design$rest, below)
  reindex(data, c(sp$take, drawn_below[[key]]), sort = TRUE)
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_certainty <- function(design, data) {
  sp <- certainty_split(design, data)
  out <- numeric(nrow(data))                 # dropped rows keep their 0
  out[sp$keep[sp$take]] <- 1
  if (length(sp$rest)) {
    out[sp$keep[sp$rest]] <- exact_inclusion(design$rest,
                                             sp$data[sp$rest, , drop = FALSE])
  }
  out
}

#' @noRd
joint_inclusion.drawn_design_certainty <- function(design, data, rows) {
  sp <- certainty_split(design, data)
  pi_all <- exact_inclusion(design, data)
  k <- length(rows)
  out <- matrix(0, k, k)

  is_certain <- rows %in% sp$keep[sp$take]
  # A certainty row is in every sample, so a pair containing one co-occurs
  # exactly when the other row is drawn.
  if (any(is_certain)) {
    out[is_certain, ] <- rep(pi_all[rows], each = sum(is_certain))
    out[, is_certain] <- rep(pi_all[rows], times = sum(is_certain))
    out[is_certain, is_certain] <- 1
  }

  # A row the design dropped is in no sample, so every joint probability
  # involving it is 0 -- which is the answer, not an error. `rows` defaults to
  # the whole frame, so erroring here made the documented default unusable.
  free <- which(!is_certain)
  pos <- match(rows[free], sp$keep[sp$rest])
  live <- free[!is.na(pos)]
  if (length(live)) {
    out[live, live] <- joint_inclusion(design$rest,
                                       sp$data[sp$rest, , drop = FALSE],
                                       pos[!is.na(pos)])
  }
  gone <- free[is.na(pos)]
  if (length(gone)) {
    out[gone, ] <- 0
    out[, gone] <- 0
  }
  diag(out) <- pi_all[rows]
  out
}

#' Variance of a certainty design
#'
#' Rows taken with certainty are in every possible sample, so they contribute a
#' fixed amount to the total and nothing at all to its variance. What is left is
#' the variance of `rest` over the rows below the threshold -- which is where
#' this hands off, so that `rest`'s own special cases (Poisson's independent
#' units, systematic's refusal) apply instead of being silently bypassed.
#'
#' @noRd
ht_variance.drawn_design_certainty <- function(design, data, rows, y, pi_i) {
  sp <- certainty_split(design, data)
  free <- !(rows %in% sp$keep[sp$take])
  if (!any(free)) {
    return(list(variance = 0,
                note = paste0("Every sampled row was taken with certainty, so ",
                              "the total is exact rather than estimated.")))
  }
  pos <- match(rows[free], sp$keep[sp$rest])
  if (anyNA(pos)) {
    return(list(variance = NA_real_, note = paste0(
      "Some sampled rows were dropped by this design's `na_rm`, so they have ",
      "no place in the frame the variance is computed over.")))
  }
  ht_variance(design$rest, sp$data[sp$rest, , drop = FALSE], pos,
              y[free], pi_i[free])
}
