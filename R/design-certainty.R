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
#' A `rest` with a random sample size, such as
#' `design_weighted(method = "poisson")` with a small `n`, can draw none of the
#' rows below the threshold. Those rows then went unobserved, so [ht_total()]
#' and [ht_mean()] report no variance rather than an exact total.
#' `draw(weights = TRUE)` refuses a `rest` that can never draw any of them,
#' such as one with `n = 0`; without weights the draw is a cut-off sample of
#' the rows above the threshold.
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
    stop("`rest` must be a design for the rows below the threshold, not ",
         a_class(rest), ". See ?designs.", call. = FALSE)
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
  if (reserve_size(rest) > 0L) {
    stop("`rest` has `reserve = ", reserve_size(rest), "`, and a certainty ",
         "design lists no reserves. Leave `reserve` out of `rest`.",
         call. = FALSE)
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

#' The certainty split, with the frame handed to `rest`, for the variance
#'
#' `certainty_split()` plus `rest_data`, the rows below the threshold, which
#' the variance estimators and the jackknife pass to `rest`. Kept in `ctx`,
#' a `variance_context()`, so the domains of one sample split the frame once.
#'
#' @noRd
certainty_parts <- function(design, data, ctx = NULL) {
  from_context(ctx, "certainty", {
    sp <- certainty_split(design, data)
    sp$rest_data <- sp$data[sp$rest, , drop = FALSE]
    sp
  })
}

#' Describe the rows below the threshold, for messages
#' @noRd
rows_below <- function(design, n) {
  paste0(n, if (n == 1L) " row" else " rows", " with `", design$above,
         "` below ", format(design$threshold, big.mark = ",",
                            scientific = FALSE, digits = 15))
}

#' Evaluate a call on `rest`, saying so in any error it raises
#'
#' `rest` sees only the rows below the threshold, so its own messages count
#' those, as in "cannot exceed the number of rows (230)" on a frame of 240.
#' The prefix says which rows it was given. The condition keeps its class.
#'
#' @noRd
in_rest <- function(design, n_below, code) {
  tryCatch(code, error = function(e) {
    e$message <- paste0("In `rest`, applied to the ",
                        rows_below(design, n_below), ": ",
                        conditionMessage(e))
    e$call <- NULL
    stop(e)
  })
}

#' The variance when `rest` drew none of the rows below the threshold
#'
#' A `rest` with a random size, such as Poisson, can come up empty. The rows
#' below the threshold then went unobserved, and their share of the total has
#' no variance that can be measured. Calling the total exact would hide them.
#'
#' @noRd
empty_rest <- function(design, n_below) {
  list(variance = NA_real_, note = paste0(
    "No variance can be estimated for the ", rows_below(design, n_below),
    ", because `rest` drew none of them. Use a `rest` with a fixed sample ",
    "size, or a larger expected one."))
}

#' Refuse weights when `rest` can never draw a row below the threshold
#'
#' With every row below the threshold at probability 0, a weighted total
#' leaves those rows out and is biased low by their whole share. A cut-off
#' sample of the rows above the threshold is still a legitimate selection, so
#' only `weights = TRUE` is refused.
#'
#' @noRd
check_rest_reach <- function(design, data, p_all) {
  sp <- certainty_split(design, data)
  below <- sp$keep[sp$rest]
  if (length(below) && all(p_all[below] == 0)) {
    stop("`rest` gives all ", rows_below(design, length(below)),
         " an inclusion probability of 0, so a weighted total would leave ",
         "them out. Give `rest` a sample size above 0, or draw with ",
         "weights = FALSE to take the rows above the threshold alone.",
         call. = FALSE)
  }
  invisible(design)
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
  drawn_below <- in_rest(design, length(sp$rest),
                         draw_design(design$rest, below))
  reindex(data, c(sp$take, drawn_below[[key]]), sort = TRUE)
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_certainty <- function(design, data) {
  sp <- certainty_split(design, data)
  out <- numeric(nrow(data))                 # dropped rows keep their 0
  out[sp$keep[sp$take]] <- 1
  if (length(sp$rest)) {
    out[sp$keep[sp$rest]] <- in_rest(
      design, length(sp$rest),
      exact_inclusion(design$rest, sp$data[sp$rest, , drop = FALSE])
    )
  }
  out
}

#' @noRd
joint_inclusion.drawn_design_certainty <- function(design, data, rows) {
  sp <- certainty_split(design, data)
  pi_all <- exact_inclusion(design, data)
  k <- length(rows)
  out <- matrix(0, k, k)

  above <- rows %in% sp$keep[sp$take]
  # A certainty row is in every sample, so a pair containing one co-occurs
  # exactly when the other row is drawn.
  if (any(above)) {
    out[above, ] <- rep(pi_all[rows], each = sum(above))
    out[, above] <- rep(pi_all[rows], times = sum(above))
    out[above, above] <- 1
  }

  # A row the design dropped is in no sample, so every joint probability
  # involving it is 0 -- which is the answer, not an error. `rows` defaults to
  # the whole frame, so erroring here made the documented default unusable.
  free <- which(!above)
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
  out[diag_index(k)] <- pi_all[rows]
  out
}

#' Variance of a certainty design
#'
#' Rows taken with certainty are in every possible sample, so they contribute a
#' fixed amount to the total and nothing at all to its variance. What is left is
#' the variance of `rest` over the rows below the threshold -- which is where
#' this hands off, so that `rest`'s own special cases (Poisson's independent
#' units, systematic's successive-difference estimator and its refusal of a
#' fixed start) apply instead of being silently bypassed.
#'
#' The total is exact only when the frame has no rows below the threshold. A
#' `rest` that drew none of them leaves their share unobserved.
#'
#' @noRd
ht_variance.drawn_design_certainty <- function(design, data, rows, y, pi_i,
                                               domain = NULL, ctx = NULL) {
  sp <- certainty_parts(design, data, ctx)
  free <- !(rows %in% sp$keep[sp$take])
  if (!any(free)) {
    if (length(sp$rest)) return(empty_rest(design, length(sp$rest)))
    return(all_certain())
  }
  pos <- match(rows[free], sp$keep[sp$rest])
  if (anyNA(pos)) {
    return(list(variance = NA_real_, note = paste0(
      "Some sampled rows were dropped by this design's `na_rm`, so they have ",
      "no place in the frame the variance is computed over.")))
  }
  ht_variance(design$rest, sp$rest_data, pos, y[free], pi_i[free],
              rest_domain(domain, sp), sub_context(ctx, "rest"))
}

#' A domain's frame rows, restricted to the rows handed to `rest`
#'
#' `domain` is `NULL`, `NA` or a logical vector over the frame, as
#' `ht_variance_dispatch()` takes it.
#'
#' @noRd
rest_domain <- function(domain, sp) {
  if (is.null(domain) || identical(domain, NA)) return(domain)
  domain[sp$keep[sp$rest]]
}
