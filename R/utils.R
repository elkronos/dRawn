# General-purpose helpers with no design-specific knowledge.
# Argument checking lives in validate.R, allocation in allocation.R,
# cluster selection in clusters.R.

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Is an inclusion probability 1?
#'
#' The one test for certainty. `n * w / sum(w)` for a unit that is certain in
#' exact arithmetic can come out a rounding error below 1, so a probability
#' within `1e-9` of 1 counts as 1.
#'
#' @noRd
is_certain <- function(p) {
  p >= 1 - 1e-9
}

#' Signal a warning with a class of its own
#'
#' The class is `c(class, "drawn_warning", "warning", "condition")`, so a
#' caller can muffle or catch one kind of warning with
#' `withCallingHandlers()` or `tryCatch()` without matching on its text.
#' `class` is named `drawn_<what>`. The warning carries no call, as the
#' package's errors do not.
#'
#' @noRd
warn_drawn <- function(message, class) {
  cond <- structure(list(message = message, call = NULL),
                    class = c(class, "drawn_warning", "warning", "condition"))
  warning(cond)
}

#' Refuse a probability or variance that the design does not define
#'
#' The condition has class `c("drawn_declined", "error", "condition")` and no
#' call. Called directly, as from [joint_prob()], it is an ordinary error.
#' `ht_variance_dispatch()` catches this class and no other, so the estimator
#' declining falls back to the jackknife or reports no variance, while any
#' other error, from running out of memory to a bug, reaches the caller.
#'
#' @noRd
stop_declined <- function(message) {
  cond <- structure(list(message = message, call = NULL),
                    class = c("drawn_declined", "error", "condition"))
  stop(cond)
}

#' Evaluate an expression under a seed without disturbing the caller's RNG
#'
#' Saves `.Random.seed`, seeds, evaluates, then restores the previous state on
#' exit. Calling [set.seed()] directly would permanently move the caller's
#' random number stream.
#'
#' The seed is applied with R's default generator, whatever the session has
#' selected, so a seed draws the same rows after `RNGkind("L'Ecuyer-CMRG")`
#' or under `sample.kind = "Rounding"`. The first element of `.Random.seed`
#' records the generator, so restoring it restores the caller's generator
#' too. With no `.Random.seed` to restore, the generator is set back with
#' RNGkind(), which seeds it afresh, and the `.Random.seed` that leaves is
#' removed, as the caller had none. RNGkind() warns when it sets the
#' "Rounding" sampler, and that warning is muffled, since the caller chose
#' that sampler and was warned then.
#'
#' @noRd
with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }
  check_seed(seed)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    old_kind <- RNGkind()
    on.exit({
      suppressWarnings(do.call(RNGkind, as.list(old_kind)))
      suppressWarnings(rm(".Random.seed", envir = globalenv()))
    }, add = TRUE)
  }
  kind <- seed_kind()
  set.seed(seed, kind = kind[[1L]], normal.kind = kind[[2L]],
           sample.kind = kind[[3L]])
  force(code)
}

#' The generator a seed is applied with
#'
#' R's defaults since 3.6.0, so in a session left at the defaults a seed
#' draws what `set.seed(seed)` alone would.
#'
#' @noRd
seed_kind <- function() c("Mersenne-Twister", "Inversion", "Rejection")

#' A seed set.seed() takes: one whole number in the integer range
#'
#' Past `.Machine$integer.max` set.seed() gives its own error after a
#' coercion warning, so the range is checked here, with the value.
#'
#' @noRd
check_seed <- function(seed) {
  if (is.null(seed)) return(invisible(NULL))
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != trunc(seed)) {
    stop("`seed` must be a single whole number, not ",
         format_bad(seed, value = TRUE), ".", call. = FALSE)
  }
  if (abs(seed) > .Machine$integer.max) {
    stop("`seed` must lie between -", .Machine$integer.max, " and ",
         .Machine$integer.max, ", not ", format(seed, scientific = FALSE),
         ". Use a smaller whole number.", call. = FALSE)
  }
  invisible(seed)
}

#' Seconds since an arbitrary origin, for timing simulated draws
#' @noRd
elapsed_time <- function() proc.time()[["elapsed"]]

#' Is R running interactively?
#'
#' A function of the package's own, so a test can stand in for it.
#'
#' @noRd
is_interactive <- function() interactive()

#' A duration in seconds, in words
#' @noRd
format_duration <- function(seconds) {
  unit <- if (seconds < 90) {
    c(1, "second")
  } else if (seconds < 90 * 60) {
    c(60, "minute")
  } else {
    c(3600, "hour")
  }
  x <- signif(seconds / as.numeric(unit[[1L]]), 2)
  paste0(format(x, big.mark = ",", scientific = FALSE), " ", unit[[2L]],
         if (x != 1) "s")
}

#' Subset rows of a frame
#'
#' `drop = FALSE` matters: without it a one-column frame collapses to a vector.
#'
#' Row names of the frame's own, such as the car names of `mtcars`, stay with
#' their rows. A row taken more than once gets names made unique by
#' make.unique(), "Mazda RX4" then "Mazda RX4.1", as `[` gives them. The
#' automatic row names 1, 2, ... are numbered afresh, since they name
#' positions rather than rows. `.row_names_info()` is negative for automatic
#' row names and positive for names of the frame's own. A tibble's `[` drops
#' row names, and that is kept.
#'
#' A plain data.frame is subset column by column. `[.data.frame` gives the
#' same result but also builds unique row names for the repeated rows, which
#' are then dropped. For the 5 million rows of a default bootstrap of a
#' 5,000-row frame with four columns, that took 8 seconds, and this takes 0.1.
#' Giving those rows unique names made from row names of the frame's own
#' takes about 8 seconds too, which only such a frame pays.
#'
#' @noRd
reindex <- function(data, idx, sort = FALSE) {
  if (sort) idx <- sort(idx)
  own <- .row_names_info(data) > 0L
  if (identical(class(data), "data.frame")) {
    out <- lapply(unclass(data), function(col) {
      if (length(dim(col)) == 2L) col[idx, , drop = FALSE] else col[idx]
    })
    attributes(out) <- attributes(data)
    attr(out, "row.names") <- if (own) {
      kept_names(attr(data, "row.names"), idx)
    } else {
      .set_row_names(length(idx))
    }
    return(out)
  }
  out <- data[idx, , drop = FALSE]
  if (!own) rownames(out) <- NULL
  out
}

#' Row names of the frame's own at `idx`, made unique when a row repeats
#' @noRd
kept_names <- function(names, idx) {
  out <- names[idx]
  if (anyDuplicated(idx)) out <- make.unique(as.character(out))
  out
}

#' Draw `n[[i]]` positions from each group of indices
#' @noRd
take_within <- function(idx_by_group, n, replace) {
  unlist(
    Map(function(idx, k) {
      if (k == 0L) return(integer(0))
      idx[sample.int(length(idx), size = k, replace = replace)]
    }, idx_by_group, n),
    use.names = FALSE
  )
}

#' An empty result that still has the input's schema
#'
#' `do.call(rbind, list())` is `NULL`, which would turn an empty window into a
#' 0x0 frame with no columns at all.
#'
#' @noRd
empty_like <- function(data) {
  data[0L, , drop = FALSE]
}

#' Require a suggested package at run time
#' @noRd
require_suggested <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package `", pkg, "` is required for ", what, ".\n",
         "Install it with install.packages(\"", pkg, "\").", call. = FALSE)
  }
  invisible(TRUE)
}

#' Designs that simulation cannot rescue
#'
#' Counting how often a row appears across repeated draws estimates an inclusion
#' probability only if the design is a probability sample of the frame. A
#' bootstrap resamples the sample, so every row appears in some replicate and
#' the count converges to 1 for every row -- a confident-looking answer that
#' means nothing. Refuse rather than return it.
#'
#' @noRd
refuse_simulation <- function(design) {
  if (inherits(design, "drawn_design_bootstrap")) {
    stop("A bootstrap is not a probability sample of a finite population, so ",
         "simulating it does not estimate an inclusion probability: every row ",
         "appears in some replicate, and the count converges to 1 for all of ",
         "them.\n", bootstrap_counts(design), call. = FALSE)
  }
  invisible(design)
}
