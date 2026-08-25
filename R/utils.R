# General-purpose helpers with no design-specific knowledge.
# Argument checking lives in validate.R, allocation in allocation.R,
# cluster selection in clusters.R.

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Evaluate an expression under a seed without disturbing the caller's RNG
#'
#' Saves `.Random.seed`, seeds, evaluates, then restores the previous state on
#' exit. Calling [set.seed()] directly would permanently move the caller's
#' random number stream.
#'
#' @noRd
with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != trunc(seed)) {
    stop("`seed` must be a single whole number, not ", format_bad(seed), ".",
         call. = FALSE)
  }
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())),
            add = TRUE)
  }
  set.seed(seed)
  force(code)
}

#' Subset rows and reset row names
#'
#' `drop = FALSE` matters: without it a one-column frame collapses to a vector.
#'
#' @noRd
reindex <- function(data, idx, sort = FALSE) {
  if (sort) idx <- sort(idx)
  out <- data[idx, , drop = FALSE]
  rownames(out) <- NULL
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

#' Sample `size` values from the elements of `x`
#'
#' `sample(x, size)` reinterprets a length-1 numeric `x` as `seq_len(x)`, so a
#' frame with a single numeric cluster id would sample integers that were never
#' in the data.
#'
#' @noRd
sample_values <- function(x, size, replace = FALSE) {
  x[sample.int(length(x), size = size, replace = replace)]
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
         "them.\nIf you want the expected number of times each row appears, ",
         "that is n_replicates * n / nrow(data).", call. = FALSE)
  }
  invisible(design)
}
