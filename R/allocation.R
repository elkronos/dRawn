# Splitting a total across groups. Used by the stratified and multi-stage
# designs and by their inclusion probabilities, so it lives on its own rather
# than inside whichever design happened to need it first.

#' Largest-remainder allocation
#'
#' Independent per-group `round()` neither hits the requested total (banker's
#' rounding takes four strata of 25 down to 8 rows for a requested 10) nor
#' leaves the rare-group floor optional.
#'
#' @param n Total to allocate.
#' @param sizes Group sizes, named.
#' @param allocation `"proportional"` or `"equal"`.
#' @param min_per_stratum Floor per group.
#' @param cap When `TRUE`, no group may be allocated more rows than it holds.
#'   `FALSE` when sampling with replacement, where over-allocation is legal.
#' @param spread Per-group standard deviation, for Neyman allocation. Groups
#'   get shares proportional to `size * spread`, which minimises the variance
#'   of a total for a fixed `n`: a group whose values barely vary needs fewer
#'   rows to pin down than an equally sized group that varies a lot.
#' @noRd
allocate <- function(n, sizes, allocation, min_per_stratum, cap = TRUE,
                     spread = NULL) {
  k <- length(sizes)
  if (k == 0L) return(integer(0))

  if (min_per_stratum * k > n) {
    stop("min_per_stratum = ", min_per_stratum, " across ", k,
         " group(s) needs at least ", min_per_stratum * k,
         " rows, but `n` is ", n, ".", call. = FALSE)
  }

  raw <- switch(allocation,
    equal = rep(n / k, k),
    neyman = {
      if (is.null(spread)) {
        stop("Neyman allocation needs `allocation_by`.", call. = FALSE)
      }
      sp <- spread[names(sizes)]
      sp[!is.finite(sp)] <- 0
      w <- sizes * sp
      # Every stratum with no spread still needs representation; fall back to
      # proportional if the auxiliary variable is constant throughout.
      if (sum(w) <= 0) n * sizes / sum(sizes) else n * w / sum(w)
    },
    n * sizes / sum(sizes)
  )

  base <- pmax(min_per_stratum, floor(raw))
  if (cap) base <- pmin(base, sizes)
  short <- n - sum(base)

  if (short > 0) {
    room <- if (cap) sizes - base else rep(Inf, k)
    ord <- order(raw - floor(raw), room, decreasing = TRUE)
    ord <- ord[room[ord] > 0]
    # Deal the remainder out one row at a time so no group overflows.
    i <- 1L
    while (short > 0 && length(ord) > 0) {
      j <- ord[(i - 1L) %% length(ord) + 1L]
      if (!cap || base[j] < sizes[j]) {
        base[j] <- base[j] + 1L
        short <- short - 1L
      } else {
        ord <- ord[ord != j]
        i <- i - 1L
      }
      i <- i + 1L
    }
  } else if (short < 0) {
    # Mirror the branch above: take rows back one at a time, cycling, so a
    # group can give up more than one. A single pass could only remove as many
    # rows as there are groups, which left the total over target whenever
    # min_per_stratum pushed the floor up by more than one row per group.
    ord <- order(raw - floor(raw), decreasing = FALSE)
    i <- 1L
    while (short < 0 && length(ord) > 0) {
      j <- ord[(i - 1L) %% length(ord) + 1L]
      if (base[j] > min_per_stratum) {
        base[j] <- base[j] - 1L
        short <- short + 1L
      } else {
        ord <- ord[ord != j]
        i <- i - 1L
      }
      i <- i + 1L
    }
  }

  stats::setNames(as.integer(base), names(sizes))
}
