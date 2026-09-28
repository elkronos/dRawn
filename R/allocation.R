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
#' @param allocation `"proportional"`, `"equal"` or `"neyman"`.
#' @param min_per_stratum Floor per group.
#' @param cap When `TRUE`, no group may be allocated more rows than it holds.
#'   `FALSE` when sampling with replacement, where over-allocation is legal.
#' @param spread Per-group standard deviation, for Neyman allocation. Groups
#'   get shares proportional to `size * spread`, which minimises the variance
#'   of a total for a fixed `n`: a group whose values barely vary needs fewer
#'   rows to pin down than an equally sized group that varies a lot.
#' @noRd
allocate <- function(n, sizes, allocation, min_per_stratum, cap = TRUE,
                     spread = NULL, min_arg = "min_per_stratum") {
  k <- length(sizes)
  if (k == 0L) return(integer(0))

  if (min_per_stratum * k > n) {
    stop(min_arg, " = ", min_per_stratum, " across ", k,
         " group(s) needs at least ", min_per_stratum * k,
         " rows, but `n` is ", n, ".", call. = FALSE)
  }
  # `n` and `sizes` are both integer, so `n * sizes` is evaluated in 32 bits and
  # overflows to NA once the product passes 2^31 -- 3,000 rows out of 2.4
  # million was enough. Everything downstream is fractional anyway.
  n <- as.numeric(n)
  sizes <- stats::setNames(as.numeric(sizes), names(sizes))

  share <- switch(allocation,
    equal = rep(1, k),
    neyman = {
      if (is.null(spread)) {
        stop("Neyman allocation needs `allocation_by`.", call. = FALSE)
      }
      sp <- spread[names(sizes)]
      sp[!is.finite(sp)] <- 0
      # A constant auxiliary variable has no spread to allocate by; fall back
      # to proportional rather than dividing by zero.
      if (sum(sizes * sp) <= 0) sizes else sizes * sp
    },
    sizes
  )
  raw <- capped_shares(n, share, sizes, cap)

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

#' Continuous allocation with take-all capping
#'
#' Splitting `n` in proportion to `share` can ask a group for more rows than it
#' holds -- routinely under Neyman allocation, where a small, highly variable
#' stratum attracts a large share, and under equal allocation. The textbook
#' remedy (Cochran 1977, section 5.9) takes such a group whole and re-splits
#' what is left over the others *in the same proportions*. Dealing the surplus
#' out one row at a time instead, as rounding does, spreads it evenly and
#' quietly abandons the allocation rule that was asked for.
#'
#' @noRd
capped_shares <- function(n, share, sizes, cap) {
  k <- length(sizes)
  raw <- numeric(k)
  free <- rep(TRUE, k)
  left <- n
  repeat {
    w <- share[free]
    raw[free] <- if (sum(w) > 0) left * w / sum(w) else left / sum(free)
    if (!cap) break
    over <- free & raw > sizes
    if (!any(over)) break
    raw[over] <- sizes[over]
    free[over] <- FALSE
    left <- n - sum(raw[!free])
    if (!any(free) || left <= 0) {
      raw[free] <- 0
      break
    }
  }
  raw
}

#' Warn when a group with rows in the frame is allocated none of them
#'
#' A stratum allocated zero rows makes every one of its frame rows unreachable:
#' inclusion probability 0, so a Horvitz-Thompson total omits them and is
#' biased by exactly their total. Proportional allocation does this silently to
#' any stratum smaller than `N / n`, and Neyman allocation to any stratum whose
#' auxiliary variable is constant. It is almost never what was meant.
#'
#' The warning carries its own class so that loops which draw on purpose many
#' times -- the Monte Carlo probabilities -- can silence the repeats.
#'
#' @noRd
warn_empty_groups <- function(n_alloc, sizes, noun = "stratum",
                              arg = "min_per_stratum") {
  empty <- names(n_alloc)[n_alloc == 0L & sizes > 0]
  if (!length(empty)) return(invisible(NULL))
  shown <- group_label(utils::head(empty, 5L))
  more <- if (length(empty) > 5L) paste0(" and ", length(empty) - 5L, " more")
          else ""
  plural <- length(empty) > 1L
  msg <- paste0(
    length(empty), " ", if (plural) paste0(sub("um$", "a", noun)) else noun,
    " (", paste0("`", shown, "`", collapse = ", "), more, ") ",
    if (plural) "are" else "is", " allocated no rows, so ",
    if (plural) "their " else "its ", fmt_n(sum(sizes[empty])),
    " frame row(s) can never be drawn. A total estimated from this sample ",
    "leaves them out entirely. Raise `n`, or set ", arg, " = 1 (or 2, so ",
    "each ", noun, " also supports a variance)."
  )
  warning(structure(class = c("drawn_empty_stratum", "warning", "condition"),
                    list(message = msg, call = NULL)))
}

#' Evaluate `code`, silencing the empty-stratum warning
#'
#' For loops that draw the same design many times on purpose; the result they
#' return already shows the unreachable rows as zeros.
#'
#' @noRd
quietly_empty <- function(code) {
  withCallingHandlers(code, drawn_empty_stratum = function(w) {
    invokeRestart("muffleWarning")
  })
}
