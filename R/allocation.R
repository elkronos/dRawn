# Splitting a total across groups. Used by the stratified and multi-stage
# designs and by their inclusion probabilities, so it lives on its own rather
# than inside whichever design happened to need it first.

#' Largest-remainder allocation
#'
#' Independent per-group `round()` neither hits the requested total (banker's
#' rounding takes four strata of 25 down to 8 rows for a requested 10) nor
#' leaves the rare-group floor optional.
#'
#' The split itself is continuous and comes from `capped_shares()`, which
#' already respects the floor and the cap. Only the rounding happens here:
#' each group gets the whole part of its share, and the rows left over go one
#' at a time to the largest fractional parts.
#'
#' @param n Total to allocate.
#' @param sizes Group sizes. Everything else is matched to them by position;
#'   names, when present, are carried to the result.
#' @param allocation `"proportional"`, `"equal"` or `"neyman"`.
#' @param min_per_stratum Floor per group. With `cap`, a group smaller than
#'   the floor is taken whole.
#' @param cap When `TRUE`, no group may be allocated more rows than it holds.
#'   `FALSE` when sampling with replacement, where over-allocation is legal.
#' @param spread Per-group standard deviation, in the order of `sizes`, for
#'   Neyman allocation. Groups
#'   get shares proportional to `size * spread`, which minimises the variance
#'   of a total for a fixed `n`: a group whose values barely vary needs fewer
#'   rows to pin down than an equally sized group that varies a lot.
#' @param cost Per-group cost of one row, in the order of `sizes`, or `NULL`.
#'   Neyman shares become `size * spread / sqrt(cost)`.
#' @param noun The group, singular and plural, for the message when the
#'   floor needs more than `n`.
#' @noRd
allocate <- function(n, sizes, allocation, min_per_stratum, cap = TRUE,
                     spread = NULL, min_arg = "min_per_stratum",
                     noun = c("stratum", "strata"), cost = NULL) {
  k <- length(sizes)
  if (k == 0L) return(integer(0))

  # `n` and `sizes` are both integer, so `n * sizes` is evaluated in 32 bits and
  # overflows to NA once the product passes 2^31 -- 3,000 rows out of 2.4
  # million was enough. Everything downstream is fractional anyway.
  n <- as.numeric(n)
  sizes <- stats::setNames(as.numeric(sizes), names(sizes))
  lower <- group_floor(min_per_stratum, sizes, cap)
  if (sum(lower) > n) {
    small <- sum(sizes < min_per_stratum)
    where <- if (k == 1L) {
      paste0(" rows in the one ", noun[1])
    } else {
      paste0(" rows across these ", k, " ", noun[2])
    }
    stop("`", min_arg, "` = ", min_per_stratum, " needs at least ", sum(lower),
         where,
         if (cap && small > 0L) {
           paste0(", with the ", small, " ",
                  if (small == 1L) noun[1] else noun[2], " smaller than ",
                  min_per_stratum, " taken whole")
         }, ", but `n` is ", n, ". Raise `n` or lower `", min_arg, "`.",
         call. = FALSE)
  }

  share <- allocation_share(allocation, sizes, spread, cost, noun)
  raw <- capped_shares(n, share, sizes, cap, lower)

  base <- floor(raw)
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
  }

  stats::setNames(as.integer(base), names(sizes))
}

#' Each group's floor
#'
#' `min_per_stratum` for every group, except that without replacement a
#' group smaller than the floor is taken whole, so it needs only its own
#' rows.
#'
#' @noRd
group_floor <- function(min_per_stratum, sizes, cap) {
  lower <- rep(as.numeric(min_per_stratum), length(sizes))
  if (cap) pmin(lower, as.numeric(sizes)) else lower
}

#' The share of `n` each group gets under an allocation rule
#'
#' Before the floor and the cap. Under Neyman allocation it is
#' `size * spread`, or `size * spread / sqrt(cost)` with costs, the split of
#' least variance for what the sample costs (Cochran 1977, section 5.5).
#' When no group has any spread, every split has the same variance, and the
#' shares are those for equal spreads: `size`, or `size / sqrt(cost)`.
#'
#' @noRd
allocation_share <- function(allocation, sizes, spread = NULL, cost = NULL,
                             noun = c("stratum", "strata")) {
  k <- length(sizes)
  sizes <- as.numeric(sizes)
  switch(allocation,
    equal = rep(1, k),
    neyman = {
      if (is.null(spread)) {
        stop("Neyman allocation needs `allocation_by`.", call. = FALSE)
      }
      sp <- as.numeric(spread)
      if (length(sp) != k || anyNA(sp) || any(!is.finite(sp) | sp < 0)) {
        stop("Neyman allocation needs a finite spread of 0 or more for ",
             "every ", noun[1], ".", call. = FALSE)
      }
      base <- sizes
      if (!is.null(cost)) {
        cs <- as.numeric(cost)
        if (length(cs) != k || anyNA(cs) || any(!is.finite(cs) | cs <= 0)) {
          stop("Neyman allocation needs a finite cost above 0 for every ",
               noun[1], ".", call. = FALSE)
        }
        base <- sizes / sqrt(cs)
      }
      # A constant auxiliary variable has no spread to allocate by; fall back
      # to the shares for equal spreads rather than dividing by zero.
      if (!any(sizes > 0 & sp > 0)) return(base)
      pow2_scale(base) * pow2_scale(sp)
    },
    sizes
  )
}

#' `x` times the power of two that brings its largest value near 1
#'
#' Neyman shares depend only on the ratios of the spreads and of the costs.
#' A product of spreads and costs near the limits of a double can overflow
#' to `Inf` or underflow to 0. Scaled first, each lies between 0 and 2, so
#' their product does neither while its ratios are within the range of a
#' double. A power of two changes no digit of a value in that range, so a
#' split whose products stay in range is the same to the last bit. The
#' scaling is done in two steps because the power can exceed the largest
#' double for a subnormal `x`.
#'
#' @noRd
pow2_scale <- function(x) {
  top <- max(x, 0)
  if (!(top > 0) || !is.finite(top)) return(x)
  e <- -floor(log2(top))
  half <- trunc(e / 2)
  x * 2^half * 2^(e - half)
}

#' Continuous allocation between a floor and a cap
#'
#' Splitting `n` in proportion to `share` can ask a group for more rows than it
#' holds, routinely under Neyman allocation, where a small, highly variable
#' stratum attracts a large share, and under equal allocation. It can also
#' give a group fewer rows than the floor `lower`. The textbook remedy
#' (Cochran 1977, section 5.9) fixes such a group at its bound and re-splits
#' what is left over the others *in the same proportions*. Raising a group to
#' its floor by taking rows from the others one at a time instead spreads the
#' cost evenly and abandons the allocation rule that was asked for.
#'
#' Each group gets `lambda * share`, held between `lower` and its cap, with
#' the one `lambda` that makes the total `n`; `split_level()` finds it. Under
#' Neyman allocation that is the allocation of least variance within the
#' bounds. Groups whose share is 0, Neyman strata whose values do not vary,
#' get their floor. When even taking every group with a share whole leaves
#' rows over, those go to the groups without a share in proportion to their
#' size, as `allocate()` does when no group has a share.
#'
#' @noRd
capped_shares <- function(n, share, sizes, cap,
                          lower = numeric(length(sizes))) {
  upper <- if (cap) sizes else rep(Inf, length(sizes))
  pos <- share > 0
  raw <- lower
  if (sum(upper[pos]) + sum(lower[!pos]) < n) {
    raw[pos] <- upper[pos]
    raw[!pos] <- bounded_split(n - sum(upper[pos]), sizes[!pos], lower[!pos],
                               upper[!pos])
  } else {
    raw[pos] <- bounded_split(n - sum(lower[!pos]), share[pos], lower[pos],
                              upper[pos])
  }
  raw
}

#' What the rule asks of each group, at the level the bounds settle on
#'
#' The `lambda * share` that `capped_shares()` holds between the floor and
#' the cap, before the holding: above a group's size when the cap keeps it
#' from more rows, below its floor when the floor gives it more than its
#' share. Groups with a share get `Inf` when even taking all of them whole
#' leaves rows over, and the groups without one then claim in proportion to
#' size, as `capped_shares()` splits those rows.
#'
#' @noRd
share_claim <- function(n, share, sizes, cap, lower) {
  sizes <- as.numeric(sizes)
  upper <- if (cap) sizes else rep(Inf, length(sizes))
  pos <- share > 0
  claim <- numeric(length(sizes))
  if (sum(upper[pos]) + sum(lower[!pos]) < n) {
    claim[pos] <- Inf
    if (any(!pos)) {
      claim[!pos] <- split_level(n - sum(upper[pos]), sizes[!pos],
                                 lower[!pos], upper[!pos]) * sizes[!pos]
    }
  } else if (any(pos)) {
    claim[pos] <- split_level(n - sum(lower[!pos]), share[pos], lower[pos],
                              upper[pos]) * share[pos]
  }
  claim
}

#' Why each group was allocated the rows it was
#'
#' `"empty"` for a group allocated no rows. Without replacement,
#' `"capped"` for one the rule asked more rows of than it holds, and
#' `"take-all"` for one taken whole for another reason: a share that comes to
#' its size, or a size below the floor. `"floor"` for a group raised to
#' `min_per_stratum` from a smaller share, and `"rule"` for the rest.
#' `share` and `lower` are as `allocate()` computes them.
#'
#' @noRd
allocation_status <- function(n, n_alloc, sizes, share, lower, cap) {
  sizes <- as.numeric(sizes)
  claim <- share_claim(as.numeric(n), share, sizes, cap, lower)
  # The level is found in floating point, so a share that comes to a bound
  # exactly can land a rounding error either side of it.
  tol <- 1e-9 * pmax(1, sizes)
  status <- rep("rule", length(sizes))
  status[claim < lower - tol] <- "floor"
  if (cap) {
    status[n_alloc == sizes] <- "take-all"
    status[claim > sizes + tol] <- "capped"
  }
  status[n_alloc == 0L] <- "empty"
  status
}

#' Split `n` in proportion to positive `share`, between `lower` and `upper`
#'
#' Needs `sum(lower) <= n <= sum(upper)`. The groups held at a bound are found
#' from `split_level()`, and the others share what is left as
#' `left * share / sum(share)`. A split that no bound touches is therefore
#' `n * share / sum(share)` to the last bit, and equal remainders tie.
#'
#' @noRd
bounded_split <- function(n, share, lower, upper) {
  if (!length(share)) return(numeric(0))
  at <- split_level(n, share, lower, upper) * share
  low <- at <= lower
  high <- !low & at >= upper
  free <- !low & !high
  raw <- ifelse(low, lower, upper)
  if (any(free)) {
    left <- n - sum(raw[!free])
    raw[free] <- pmin(pmax(left * share[free] / sum(share[free]),
                           lower[free]), upper[free])
  }
  raw
}

#' The `lambda` at which `sum(pmin(pmax(lambda * share, lower), upper))` is `n`
#'
#' The sum is piecewise linear in `lambda` and never decreasing. A group starts
#' to grow at `lower / share` and stops at `upper / share`. Bisection over
#' those knots finds the last one at which the sum is at most `n`, in time
#' that grows as `k log k`, and the level lies on the segment after it, where
#' the sum is the rows of the groups held at a bound plus `lambda` times the
#' shares of the others. Both are summed afresh over the groups. A running
#' sum of slopes, raised by a share at one knot and lowered by it at
#' another, loses a small share added to a much larger one, and with shares
#' many orders of magnitude apart its level can lie far past `n`.
#'
#' @noRd
split_level <- function(n, share, lower, upper) {
  from <- lower / share
  to <- upper / share
  knot <- c(from, to)
  knot <- sort(unique(knot[is.finite(knot)]), method = "radix")
  if (!length(knot)) return(0)
  total <- function(lambda) sum(pmin(pmax(lambda * share, lower), upper))
  a <- 1L
  b <- length(knot)
  while (a < b) {
    m <- (a + b + 1L) %/% 2L
    if (total(knot[m]) <= n) a <- m else b <- m - 1L
  }
  at <- knot[a]
  free <- from <= at & to > at
  slope <- sum(share[free])
  if (slope <= 0) return(at)
  held <- sum(lower[from > at]) + sum(upper[to <= at])
  max(at, (n - held) / slope)
}

#' Warn when a group with rows in the frame is allocated none of them
#'
#' A stratum allocated zero rows makes every one of its frame rows unreachable:
#' inclusion probability 0, so a Horvitz-Thompson total omits them and is
#' biased by exactly their total. Proportional allocation does this silently to
#' any stratum smaller than `N / n`, and Neyman allocation to any stratum whose
#' auxiliary variable is constant. It is almost never what was meant. An `n`
#' of 0 asks for no rows at all and is not warned about.
#'
#' `flat` marks the groups that Neyman allocation gave no share because the
#' column `by` does not vary within them, or, with `by` `NULL`, because the
#' standard deviations given for them are 0. Raising `n` reaches those only
#' once every other group is taken whole, so their message does not suggest
#' it.
#' `fix` replaces the closing advice, for counts given per stratum.
#'
#' The warning carries its own class so that loops which draw on purpose many
#' times -- the Monte Carlo probabilities -- can silence the repeats.
#'
#' @noRd
warn_empty_groups <- function(n_alloc, sizes, labels, noun = "stratum",
                              arg = "min_per_stratum", flat = NULL,
                              by = NULL, fix = NULL) {
  if (sum(n_alloc) == 0) return(invisible(NULL))
  empty <- which(n_alloc == 0L & sizes > 0)
  if (!length(empty)) return(invisible(NULL))
  plural <- length(empty) > 1L
  nouns <- c(noun, sub("um$", "a", noun))
  flat <- if (is.null(flat)) logical(length(empty)) else flat[empty]
  cause <- ""
  if (any(flat)) {
    some <- !all(flat)
    many <- sum(flat) > 1L
    who <- if (some) quote_labels(labels[empty[flat]]) else if (many) "them"
      else "it"
    cause <- paste0(
      if (is.null(by)) {
        paste0("`allocation_by` gives ", who, " a standard deviation of 0")
      } else {
        paste0("`", by, "` is constant within ", who)
      },
      ", so Neyman allocation gives ", if (many) "them" else "it",
      " no share until every other ", noun, " is taken whole. ")
  }
  fix <- fix %||% paste0(
    if (all(flat)) "Set " else "Raise `n`, or set ", arg, " = 1 (or 2, so ",
    "each ", noun, " also supports a variance).")
  warn_drawn(paste0(
    length(empty), " ", nouns[if (plural) 2L else 1L],
    " (", quote_labels(labels[empty]), ") ",
    if (plural) "are" else "is", " allocated no rows, so ",
    if (plural) "their " else "its ", fmt_n(sum(sizes[empty])),
    " frame row(s) can never be drawn. A total estimated from this sample ",
    "leaves them out entirely. ", cause, fix
  ), "drawn_empty_stratum")
}

#' Up to five group labels for a message, quoted, and how many more
#' @noRd
quote_labels <- function(labels, max = 5L) {
  shown <- paste0("`", group_label(utils::head(labels, max)), "`",
                  collapse = ", ")
  if (length(labels) > max) {
    shown <- paste0(shown, " and ", length(labels) - max, " more")
  }
  shown
}

# ---- per-stratum values given by the user ---------------------------------

#' Per-stratum values as a data frame of the strata columns and `value`
#'
#' `spec` is a named vector, for a single `strata` column, whose names are
#' that column's values as text; or a data frame holding the `strata` columns
#' and the `value` column(s). Either way the result is a plain data frame of
#' those columns and nothing else, one row per stratum. A missing stratum
#' value and a stratum given twice are refused. What the values must be is
#' for the caller to check.
#'
#' @noRd
strata_frame <- function(spec, strata, arg, value = "n") {
  cols <- c(strata, value)
  if (is.data.frame(spec)) {
    if (any(strata %in% value)) {
      stop("`", arg, "` holds its values in a column named ",
           paste0("`", value, "`", collapse = ", "), ", which is also a ",
           "`strata` column. Rename the strata column in `data` and in `",
           arg, "`.", call. = FALSE)
    }
    validate_data(spec, cols, arg)
    check_key_columns(spec, strata, arg)
    out <- lapply(cols, function(col) spec[[col]])
  } else {
    if (length(strata) != 1L) {
      stop("`", arg, "` is a vector, which names one `strata` column's ",
           "values, but `strata` has ", length(strata), " columns. Pass a ",
           "data frame with columns ",
           paste0("`", utils::head(cols, -1L), "`", collapse = ", "),
           " and `", utils::tail(cols, 1L), "` instead.", call. = FALSE)
    }
    if (!length(spec) || is.null(names(spec)) || anyNA(names(spec))) {
      stop("`", arg, "` must be named by stratum, as in ",
           "c(north = 10, south = 5), or be a data frame with columns `",
           strata, "` and `", value[1], "`.", call. = FALSE)
    }
    out <- list(names(spec), unname(spec))
  }
  names(out) <- cols
  out <- structure(out, class = "data.frame",
                   row.names = seq_along(out[[1L]]))

  key <- group_key(out, strata)
  if (anyNA(key)) {
    stop(sum(is.na(key)), " row(s) of `", arg, "` have a missing stratum ",
         "value. Every row must name a stratum.", call. = FALSE)
  }
  twice <- unique(as.integer(key)[duplicated(as.integer(key))])
  if (length(twice)) {
    stop("`", arg, "` gives ",
         if (length(twice) == 1L) "stratum " else "strata ",
         quote_labels(levels(key)[twice]), " more than once. Give each ",
         "stratum once.", call. = FALSE)
  }
  out
}

#' Match per-stratum values given by the user to the strata of a frame
#'
#' `spec` is as for `strata_frame()`. A named vector is matched to the values
#' of the `strata` column as text, as `as.character()` gives them; a data
#' frame is matched on the values themselves, as `group_key2()` compares
#' them, so a number or a date is matched exactly. Returns a data frame of
#' the `value` columns with one row per stratum of `data`, in key order: the
#' levels of `key`, which is `group_key(data, strata)`. A stratum of `data`
#' with no entry is refused by name, so a misspelt stratum cannot quietly
#' get nothing. An entry for a stratum `data` does not have is refused too,
#' unless `ignore_extra`: a count promises rows from its stratum, while a
#' standard deviation or a cost for a stratum with no rows changes nothing,
#' and the same values then serve a part of the frame, such as the `rest`
#' of `design_certainty()`.
#'
#' @noRd
match_strata <- function(spec, data, strata, arg, value = "n",
                         key = group_key(data, strata),
                         ignore_extra = FALSE) {
  spec <- strata_frame(spec, strata, arg, value)
  both <- group_key2(spec, data, strata)
  given <- as.integer(both$x)
  # Each stratum of `data` on the shared codes, read at its first row.
  have <- as.integer(both$table)[match(seq_len(nlevels(key)),
                                       as.integer(key))]
  pos <- match(have, given)
  if (anyNA(pos)) {
    lost <- levels(key)[is.na(pos)]
    stop("`", arg, "` has no entry for ", length(lost), " ",
         if (length(lost) == 1L) "stratum" else "strata", " of `data`: ",
         quote_labels(lost), ". Give one for every stratum.", call. = FALSE)
  }
  extra <- which(!(given %in% have))
  if (length(extra) && !ignore_extra) {
    stop("`", arg, "` names ", length(extra), " ",
         if (length(extra) == 1L) "stratum" else "strata",
         " that `data` does not have: ",
         quote_labels(levels(both$x)[given[extra]]), ". Check the ",
         "spelling, or drop ", if (length(extra) == 1L) "it" else "them", ".",
         call. = FALSE)
  }
  out <- spec[pos, value, drop = FALSE]
  row.names(out) <- NULL
  out
}

#' Refuse a count for a stratum whose every row `na_rm` drops
#'
#' A spread or weighted design splits `n` over the strata of the rows it
#' keeps, so a stratum `na_rm` empties is not one of them, and
#' `match_strata()` would refuse its count as one for a stratum `data` does
#' not have. This names the reason instead. `key` is `group_key()` of the
#' `strata` over every row of `data`, `keep` the rows kept, and `what` why
#' the others were dropped, as in "a missing weight in `w`".
#'
#' @noRd
refuse_emptied_counts <- function(design, data, key, keep, what) {
  if (!has_counts(design)) return(invisible(NULL))
  code <- as.integer(key)
  k <- nlevels(key)
  emptied <- which(tabulate(code, k) > 0L & tabulate(code[keep], k) == 0L)
  if (!length(emptied)) return(invisible(NULL))
  spec <- strata_frame(design$allocation, design$strata, "allocation")
  both <- group_key2(spec, data, design$strata)
  first <- match(emptied, code)
  named <- emptied[as.integer(both$table)[first] %in% as.integer(both$x)]
  if (!length(named)) return(invisible(NULL))
  one <- length(named) == 1L
  stop("`allocation` gives ",
       if (one) "a count for stratum " else "counts for strata ",
       quote_labels(levels(key)[named]), ", but na_rm = TRUE drops every ",
       "row of ", if (one) "it" else "them", " for ", what, ". Fill in ",
       "those values, or leave ", if (one) "the stratum" else "those strata",
       " out of `allocation`.", call. = FALSE)
}
