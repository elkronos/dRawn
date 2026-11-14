#' Simple random sampling
#'
#' Draws `n` rows uniformly at random.
#'
#' @section Reserves:
#' Fieldwork often finds that some sampled units cannot be reached: a site
#' with no access, a business that has closed, a plot under water. `reserve`
#' draws that many more units with the main sample, as a list to replace them
#' from. [draw()] returns the `n` units of the main sample followed by the
#' `reserve` units, with an integer `.reserve` column: 0 for the main sample,
#' and 1, 2, ... for the reserves, in the order they are used. The list is
#' one random draw of `n + reserve` units in random order, so the main sample
#' is a simple random sample of `n`, the one the design without `reserve`
#' draws with the same seed, and the first `k` units of the list are a
#' simple random sample of `k` for any `k`.
#'
#' Visit every unit of the main sample, and for each one found unavailable,
#' visit the next reserve. [use_reserves()] takes the list and what the
#' visits found, and returns the sample to estimate from: the first `n` units
#' of the list found available. They are a simple random sample of `n` from
#' the units of the frame that are available, and [use_reserves()] weights
#' them as one, with probability `n / (N - u)` when `u` listed units were
#' found unavailable. Replacing units changes what is estimated: the totals
#' and means are those of the units that can be reached, not of the whole
#' frame. "What the estimates describe" in [use_reserves()] has the details.
#'
#' Until then the list is not a sample to estimate from. [inclusion_prob()]
#' and the `.prob` column give every row the main sample's probability,
#' `n / N`, and [ht_total()], [ht_mean()], [as_svydesign()] and
#' [sample_summary()] refuse the list.
#'
#' @param n Number of rows to draw. A single non-negative whole number, or a
#'   plan from [plan_size()].
#' @param replace Sample with replacement? Note that `draw(weights = TRUE)`
#'   declines a with-replacement design: `.prob` would be the chance of
#'   appearing at least once, which does not weight a sample holding
#'   duplicates. See [draw()].
#' @param reserve Number of reserve units to list after the main sample, to
#'   replace units found unavailable in the field. A single non-negative whole
#'   number; `n + reserve` cannot exceed the rows of the frame. Not available
#'   with `replace = TRUE` or `n = 0`. See "Reserves".
#'
#' @return A design object, for use with [draw()].
#'
#' @examples
#' df <- data.frame(id = 1:100)
#' draw(df, design_simple(n = 10), seed = 1)
#' draw(df, design_simple(n = 150, replace = TRUE), seed = 1)
#'
#' # Ten units and three reserves, used in the order of `.reserve`
#' draw(df, design_simple(n = 10, reserve = 3), seed = 1)
#'
#' @family designs
#' @seealso [draw()], [use_reserves()]
#' @export
design_simple <- function(n, replace = FALSE, reserve = 0) {
  n <- check_count(plan_n(n), "n")
  replace <- check_flag(replace, "replace")
  reserve <- check_reserve(reserve, n)
  if (replace && reserve > 0L) {
    stop("`reserve` lists units to visit in place of unavailable ones, and ",
         "`replace = TRUE` can draw a unit more than once. Use ",
         "replace = FALSE with `reserve`.", call. = FALSE)
  }
  new_design("simple", list(n = n, replace = replace, reserve = reserve))
}

#' Refuse a simple random sample the frame cannot supply
#'
#' The checks every path runs before it draws or answers, so
#' [inclusion_prob()] and [joint_prob()] refuse what [draw()] refuses, with
#' the same message. Returns the frame rows the design can select: every
#' row, less those use_reserves() recorded as unavailable.
#'
#' @noRd
simple_setup <- function(design, data) {
  validate_data(data)
  keep <- available_rows(design, nrow(data))
  check_draw_size(design$n, length(keep), design$replace)
  check_reserve_size(design, length(keep))
  keep
}

#' @export
draw_design.drawn_design_simple <- function(design, data) {
  keep <- simple_setup(design, data)
  r <- reserve_size(design)
  idx <- sample.int(length(keep), size = design$n, replace = design$replace)
  if (r > 0L) {
    # The reserves are drawn from the rows left, after the main sample, so
    # the main sample is the one the design without `reserve` draws, and the
    # list is one random ordering of `n + reserve` distinct rows.
    left <- setdiff(seq_along(keep), idx)
    idx <- c(idx, left[sample.int(length(left), r)])
  }
  out <- reindex(data, keep[idx])
  if (r > 0L) out <- with_reserve(out, design$n, r)
  out
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_simple <- function(design, data) {
  keep <- simple_setup(design, data)
  N <- length(keep)
  out <- numeric(nrow(data))
  out[keep] <- if (design$replace) {
    # P(selected at least once) over n independent draws.
    1 - (1 - 1 / N)^design$n
  } else {
    design$n / N
  }
  out
}

# ---- reserves ---------------------------------------------------------------

#' The number of reserves a design lists
#'
#' 0 for a design without `reserve`, including one saved before it existed.
#'
#' @noRd
reserve_size <- function(design) {
  design[["reserve"]] %||% 0L
}

#' The same design without its reserves
#'
#' Its draws are the main samples of the design's lists, with the same random
#' numbers, which is what the simulated probabilities and the plots count.
#'
#' @noRd
main_design <- function(design) {
  if (reserve_size(design) > 0L) design[["reserve"]] <- 0L
  design
}

#' Check `reserve`, which needs a main sample to replace units of
#'
#' With `n = 0` there is no unit for a reserve to replace, and the list
#' would hold the reserves alone at probability 0. Returns `reserve` as an
#' integer.
#'
#' @noRd
check_reserve <- function(reserve, n) {
  reserve <- check_count(reserve, "reserve")
  if (reserve > 0L && n == 0) {
    stop("`reserve` is ", reserve, " and `n` is 0, so there is no main ",
         "sample for the reserves to replace units of. Set `n` to 1 or more, ",
         "or leave out `reserve`.", call. = FALSE)
  }
  reserve
}

#' Refuse more rows than the frame holds for the list and its reserves
#' @noRd
check_reserve_size <- function(design, n_rows) {
  r <- reserve_size(design)
  total <- as.numeric(design$n) + r
  if (r > 0L && total > n_rows) {
    stop("`n` and `reserve` together (", total, ") cannot exceed the number ",
         "of rows (", n_rows, "). Lower `reserve`.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Frame rows a design can select, less those found unavailable
#'
#' use_reserves() records the frame rows found unavailable as `unavailable`
#' on the design of the sample it returns, and the design gives them
#' probability 0. A design built by a constructor has none.
#'
#' @noRd
available_rows <- function(design, n_rows) {
  gone <- design[["unavailable"]]
  if (is.null(gone)) return(seq_len(n_rows))
  if (length(gone) && max(gone) > n_rows) {
    stop("`design` records units found unavailable in rows up to ",
         max(gone), ", and `data` has ", n_rows, " rows. It is the design of ",
         "a sample from use_reserves(), and holds for that sample's frame ",
         "only.", call. = FALSE)
  }
  which(!seq_len(n_rows) %in% gone)
}

#' Put the `.reserve` column in front of a drawn list
#'
#' `out` holds the `n` rows of the main sample and then the `r` reserves in
#' the order they are used. The column is built in place and moved to the
#' front by position, as `.replicate` is, so the frame keeps its class.
#'
#' @noRd
with_reserve <- function(out, n, r) {
  if (".reserve" %in% names(out)) {
    stop("`data` already has a column called `.reserve`, which draw() adds ",
         "to number the reserves of a design with `reserve`. Rename it.",
         call. = FALSE)
  }
  k <- ncol(out)
  nm <- names(out)
  out[[".reserve"]] <- c(integer(n), seq_len(r))
  out <- out[, c(k + 1L, seq_len(k)), drop = FALSE]
  names(out) <- c(".reserve", nm)
  out
}
