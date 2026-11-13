# The state a weighted sample carries, and the one place it is read.
#
# draw(weights = TRUE) attaches four attributes: the design (`drawn_design`),
# the frame rows in draw order (`drawn_rows`), a copy of the frame
# (`drawn_population`) and every frame row's inclusion probability at draw
# time (`drawn_pi`). It also adds a `.row` column holding each sampled row's
# frame row. Everything that estimates from a sample reads that state through
# drawn_state(), which matches rows through `.row` rather than by position,
# checks `.prob` against the probabilities recorded at draw time, and checks
# that the design still gives the frame those probabilities. Design keys
# (strata, clusters, time) are read from the frame's copy of each sampled row,
# never from the sample's own columns, which the caller can rename or edit.

#' Copy a data frame column by column
#'
#' A data frame holds references to its column vectors, and data.table edits
#' those vectors in place (`:=`, `set()`, `setorder()`, `setkey()`). A frame
#' stored by reference therefore changes when the caller's does. Every column
#' here is a new vector with the original's attributes, so the copy keeps the
#' frame's class and attributes (tibble, data.table, sf) and later edits to
#' the original do not reach it. The frame's own attributes are copied too,
#' because `setnames()`, `setcolorder()` and deleting a column with `:=`
#' rewrite the names vector in place, and a shared one renamed or dropped
#' columns of the copy.
#'
#' @noRd
snapshot_frame <- function(data) {
  out <- lapply(unclass(data), copy_column)
  attributes(out) <- lapply(attributes(data), function(a) {
    if (is.atomic(a)) copy_column(a) else a
  })
  out
}

#' @noRd
copy_column <- function(x) {
  if (is.data.frame(x)) return(snapshot_frame(x))
  if (!is.atomic(x) && !is.list(x)) return(x)
  # Without its class, a column has no `[` or length() method that could
  # change what is copied, so every element is copied. The attributes are
  # then put back as they were.
  base <- unclass(x)
  out <- .subset(base, seq_along(base))
  attributes(out) <- attributes(x)
  out
}

#' Refuse a sample whose design cannot be recovered
#'
#' With no `reason`, the sample has lost the attributes that hold its design.
#' Otherwise `reason` completes a sentence that starts with the argument name.
#'
#' @noRd
stop_not_drawn <- function(reason = NULL, arg = "sample") {
  if (is.null(reason)) {
    reason <- paste0(
      "does not carry its design. It was drawn without `weights = TRUE`, or ",
      "a later step dropped the attributes that hold the design; merge(), ",
      "column selection (`s[, cols]`), subset(), transform() and cbind() all ",
      "do. Estimate from the result of draw(..., weights = TRUE), use `by` ",
      "for groups, and add fieldwork columns to that result with ",
      "`s$col <- value` (matched on an id with match() if they come in ",
      "another order).")
  }
  stop("`", arg, "` ", reason, call. = FALSE)
}

#' Up to five values, for a message
#' @noRd
format_rows <- function(x) {
  x <- sort(unique(x), na.last = TRUE)
  paste0(paste(utils::head(x, 5L), collapse = ", "),
         if (length(x) > 5L) ", ..." else "")
}

#' The design, frame and frame rows behind a weighted sample
#'
#' Returns `design`, `pop` (the frame as it was at draw time), `rows` (the
#' frame row of each sample row, in the sample's current order), `pi_all`
#' (every frame row's inclusion probability at draw time, `NULL` for a sample
#' saved by drawn 0.2.0) and `units` (`pop[rows, ]`, the frame's copy of each
#' sample row).
#'
#' `rows` comes from the `.row` column, so reordering a sample, by `[`,
#' dplyr::arrange() or data.table::setorder(), changes nothing. A subset, a
#' sample bound to itself or another sample, or a join that repeats rows does
#' not match the rows that were drawn and is refused. A zero-row sample passes;
#' callers decide what to do with it.
#'
#' Strata, clusters, time stamps and every other design key are read from
#' `units`, so renaming, recoding or dropping one of those columns in the
#' sample changes nothing. `y` and `by` come from the sample.
#'
#' A sample with `drawn_pi` is also refused when its `.prob` column no longer
#' holds the probabilities recorded at draw time, and when the design, applied
#' to the stored frame in this session, gives the frame other probabilities
#' than it did then. The joint probabilities and variance are recomputed in
#' this session, so they would otherwise describe another design. Samples saved
#' by drawn 0.2.0 get neither check: they record no probabilities to compare.
#'
#' @noRd
drawn_state <- function(sample, arg = "sample") {
  if (!is.data.frame(sample)) {
    stop_not_drawn(paste0("must be a data frame returned by ",
                          "draw(..., weights = TRUE), not a ",
                          class(sample)[1], "."), arg)
  }
  design <- attr(sample, "drawn_design", exact = TRUE)
  drawn <- attr(sample, "drawn_rows", exact = TRUE)
  pop <- attr(sample, "drawn_population", exact = TRUE)
  pi_all <- attr(sample, "drawn_pi", exact = TRUE)
  if (is.null(design) || is.null(drawn) || is.null(pop)) {
    stop_not_drawn(arg = arg)
  }

  if (is.null(pi_all)) {
    # Saved by drawn 0.2.0, which had no `.row` column. Rows can only be
    # matched by position, and only a change in the row count shows.
    if (length(drawn) != nrow(sample)) {
      stop_not_drawn(paste0(
        "has ", nrow(sample), " rows but was drawn with ", length(drawn),
        ", so its rows no longer match the frame rows it records. Estimate ",
        "from the sample as draw() returned it, and use `by` for groups."),
        arg)
    }
    return(list(design = design, pop = pop, rows = drawn, pi_all = NULL,
                units = pop[drawn, , drop = FALSE]))
  }

  rows <- sample_rows(sample, drawn, nrow(pop), arg)
  check_prob(sample, pi_all[rows], arg)
  check_drift(design, pop, pi_all, arg)
  list(design = design, pop = pop, rows = rows, pi_all = pi_all,
       units = pop[rows, , drop = FALSE])
}

#' Do two vectors of probabilities agree, element by element?
#'
#' To a relative tolerance of `sqrt(.Machine$double.eps)`, the one
#' `all.equal()` uses. Zeros must match exactly, and so must `NA`s.
#'
#' @noRd
same_prob <- function(a, b) {
  both_na <- is.na(a) & is.na(b)
  close <- abs(a - b) <= sqrt(.Machine$double.eps) * pmax(abs(a), abs(b))
  both_na | (!is.na(close) & close)
}

#' Refuse a `.prob` column that is not the one draw() set
#'
#' ht_total() and ht_mean() divide by `.prob`, and as_svydesign() gives it to
#' `survey` as the probabilities of a PPS design. A value changed after the
#' draw changes those estimates while the joint probabilities and the variance
#' still follow the design. `expected` is the draw-time probability of each
#' sample row, or `NULL` for a sample saved by drawn 0.2.0, which records none;
#' then only the column's type and range are checked.
#'
#' @noRd
check_prob <- function(sample, expected = NULL, arg = "sample") {
  fix <- paste0(" `.prob` holds the design's inclusion probabilities and ",
                "must stay as draw() set it. Put adjusted weights in ",
                "`.weight` and estimate from them with as_svydesign().")
  p <- sample[[".prob"]]
  if (is.null(p)) {
    stop_not_drawn(paste0("has no `.prob` column, so it was removed after ",
                          "draw().", fix), arg)
  }
  if (!is.numeric(p)) {
    stop_not_drawn(paste0("has a ", class(p)[1], " `.prob` column, so it ",
                          "was changed after draw().", fix), arg)
  }
  # A probability of 1 can carry rounding error from a sum, so the upper
  # bound allows for it.
  bad <- is.na(p) | p <= 0 | p > 1 + 1e-9
  if (any(bad)) {
    vals <- unique(signif(p[bad], 4))
    stop_not_drawn(paste0(
      "has `.prob` values that are not probabilities in (0, 1] in ",
      sum(bad), " row(s): ", paste(utils::head(vals, 5L), collapse = ", "),
      if (length(vals) > 5L) ", ...", ". `.prob` was changed after draw().",
      fix), arg)
  }
  if (!is.null(expected)) {
    off <- !same_prob(p, expected)
    if (any(off)) {
      stop_not_drawn(paste0(
        "has `.prob` values in ", sum(off), " row(s) that differ from the ",
        "inclusion probabilities recorded when it was drawn, so `.prob` was ",
        "changed after draw().", fix), arg)
    }
  }
  invisible(p)
}

#' Refuse a sample whose design gives its frame other probabilities now
#'
#' The variance, the design's degrees of freedom and the `survey` mapping are
#' all recomputed from the design and the stored frame in the current session.
#' If the design no longer gives the frame the probabilities it gave at draw
#' time, those figures describe another design. A spatial design decides which
#' rows lie inside its region under the current [sf::sf_use_s2()], so turning
#' spherical geometry on or off between drawing and estimating is enough.
#' Reading a sample saved with saveRDS() into a new session, where s2 is on
#' again, does the same.
#'
#' @noRd
check_drift <- function(design, pop, pi_all, arg = "sample") {
  now <- exact_inclusion(design, pop)
  differ <- if (length(now) != length(pi_all)) {
    rep(TRUE, max(length(now), length(pi_all)))
  } else {
    !same_prob(now, pi_all)
  }
  if (!any(differ)) return(invisible(NULL))

  fix <- if (uses_s2(design)) {
    paste0("sf::sf_use_s2() decides which rows of a spatial design fall ",
           "inside `region`, and it is ", isTRUE(sf::sf_use_s2()), " now. ",
           "Set it to the value it had when the sample was drawn, then ",
           "estimate.")
  } else {
    paste0("Estimate in a session set up as the one the sample was drawn in, ",
           "with the same version of drawn, or draw the sample again.")
  }
  stop_not_drawn(paste0(
    "has a design that gives ", sum(differ), " of its ", length(differ),
    " frame rows a different inclusion probability in this session or this ",
    "version of drawn than when the sample was drawn, so its standard errors ",
    "would describe another design. ", fix), arg)
}

#' Does a design decide membership with `sf`, directly or through `rest`?
#' @noRd
uses_s2 <- function(design) {
  inherits(design, "drawn_design_spatial") ||
    (inherits(design, "drawn_design_certainty") && uses_s2(design$rest))
}

#' Frame rows from the `.row` column, checked against the rows drawn
#' @noRd
sample_rows <- function(sample, drawn, N, arg) {
  ids <- sample[[".row"]]
  if (is.null(ids)) {
    stop_not_drawn(paste0(
      "has no `.row` column, so its rows cannot be matched to the frame rows ",
      "they were drawn from. Keep the `.prob`, `.weight` and `.row` columns ",
      "that draw(..., weights = TRUE) adds."), arg)
  }
  if (!is.numeric(ids)) {
    stop_not_drawn(paste0(
      "has a ", class(ids)[1], " `.row` column, where draw() puts the frame ",
      "row of each sampled row. Estimate from the sample as draw() returned ",
      "it."), arg)
  }
  bad <- is.na(ids) | ids != trunc(ids) | ids < 1 | ids > N
  if (any(bad)) {
    stop_not_drawn(paste0(
      "has `.row` values that are not rows of the ", N, "-row frame it was ",
      "drawn from: ", format_rows(ids[bad]), ". Estimate from the sample as ",
      "draw() returned it."), arg)
  }
  ids <- as.integer(ids)
  dup <- ids[duplicated(ids)]
  if (length(dup)) {
    stop_not_drawn(paste0(
      "holds frame row(s) ", format_rows(dup), " more than once, as rbind() ",
      "of samples or a join that matches a row twice produces. Estimate from ",
      "the sample as draw() returned it, and add fieldwork columns with one ",
      "value per sampled row."), arg)
  }
  extra <- setdiff(ids, drawn)
  if (length(extra)) {
    stop_not_drawn(paste0(
      "holds frame row(s) ", format_rows(extra), " that were not in the ",
      "sample drawn, as rbind() of two samples produces. Estimate from each ",
      "sample on its own."), arg)
  }
  if (length(ids) < length(drawn)) {
    stop_not_drawn(paste0(
      "has ", length(ids), " of the ", length(drawn), " rows that were ",
      "drawn. A subset of a sample is not a sample from the design, and its ",
      "standard errors would be wrong. Estimate from the full result of ",
      "draw() and use `by` for groups."), arg)
  }
  ids
}
