# The state a weighted sample carries, and the one place it is read.
#
# draw(weights = TRUE) attaches four attributes: the design (`drawn_design`),
# the frame rows in draw order (`drawn_rows`), a copy of the frame
# (`drawn_population`) and every frame row's inclusion probability at draw
# time (`drawn_pi`). It also adds a `.row` column holding each sampled row's
# frame row. Everything that estimates from a sample reads that state through
# drawn_state(), which matches rows through `.row` rather than by position.

#' Copy a data frame column by column
#'
#' A data frame holds references to its column vectors, and data.table edits
#' those vectors in place (`:=`, `set()`, `setorder()`, `setkey()`). A frame
#' stored by reference therefore changes when the caller's does. Every column
#' here is a new vector with the original's attributes, so the copy keeps the
#' frame's class and attributes (tibble, data.table, sf) and later edits to
#' the original do not reach it.
#'
#' @noRd
snapshot_frame <- function(data) {
  out <- lapply(unclass(data), copy_column)
  attributes(out) <- attributes(data)
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
  x <- sort(unique(x))
  paste0(paste(utils::head(x, 5L), collapse = ", "),
         if (length(x) > 5L) ", ..." else "")
}

#' The design, frame and frame rows behind a weighted sample
#'
#' Returns `design`, `pop` (the frame as it was at draw time), `rows` (the
#' frame row of each sample row, in the sample's current order) and `pi_all`
#' (every frame row's inclusion probability at draw time, `NULL` for a sample
#' saved by drawn 0.2.0).
#'
#' `rows` comes from the `.row` column, so reordering a sample, by `[`,
#' dplyr::arrange() or data.table::setorder(), changes nothing. A subset, a
#' sample bound to itself or another sample, or a join that repeats rows does
#' not match the rows that were drawn and is refused. A zero-row sample passes;
#' callers decide what to do with it.
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
    return(list(design = design, pop = pop, rows = drawn, pi_all = NULL))
  }

  list(design = design, pop = pop,
       rows = sample_rows(sample, drawn, nrow(pop), arg), pi_all = pi_all)
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
