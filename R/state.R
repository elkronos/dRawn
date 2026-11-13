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
#
# A sample weighted conditional on its clusters (see conditional_design())
# also carries the codes of the clusters selected (`drawn_selection`), and
# its `drawn_pi` holds the probabilities given that selection: `NA` for the
# rows of the clusters not selected, whose probability the selection does
# not fix. Such a row could have been drawn, so reachable() counts it.
#
# as_drawn() sets the same columns and attributes on rows read back from a
# file, through rebuild_sample(), which every function that makes a weighted
# sample other than draw() uses. compact_sample() replaces the stored frame
# with its design columns and records what it dropped (`drawn_compact`).
#
# A design with `reserve` gives a list rather than a sample: the main sample
# and then the reserves, with `drawn_rows` in that order, which is the order
# the reserves are used in. drawn_state() refuses it unless the caller reads
# lists, and use_reserves() turns it into a sample.

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
      "add_fieldwork() or `s$col <- value`. A sample read back from a file ",
      "or a database is rebuilt with as_drawn().")
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
#' saved by drawn 0.2.0), `units` (`pop[rows, ]`, the frame's copy of each
#' sample row), `selection` (the codes of the clusters selected, for a
#' sample weighted conditional on them, else `NULL`), `compact` (what
#' compact_sample() dropped from the stored frame, else `NULL`) and `reserve`
#' (for a list from a design with `reserve`, each row's place among the
#' reserves, 0 for the main sample, else `NULL`). In a conditional sample,
#' `pi_all` is `NA` for the rows of clusters not selected; count the rows a
#' design can reach with reachable().
#'
#' A list from a design with `reserve` is refused unless `reserves` is
#' `TRUE`, since it is not a sample to estimate from until use_reserves()
#' has replaced the units found unavailable. Its `.reserve` column, when it
#' has one, must agree with the order the reserves were drawn in.
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
#' A conditional sample's probabilities are worked out again from its stored
#' selection.
#'
#' @noRd
drawn_state <- function(sample, arg = "sample", reserves = FALSE) {
  if (!is.data.frame(sample)) {
    stop_not_drawn(paste0("must be a data frame returned by ",
                          "draw(..., weights = TRUE), not a ",
                          class(sample)[1], "."), arg)
  }
  design <- attr(sample, "drawn_design", exact = TRUE)
  drawn <- attr(sample, "drawn_rows", exact = TRUE)
  pop <- attr(sample, "drawn_population", exact = TRUE)
  pi_all <- attr(sample, "drawn_pi", exact = TRUE)
  compact <- attr(sample, "drawn_compact", exact = TRUE)
  if (is.null(design) || is.null(drawn) || is.null(pop)) {
    stop_not_drawn(arg = arg)
  }
  r <- reserve_size(design)
  if (r > 0L && !reserves) {
    stop_not_drawn(paste0(
      "lists the ", design$n, " units of a main sample and ", r, " reserve",
      if (r > 1L) "s", " (`reserve = ", r, "`), so it is not yet a sample to ",
      "estimate from. Pass the list to use_reserves() with the units found ",
      "unavailable marked, and estimate from the sample it returns."), arg)
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
                units = pop[drawn, , drop = FALSE], selection = NULL,
                compact = compact))
  }

  rows <- sample_rows(sample, drawn, nrow(pop), arg)
  check_prob(sample, pi_all[rows], arg)
  selection <- check_selection(sample, design, pop, arg)
  check_drift(design, pop, pi_all, arg, selection)
  list(design = design, pop = pop, rows = rows, pi_all = pi_all,
       units = pop[rows, , drop = FALSE], selection = selection,
       compact = compact,
       reserve = if (r > 0L) list_places(sample, rows, drawn, design$n, arg))
}

#' Each row's place among a list's reserves, checked against `.reserve`
#'
#' `drawn` holds the frame rows of the main sample and then of the reserves
#' in the order they are used, so a row's place in it, less `n`, is its
#' reserve number, and 0 for a row of the main sample.
#'
#' @noRd
list_places <- function(sample, rows, drawn, n, arg = "sample") {
  place <- pmax(match(rows, drawn) - n, 0L)
  col <- sample[[".reserve"]]
  if (!is.null(col)) {
    ok <- is.numeric(col) && is.null(dim(col)) && !anyNA(col) &&
      all(col == place)
    if (!ok) {
      stop_not_drawn(paste0(
        "has `.reserve` values that differ from the order its reserves were ",
        "drawn in, so `.reserve` was changed after draw(). use_reserves() ",
        "uses the reserves in the order drawn; keep `.reserve` as draw() set ",
        "it."), arg)
    }
  }
  as.integer(place)
}

#' Frame rows a sample's design could have drawn
#'
#' Rows with a draw-time probability above 0, and the rows of the clusters a
#' conditional sample did not select, whose probability is `NA` but which
#' another selection could have drawn.
#'
#' @noRd
reachable <- function(pi_all) {
  is.na(pi_all) | pi_all > 0
}

#' The clusters a conditional sample selected, checked against its frame
#'
#' `NULL` for a design that is not weighted conditional on its clusters,
#' whose probabilities do not depend on them. Otherwise the stored codes
#' must be `n_clusters` distinct clusters of the frame; a sample that lost
#' them cannot have its probabilities checked or its variance worked out.
#'
#' @noRd
check_selection <- function(sample, design, pop, arg = "sample") {
  if (!conditional_design(design, pop)) return(NULL)
  selection <- attr(sample, "drawn_selection", exact = TRUE)
  fix <- paste0(" `.prob` is conditional on the clusters draw() selected, ",
                "which it records with the sample. Estimate from the sample ",
                "as draw() returned it.")
  if (is.null(selection)) {
    stop_not_drawn(paste0("does not record which clusters were selected, ",
                          "so its probabilities cannot be checked.", fix),
                   arg)
  }
  total <- count_clusters(design, pop)$total
  ok <- is.numeric(selection) && !anyNA(selection) &&
    length(selection) == design$n_clusters &&
    all(selection == trunc(selection) & selection >= 1 &
          selection <= total) &&
    !anyDuplicated(selection)
  if (!ok) {
    stop_not_drawn(paste0("records clusters that are not ",
                          design$n_clusters, " distinct clusters of the ",
                          total, " in its frame.", fix), arg)
  }
  as.integer(selection)
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
#' again, does the same. A sample weighted conditional on its clusters is
#' checked against the probabilities given the stored `selection`.
#'
#' @noRd
check_drift <- function(design, pop, pi_all, arg = "sample",
                        selection = NULL) {
  now <- if (is.null(selection)) {
    exact_inclusion(design, pop)
  } else {
    conditional_inclusion(design, pop, selection)
  }
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

# ---- rebuilding, joining and compacting a sample ------------------------------

#' Rebuild a weighted sample from its rows
#'
#' A weighted sample carries its design, a copy of its frame and its
#' inclusion probabilities as attributes. A CSV or spreadsheet file, a
#' database table and [merge()] keep none of them, so [ht_total()],
#' [ht_mean()], [as_svydesign()] and [sample_summary()] refuse the rows that
#' come back. `as_drawn()` puts the attributes back. It matches each row of
#' `data` to its row of `frame`, works out the inclusion probabilities again
#' from `design`, and sets `.prob`, `.weight` and `.row` as [draw()] sets
#' them.
#'
#' @section Matching rows to the frame:
#' With `by = NULL`, each row is matched through its `.row` column, the frame
#' row [draw()] recorded for it. `frame` must then be the frame the sample
#' was drawn from, with its rows in the same order. With `by`, rows are
#' matched on the values of one or more id columns that `data` and `frame`
#' share. Values are compared as [match()] compares them, and as text when
#' the two columns have different classes, so an id stored as `"007"` does
#' not match the number 7 that a CSV file gives back.
#'
#' When `data` has both `.row` and the `by` columns, they must agree. A frame
#' row found through `by` that is not the one `.row` names means `frame` is
#' not in the order it had when the sample was drawn, which changes the
#' probabilities of designs that depend on that order, such as
#' [design_systematic()]; the rebuild is refused.
#'
#' A row with a missing id, an id that appears twice in `data`, an id that
#' `frame` does not have or holds more than once, and a row the design gives a
#' probability of 0 are refused, with the ids named.
#'
#' @section Every row drawn:
#' `data` must hold every row that was drawn, with or without a fieldwork
#' result. A subset, such as the rows that responded, is not a sample from
#' the design, and its standard errors would be wrong. For a design that
#' draws a fixed number of rows (simple, stratified, spread, temporal,
#' spatial and reservoir designs, systematic designs given `n`, and
#' [design_weighted()] with `method = "systematic"`, alone or as the `rest`
#' of [design_certainty()]), a `data` with another number of rows is
#' refused. A design weighted conditional on its clusters (see "Design
#' weights" in [draw()]) is rebuilt from the clusters its rows belong to, and
#' refused when they are not `n_clusters` clusters or a cluster does not hold
#' the rows the design takes from it. Keep the rows with no response, with
#' `y` missing, and estimate a mean with `ht_mean(na_rm = TRUE)`; see
#' "Missing values and nonresponse" in [ht_total()].
#'
#' A list from a design with `reserve` is rebuilt whole, the main sample and
#' every reserve, with the `.reserve` column [draw()] added, which records
#' the order the reserves are used in. Pass the result to [use_reserves()].
#'
#' @section What comes back:
#' The rows of `data` in their own order, with `.prob`, `.weight` and `.row`
#' first and the other columns of `data` after them. `.prob` and `.weight`
#' come from `design`. A numeric `.prob` in `data` must agree with them, or
#' the design or the frame is not the one the sample was drawn with. A
#' `.weight` adjusted after the draw is replaced by `1 / .prob`, so the
#' adjustment has to be made again.
#'
#' Strata, clusters and the other design columns are read from the copy of
#' `frame` stored with the result, as for a sample from [draw()], so `data`
#' needs only the columns to estimate from: `y` and the `by` columns of the
#' estimates. [add_fieldwork()] joins results that come back keyed by an id.
#'
#' @param data The sampled rows as read back: a data frame with the `.row`
#'   column [draw()] added, or with the id column or columns named in `by`.
#' @param design The design the sample was drawn with.
#' @param frame The data frame the sample was drawn from, as it was then.
#' @param by `NULL` to match rows through `.row`, or the names of one or more
#'   columns of `data` and `frame` that identify a frame row.
#'
#' @return A data frame of the class of `data`, ready for [ht_total()],
#'   [ht_mean()], [as_svydesign()] and [sample_summary()].
#'
#' @examples
#' pop <- data.frame(id = 1:200, site = rep(c("a", "b"), times = c(150, 50)))
#' by_site <- design_stratified("site", n = 40)
#' s <- draw(pop, by_site, seed = 1, weights = TRUE)
#'
#' # A file keeps the rows and loses the design
#' path <- tempfile(fileext = ".csv")
#' write.csv(s, path, row.names = FALSE)
#' back <- as_drawn(read.csv(path), by_site, pop)
#' sample_summary(back)
#'
#' # Without `.row`, match on an id column
#' ids <- read.csv(path)[, c("id", "site")]
#' back <- as_drawn(ids, by_site, pop, by = "id")
#' unlink(path)
#'
#' @seealso [add_fieldwork()] to join fieldwork results, [compact_sample()]
#'   to drop frame columns before saving or sharing a sample, and [draw()].
#' @export
as_drawn <- function(data, design, frame, by = NULL) {
  if (!is_design(design)) {
    stop("`design` must come from one of the design_*() constructors, not a ",
         class(design)[1], ". See ?designs.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame of the sampled rows, not a ",
         class(data)[1], ".", call. = FALSE)
  }
  if (!is.data.frame(frame)) {
    stop("`frame` must be the data frame the sample was drawn from, not a ",
         class(frame)[1], ".", call. = FALSE)
  }
  if (samples_with_replacement(design)) {
    stop("`design` samples with replacement, and draw() gives such a sample ",
         "no weights to rebuild; see \"Design weights\" in ?draw.",
         call. = FALSE)
  }
  clash <- intersect(c(".prob", ".weight", ".row"), names(frame))
  if (length(clash)) {
    stop("`frame` has ", quote_names(clash), ", which a weighted sample has ",
         "and a frame does not. Pass the data frame the sample was drawn ",
         "from as `frame`, and the sample as `data`.", call. = FALSE)
  }
  if (anyDuplicated(names(data))) {
    dup <- unique(names(data)[duplicated(names(data))])
    stop("`data` has duplicated column name(s): ", format_bad(dup),
         ". Rename them before rebuilding the sample.", call. = FALSE)
  }
  warn_region(design, frame)

  rows <- frame_rows(data, frame, by)
  selection <- NULL
  if (conditional_design(design, frame)) {
    selection <- rebuilt_selection(design, frame, rows)
    pi_all <- conditional_inclusion(design, frame, selection)
  } else {
    pi_all <- exact_inclusion(design, frame)
    if (inherits(design, "drawn_design_certainty")) {
      check_rest_reach(design, frame, pi_all)
    }
  }
  p <- pi_all[rows]
  zero <- is.na(p) | p <= 0
  if (any(zero)) {
    stop("`data` holds ", sum(zero), " row(s) that `design` gives a ",
         "probability of 0 in `frame`, so it cannot have drawn them: frame ",
         "row(s) ", format_rows(rows[zero]), ". Pass the design and the ",
         "frame the sample was drawn with.", call. = FALSE)
  }
  r <- reserve_size(design)
  if (r > 0L) {
    # A list: the main sample and its reserves, whose order comes from
    # `.reserve`.
    listed <- as.numeric(design$n) + r
    if (nrow(data) != listed) {
      stop("`data` has ", nrow(data), " rows, and `design` lists ",
           listed, ": ", design$n, " in the main sample and ", r,
           " reserve", if (r > 1L) "s", ". Rebuild the whole list as draw() ",
           "returned it, then pass it to use_reserves().", call. = FALSE)
    }
    reserve <- rebuilt_reserve(data, design$n, r)
  } else if (fixed_size(design)) {
    expected <- sum(pi_all)
    if (abs(nrow(data) - expected) > 1e-6 * max(1, expected)) {
      stop("`data` has ", nrow(data), " rows, and `design` draws ",
           format(signif(expected, 6)), " from `frame`. A sample rebuilt ",
           "from some of the rows drawn, such as those that responded, is ",
           "not a sample from the design, and its standard errors would be ",
           "wrong. Keep every row drawn, with `y` missing where there is no ",
           "result, and estimate a mean with ht_mean(na_rm = TRUE).",
           call. = FALSE)
    }
  }
  old <- data[[".prob"]]
  if (is.numeric(old)) {
    off <- !same_prob(old, p)
    if (any(off)) {
      stop("`data` has `.prob` values in ", sum(off), " row(s) that differ ",
           "from the probabilities `design` gives `frame`, such as ",
           format(signif(old[off][1], 4)), " where the design gives ",
           format(signif(p[off][1], 4)), ". Pass the design and the frame ",
           "the sample was drawn with.", call. = FALSE)
    }
  }
  if (r > 0L) {
    return(rebuild_sample(data, design, frame, rows, pi_all, selection,
                          reserve = reserve, drawn = rows[order(reserve)]))
  }
  rebuild_sample(data, design, frame, rows, pi_all, selection)
}

#' Each row's reserve number, from the `.reserve` column of a list
#'
#' A list drawn with `reserve = r` numbers its `n` main rows 0 and its
#' reserves 1 to `r`, the order they are used in, which nothing else in a
#' file records.
#'
#' @noRd
rebuilt_reserve <- function(data, n, r) {
  col <- data[[".reserve"]]
  if (is.null(col)) {
    stop("`data` has no `.reserve` column, which records the order of the ",
         "reserves `design` lists. Keep the `.reserve` column draw() added.",
         call. = FALSE)
  }
  ok <- is.numeric(col) && is.null(dim(col)) && !anyNA(col) &&
    all(col == trunc(col)) &&
    identical(sort(as.integer(col)), c(integer(n), seq_len(r)))
  if (!ok) {
    stop("`data` has a `.reserve` column that does not hold 0 for ", n,
         " rows and each of 1 to ", r, " once, as draw() set it. Restore it ",
         "from the list draw() returned.", call. = FALSE)
  }
  as.integer(col)
}

#' Set the columns and attributes of a weighted sample on rows of a frame
#'
#' `rows` are the frame rows of `data`'s rows, in `data`'s order, and
#' `pi_all` every frame row's inclusion probability, as draw() records them.
#' `selection` holds the codes of the clusters selected, for a design
#' weighted conditional on them. Any `.prob`, `.weight` and `.row` columns and
#' `drawn_*` attributes `data` already has are replaced, so drawn_state()
#' accepts the result on the same terms as a sample from draw().
#'
#' `reserve`, when given, is each row's reserve number, set as a `.reserve`
#' column after `.row`. `drawn` is the frame rows in the order recorded as
#' `drawn_rows`, which for a list is the main sample and then the reserves
#' in the order they are used. `compact` is what compact_sample() dropped
#' from `frame`, kept as `drawn_compact`.
#'
#' @noRd
rebuild_sample <- function(data, design, frame, rows, pi_all,
                           selection = NULL, reserve = NULL, drawn = rows,
                           compact = NULL) {
  added <- c(".prob", ".weight", ".row", if (!is.null(reserve)) ".reserve")
  out <- data
  for (col in intersect(added, names(out))) out[[col]] <- NULL
  p <- pi_all[rows]
  out[[".prob"]] <- p
  out[[".weight"]] <- 1 / p
  out[[".row"]] <- rows
  if (!is.null(reserve)) out[[".reserve"]] <- as.integer(reserve)
  out <- out[, c(added, setdiff(names(out), added)), drop = FALSE]
  for (a in grep("^drawn_", names(attributes(out)), value = TRUE)) {
    attr(out, a) <- NULL
  }
  attr(out, "drawn_design") <- design
  # `.row` holds `rows` itself, so the attribute gets its own copy, as in
  # draw().
  attr(out, "drawn_rows") <- copy_column(drawn)
  attr(out, "drawn_population") <- snapshot_frame(frame)
  attr(out, "drawn_pi") <- pi_all
  if (!is.null(selection)) attr(out, "drawn_selection") <- selection
  if (!is.null(compact)) attr(out, "drawn_compact") <- compact
  table_ready(out)
}

#' The frame row of each row of `data`, through `.row` or the `by` columns
#' @noRd
frame_rows <- function(data, frame, by) {
  N <- nrow(frame)
  stored <- data[[".row"]]
  if (is.null(by)) {
    if (is.null(stored)) {
      stop("`data` has no `.row` column, which matches each row to the ",
           "frame row it was drawn from. Name an id column of `data` and ",
           "`frame` in `by`.", call. = FALSE)
    }
    if (!is.numeric(stored)) {
      stop("`data` has a ", class(stored)[1], " `.row` column, where draw() ",
           "puts the frame row of each sampled row. Name an id column of ",
           "`data` and `frame` in `by` instead.", call. = FALSE)
    }
    bad <- is.na(stored) | stored != trunc(stored) | stored < 1 | stored > N
    if (any(bad)) {
      stop("`data` has `.row` values that are not rows of the ", N,
           "-row `frame`: ", format_rows(stored[bad]), ". Pass the frame the ",
           "sample was drawn from, or name an id column in `by`.",
           call. = FALSE)
    }
    rows <- as.integer(stored)
    dup <- rows[duplicated(rows)]
    if (length(dup)) {
      stop("`data` holds frame row(s) ", format_rows(dup), " more than ",
           "once. A sample holds each row it drew once; remove the repeats.",
           call. = FALSE)
    }
    return(rows)
  }

  by <- check_columns(by, "by")
  for (where in c("data", "frame")) {
    absent <- setdiff(by, names(if (where == "data") data else frame))
    if (length(absent)) {
      stop("`by` names ", quote_names(absent), ", which `", where, "` does ",
           "not have.", call. = FALSE)
    }
  }
  check_key_columns(data, by, "by")
  check_key_columns(frame, by, "by")
  keys <- group_key2(data, frame, by)
  kd <- as.integer(keys$x)
  kf <- as.integer(keys$table)
  labels <- group_label(levels(keys$x))
  id <- quote_names(by)
  if (anyNA(kd)) {
    stop("`data` has no value of ", id, " in ", sum(is.na(kd)), " row(s), ",
         "so they cannot be matched to `frame`. Every sampled row needs one.",
         call. = FALSE)
  }
  dup <- kd[duplicated(kd)]
  if (length(dup)) {
    stop("`data` has ", id, " ", format_rows(labels[dup]), " in more than ",
         "one row. A sample holds each row it drew once; remove the repeats.",
         call. = FALSE)
  }
  rows <- match(kd, kf)
  if (anyNA(rows)) {
    stop("`data` has ", id, " ", format_rows(labels[kd[is.na(rows)]]),
         ", which `frame` does not have. Pass the frame the sample was ",
         "drawn from.", call. = FALSE)
  }
  twice <- intersect(kd, kf[duplicated(kf)])
  if (length(twice)) {
    stop("`frame` has ", id, " ", format_rows(labels[twice]), " in more ",
         "than one row, so ", id, " does not identify a frame row. Name ",
         "columns in `by` that do, or match through `.row` with by = NULL.",
         call. = FALSE)
  }
  if (is.numeric(stored)) {
    off <- is.na(stored) | stored != rows
    if (any(off)) {
      stop("`data` has a `.row` column that disagrees with ", id, " in ",
           sum(off), " row(s), so `frame` does not have its rows in the ",
           "order they had when the sample was drawn, and some designs ",
           "depend on that order. Pass `frame` in that order.", call. = FALSE)
    }
  }
  rows
}

#' The clusters a conditionally weighted sample selected, from its rows
#'
#' A file keeps the rows and not the selection. The clusters are those the
#' rows belong to, and the probabilities given them need each cluster to
#' hold the rows the design takes from it. Which clusters a largest-remainder
#' split gives one row more depends on the order they were selected in, which
#' is lost too. Listing first the clusters holding more rows than the split
#' in code order gives them reproduces the split the draw made.
#'
#' @noRd
rebuilt_selection <- function(design, frame, rows) {
  cl <- conditional_setup(design, frame)
  code <- cl$code[rows]
  chosen <- sort(unique(code[!is.na(code)]))
  if (length(chosen) != design$n_clusters) {
    stop("`data` holds rows of ", length(chosen), " cluster(s), and ",
         "`design` selects ", design$n_clusters, ". A sample rebuilt from ",
         "some of the rows drawn is not a sample from the design; keep every ",
         "row drawn.", call. = FALSE)
  }
  sizes <- tabulate(cl$code, cl$total)
  got <- tabulate(match(code, chosen), length(chosen))
  take <- cluster_takes(design, sizes[chosen])
  if (any(take != got)) {
    chosen <- chosen[order(take - got)]
    got <- tabulate(match(code, chosen), length(chosen))
    take <- cluster_takes(design, sizes[chosen])
  }
  off <- which(take != got)
  if (length(off)) {
    label <- cl$labels[match(chosen[off[1]], cl$code)]
    stop("`data` holds ", got[off[1]], " row(s) of cluster `", format(label),
         "`, where `design` takes ", take[off[1]], " from it given the ",
         "clusters selected. A sample rebuilt from some of the rows drawn is ",
         "not a sample from the design; keep every row drawn.", call. = FALSE)
  }
  chosen
}

#' Does a design draw the same number of rows from a frame every time?
#'
#' Only designs certain to are listed, since as_drawn() refuses a sample of
#' another size from them. Random-size designs, such as Poisson sampling,
#' whole clusters of unequal size or an interval that does not divide the
#' frame, are not.
#'
#' @noRd
fixed_size <- function(design) {
  if (inherits(design, "drawn_design_certainty")) {
    return(fixed_size(design$rest))
  }
  if (inherits(design, "drawn_design_systematic")) {
    return(!is.null(design$n))
  }
  if (inherits(design, "drawn_design_weighted")) {
    return(identical(design$method, "systematic"))
  }
  inherits(design, c("drawn_design_simple", "drawn_design_stratified",
                     "drawn_design_spread", "drawn_design_temporal",
                     "drawn_design_spatial", "drawn_design_reservoir"))
}

#' Join fieldwork results to a weighted sample
#'
#' Fieldwork results usually come back as a table of their own, keyed by an
#' id and in another order. `add_fieldwork()` adds its columns to `sample`,
#' row by row through the id, and keeps the sample's rows, their order and
#' the attributes that hold its design. [merge()] reorders the rows and drops
#' those attributes, and a join that matches a row twice would repeat it.
#'
#' A sampled row with no row in `field` gets missing values, as a
#' nonrespondent should: [ht_total()] refuses them, and `ht_mean(na_rm = TRUE)`
#' estimates the mean over the rows with a result (see "Missing values and
#' nonresponse" in [ht_total()]). An id that appears twice in `field`, an id
#' that is not in `sample`, a missing id, and a column of `field` that
#' `sample` already has are refused. Ids are compared as in [as_drawn()].
#'
#' A list from a design with `reserve` takes fieldwork columns too, such as
#' which units were found unavailable, before [use_reserves()] turns it into
#' the sample to estimate from.
#'
#' @param sample A data frame returned by [draw()] with `weights = TRUE`, or
#'   by [as_drawn()].
#' @param field A data frame of results: the `by` columns and the columns to
#'   add.
#' @param by The names of one or more columns of `sample` and `field` that
#'   identify a sampled row.
#'
#' @return `sample` with the other columns of `field` added after its own.
#'
#' @examples
#' pop <- data.frame(id = 1:200, site = rep(c("a", "b"), times = c(150, 50)))
#' s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
#'
#' # Results for 36 of the 40 rows, in the order they came in
#' field <- data.frame(id = rev(s$id)[1:36], spend = seq(20, 370, by = 10))
#' s <- add_fieldwork(s, field, by = "id")
#' ht_mean(s, "spend", na_rm = TRUE)
#'
#' @seealso [as_drawn()], [ht_mean()], [use_reserves()]
#' @export
add_fieldwork <- function(sample, field, by) {
  drawn_state(sample, reserves = TRUE)
  if (!is.data.frame(field)) {
    stop("`field` must be a data frame of fieldwork results, not a ",
         class(field)[1], ".", call. = FALSE)
  }
  by <- check_columns(by, "by")
  for (where in c("sample", "field")) {
    absent <- setdiff(by, names(if (where == "sample") sample else field))
    if (length(absent)) {
      stop("`by` names ", quote_names(absent), ", which `", where, "` does ",
           "not have.", call. = FALSE)
    }
  }
  check_key_columns(sample, by, "by")
  check_key_columns(field, by, "by")
  new <- setdiff(names(field), by)
  if (!length(new)) {
    stop("`field` has no column besides ", quote_names(by), ", so there is ",
         "nothing to add.", call. = FALSE)
  }
  clash <- intersect(new, names(sample))
  if (length(clash)) {
    them <- if (length(clash) > 1L) "them" else "it"
    stop("`field` has ", quote_names(clash), ", which `sample` already has. ",
         "Rename ", them, " in `field`, or remove ", them, " from `sample` ",
         "with `s$col <- NULL` first.", call. = FALSE)
  }

  keys <- group_key2(sample, field, by)
  ks <- as.integer(keys$x)
  kf <- as.integer(keys$table)
  labels <- group_label(levels(keys$x))
  id <- quote_names(by)
  if (anyNA(ks)) {
    stop("`sample` has no value of ", id, " in ", sum(is.na(ks)), " row(s), ",
         "so `field` cannot be matched to them.", call. = FALSE)
  }
  dup <- ks[duplicated(ks)]
  if (length(dup)) {
    stop("`sample` has ", id, " ", format_rows(labels[dup]), " in more ",
         "than one row, so ", id, " does not identify its rows. Join on ",
         "columns that do.", call. = FALSE)
  }
  if (anyNA(kf)) {
    stop("`field` has no value of ", id, " in ", sum(is.na(kf)), " row(s), ",
         "so they cannot be matched to the sample.", call. = FALSE)
  }
  dup <- kf[duplicated(kf)]
  if (length(dup)) {
    stop("`field` has ", id, " ", format_rows(labels[dup]), " in more than ",
         "one row. Give each sampled row one result.", call. = FALSE)
  }
  unknown <- setdiff(kf, ks)
  if (length(unknown)) {
    stop("`field` has ", id, " ", format_rows(labels[unknown]), ", which ",
         if (length(unknown) > 1L) "are" else "is", " not in `sample`. ",
         "Remove results for rows the sample did not draw.", call. = FALSE)
  }

  at <- match(ks, kf)
  out <- sample
  for (col in new) {
    x <- field[[col]]
    out[[col]] <- if (is.null(dim(x))) x[at] else x[at, , drop = FALSE]
  }
  table_ready(out)
}

#' Replace unavailable units from a reserve list
#'
#' [design_simple()] and [design_spread()] with `reserve` give a list: the
#' `n` units of the main sample, then the reserves, with a `.reserve` column
#' that numbers them 1, 2, ... in the order they are used. Visit every unit
#' of the main sample, and for each one found unavailable, visit the next
#' reserve, until the units found available make up `n`.
#' `use_reserves()` takes the list and what the visits found, and returns the
#' sample to estimate from: the units of the main sample found available,
#' and as many reserves, taken in order, as replace the others.
#'
#' @section What the estimates describe:
#' Replacing units changes the population being estimated. The sample that
#' comes back stands for the units of the frame that can be reached, and its
#' totals and means are theirs, not the whole frame's. It is weighted as a
#' sample of `n` drawn at equal probability from the frame without the `u`
#' listed units found unavailable: each unit used has probability
#' `n / (N - u)`, and the units found unavailable have probability 0, so
#' [ht_total()] and [ht_mean()] count them among the frame rows the estimate
#' does not cover.
#'
#' The probabilities assume that every unit of the frame nobody visited is
#' available. A frame that holds more unavailable units than were found has
#' fewer than `N - u` units that can be reached, and the totals are too high
#' by about the share of the frame those units make up. The Hajek mean of
#' [ht_mean()] does not depend on `N - u` and is not affected.
#'
#' For [design_simple()] the weighting is exact under that assumption, when
#' the reserves suffice. The list is one random draw in random order, so its
#' first `n` available units are a simple random sample of `n` from the
#' available units of the frame, and the variance is that of a simple random
#' sample from `N - u`. For
#' [design_spread()] the reserves are a second spread sample from the rows
#' the main sample left, in random order. The units used are spread over the
#' space, and under the same assumption their probabilities average
#' `n / (N - u)`, but each unit's own probability is not known exactly, so
#' the weights are an approximation. The variance is the local mean
#' estimator over the units used.
#'
#' Reserves are for units that cannot be reached or measured at all. A unit
#' that is reached but gives no answer is nonresponse: keep it in the sample
#' with `y` missing, as "Missing values and nonresponse" in [ht_total()]
#' describes.
#'
#' @section Marking the visits:
#' `unavailable` holds one value for each row of `sample`: `TRUE` for a unit
#' found unavailable, `FALSE` for one found available, and `NA` for one not
#' visited. Every unit of the main sample needs `TRUE` or `FALSE`, and so
#' does each reserve up to the one that completes the sample; a missing
#' value there is refused, with the frame row named. Reserves visited past
#' that one are not used, and those of them found unavailable count among
#' the `u` units outside the population. When the reserves run out first,
#' the sample holds fewer than `n` units and is weighted as a sample of that
#' many, with a warning.
#'
#' @section What comes back:
#' The units used, in the order of `sample`, with `.prob`, `.weight`, `.row`
#' and `.reserve` first. `.reserve` is 0 for a unit of the main sample and
#' the reserve's number for a unit that replaced one. The design stored with
#' the result is the list's design with `n` the number of units used, no
#' `reserve`, and the frame rows found unavailable recorded as
#' `unavailable`; [inclusion_prob()] of it gives those rows probability 0.
#'
#' @param sample A list from [draw()] with `weights = TRUE` and a design
#'   with `reserve`, or from [as_drawn()].
#' @param unavailable A logical vector with one value for each row of
#'   `sample`, in its order, or the name of such a column of `sample`. See
#'   "Marking the visits".
#'
#' @return A weighted sample, ready for [ht_total()], [ht_mean()],
#'   [as_svydesign()] and [sample_summary()].
#'
#' @examples
#' pop <- data.frame(id = 1:200, spend = (1:200) %% 37)
#' s <- draw(pop, design_simple(n = 20, reserve = 5), seed = 1, weights = TRUE)
#' table(s$.reserve)
#'
#' # Two units of the main sample could not be reached, nor could the first
#' # reserve; reserves 2 and 3 replaced them, and 4 and 5 were not visited
#' gone <- c(s$id[c(3, 11)], s$id[s$.reserve == 1])
#' found <- ifelse(s$.reserve <= 3, s$id %in% gone, NA)
#' used <- use_reserves(s, found)
#' table(used$.reserve)
#' ht_total(used, "spend")
#'
#' @seealso [design_simple()], [design_spread()], [add_fieldwork()] to add
#'   what the visits found to the list.
#' @export
use_reserves <- function(sample, unavailable) {
  state <- drawn_state(sample, reserves = TRUE)
  design <- state$design
  r <- reserve_size(design)
  if (r == 0L || is.null(state$reserve)) {
    stop("`sample` was drawn with no `reserve`, so it lists no reserves to ",
         "use. Estimate from it as it is.", call. = FALSE)
  }
  gone <- visit_findings(sample, unavailable)
  place <- state$reserve
  rows <- state$rows
  main <- place == 0L
  open <- main & is.na(gone)
  if (any(open)) {
    stop("`unavailable` is missing for ", sum(open), " unit(s) of the main ",
         "sample, at frame row(s) ", format_rows(rows[open]), ". Visit every ",
         "unit of the main sample and mark it TRUE or FALSE.", call. = FALSE)
  }

  # Each unavailable unit of the main sample is replaced by the next reserve
  # found available, in the order of `.reserve`.
  used <- main & !gone
  need <- sum(gone[main])
  at <- match(seq_len(r), place)
  for (k in seq_len(r)) {
    if (need == 0L) break
    i <- at[k]
    if (is.na(gone[i])) {
      stop("`unavailable` is missing for reserve ", k, ", at frame row ",
           rows[i], ", which is needed to replace ", need, " more ",
           "unavailable unit(s). Visit the reserves in order and mark each ",
           "one TRUE or FALSE.", call. = FALSE)
    }
    if (!gone[i]) {
      used[i] <- TRUE
      need <- need - 1L
    }
  }
  m <- sum(used)
  if (m == 0L) {
    stop("Every unit `sample` lists was found unavailable, so there is no ",
         "sample to estimate from.", call. = FALSE)
  }
  if (need > 0L) {
    warn_drawn(paste0(
      "The ", r, " reserve(s) of `sample` ran out with ", need, " unavailable ",
      "unit(s) of the main sample not replaced, so the sample holds ", m,
      " units rather than ", design$n, ". It is weighted as a sample of ", m,
      " from the units not found unavailable."), "drawn_reserves_short")
  }

  resolved <- design
  resolved$n <- as.integer(m)
  resolved$reserve <- 0L
  resolved$unavailable <- sort(rows[gone %in% TRUE])
  pi_all <- exact_inclusion(resolved, state$pop)
  keep <- which(used)
  rebuild_sample(reindex(sample, keep), resolved, state$pop, rows[keep],
                 pi_all, reserve = place[keep], compact = state$compact)
}

#' What the visits to a list found, as a logical vector over its rows
#' @noRd
visit_findings <- function(sample, unavailable) {
  if (is.character(unavailable) && length(unavailable) == 1L &&
      !is.na(unavailable)) {
    if (!unavailable %in% names(sample)) {
      stop("`sample` has no column `", unavailable, "`.", call. = FALSE)
    }
    unavailable <- sample[[unavailable]]
  }
  if (!is.logical(unavailable) || !is.null(dim(unavailable)) ||
      length(unavailable) != nrow(sample)) {
    stop("`unavailable` must be a logical vector with one value for each of ",
         "the ", nrow(sample), " rows of `sample`, or the name of such a ",
         "column, not ", format_bad(unavailable), ". Mark a unit TRUE when it ",
         "was found unavailable, FALSE when it was available and NA when it ",
         "was not visited; `sample$id %in% ids` turns a list of ids into ",
         "such a vector.", call. = FALSE)
  }
  as.vector(unavailable)
}

#' Drop the frame columns a sample's design does not need
#'
#' A weighted sample holds its own copy of the whole frame, every column of
#' it, so that [ht_total()], [ht_mean()], [as_svydesign()] and
#' [sample_summary()] can work from the sample alone. That copy is kept in
#' memory and written to every file the sample is saved to, which costs as
#' much as the frame each time, and hands the frame, with any names,
#' addresses or other personal data in it, to whoever receives the file.
#' `compact_sample()` replaces the copy with one that holds only the columns
#' the design reads (strata, clusters, size measures, `order_by`, time
#' stamps, coordinates and the like) and the columns named in `keep`.
#'
#' Estimates, standard errors, [as_svydesign()] and [sample_summary()] are
#' unchanged, with two exceptions that read other frame columns. A `by`
#' domain's size in the frame, which `ht_mean(estimator = "ht")` with `by`
#' divides by, needs the `by` columns in the stored frame; without them that
#' estimator is refused with a message naming `keep`. And a stratum with a
#' single sampled row leaves a domain's variance unmeasured unless the
#' stored frame's `by` columns show the stratum holds none of the domain (see
#' "Domains" in [ht_total()]). Keep the `by` columns of the estimates in
#' `keep`, and the columns [frame_totals()] is asked for. Whether a column of
#' the frame left out holds only 0 and 1, which makes [ht_mean()] treat it as
#' a proportion, is recorded.
#'
#' The sample's own columns are not changed: drop columns of the sample with
#' `s$col <- NULL`. The inclusion probabilities of every frame row are kept,
#' one number per row.
#'
#' @param sample A data frame returned by [draw()] with `weights = TRUE`, or
#'   by [as_drawn()].
#' @param keep Names of other columns of the frame to keep, such as the `by`
#'   columns of the estimates.
#'
#' @return `sample`, with a smaller copy of the frame.
#'
#' @examples
#' pop <- data.frame(id = 1:5000, site = rep(c("a", "b"), each = 2500),
#'                   note = strrep("x", 40), spend = (1:5000) %% 97)
#' s <- draw(pop, design_stratified("site", n = 100), seed = 1, weights = TRUE)
#' small <- compact_sample(s)
#'
#' # The same estimate from a fraction of the bytes
#' c(full = length(serialize(s, NULL)),
#'   compact = length(serialize(small, NULL)))
#' all.equal(ht_total(s, "spend")$se, ht_total(small, "spend")$se)
#'
#' @seealso [draw()], [as_drawn()]
#' @export
compact_sample <- function(sample, keep = NULL) {
  state <- drawn_state(sample, reserves = TRUE)
  pop <- state$pop
  if (!is.null(keep)) {
    keep <- check_columns(keep, "keep")
    unknown <- setdiff(keep, names(pop))
    if (length(unknown)) {
      gone <- intersect(unknown, state$compact$dropped)
      stop("`keep` names ", quote_names(unknown), ", which ",
           if (length(unknown) > 1L) "are not columns" else "is not a column",
           " of the frame stored with `sample`. ",
           if (length(gone)) {
             paste0("An earlier compact_sample() dropped ",
                    quote_names(gone), "; compact the sample draw() returned ",
                    "instead.")
           } else {
             "Name columns of the frame the sample was drawn from."
           }, call. = FALSE)
    }
  }
  cols <- names(pop)[names(pop) %in% c(design_columns(state$design), keep)]
  small <- frame_columns(pop, cols)
  dropped <- setdiff(names(pop), names(small))
  binary <- dropped[vapply(dropped, function(nm) binary_column(pop[[nm]]),
                           logical(1))]
  out <- sample
  attr(out, "drawn_population") <- small
  attr(out, "drawn_compact") <- list(
    dropped = union(state$compact$dropped, dropped),
    binary = union(state$compact$binary, binary)
  )
  table_ready(out)
}

#' The frame columns a design reads
#'
#' Strata, the column of standard deviations a Neyman allocation reads,
#' clusters, weights and size measures, the columns a spread design spreads
#' across, time stamps, `order_by` (of a systematic design, or of the walk
#' within each stratum) and coordinates. A certainty design adds `above` to
#' those of its `rest`.
#'
#' @noRd
design_columns <- function(design) {
  if (inherits(design, "drawn_design_certainty")) {
    return(unique(c(design[["above"]], design_columns(design[["rest"]]))))
  }
  allocation_by <- design[["allocation_by"]]
  within <- design[["within"]]
  unique(c(design[["strata"]],
           if (is.character(allocation_by)) allocation_by,
           design[["clusters"]], design[["weights"]], design[["across"]],
           design[["size"]], design[["time"]], design[["order_by"]],
           design[["coords"]],
           if (is_design(within)) design_columns(within)))
}

#' Columns `cols` of a frame, as a plain data frame
#'
#' An `sf` frame whose geometry is one of `cols` keeps its class, since a
#' design reads coordinates and their CRS from it. Any other frame becomes a
#' plain data frame, which estimates read the same way, and an `sf` frame
#' whose geometry the design does not read loses it.
#'
#' @noRd
frame_columns <- function(pop, cols) {
  geometry <- attr(pop, "sf_column", exact = TRUE)
  if (!is.null(geometry) && geometry %in% cols) {
    return(pop[, cols, drop = FALSE])
  }
  out <- unclass(pop)[cols]
  structure(out, names = cols, row.names = .set_row_names(nrow(pop)),
            class = "data.frame")
}

#' Does a column hold a proportion: logical, or numbers that are all 0 or 1?
#' @noRd
binary_column <- function(x) {
  if (is.logical(x)) return(TRUE)
  if (!is.numeric(x)) return(FALSE)
  x <- x[!is.na(x)]
  length(x) > 0L && all(x == 0 | x == 1)
}
