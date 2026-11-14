#' Draw a sample
#'
#' Applies a sampling design to data.
#'
#' @section What you get back:
#' A data frame with the same class and the same columns, in the same order, as
#' `data`. `weights = TRUE` prepends three columns (see "Design weights"), and
#' some designs qualify it further:
#'
#' * [design_bootstrap()] prepends a `.replicate` column identifying which
#'   replicate each row belongs to, so all replicates come back in one frame.
#'   Split them with `split(out, out$.replicate)` if you need a list. With
#'   `output = "indices"` or `"weights"` it returns a matrix of row numbers or
#'   counts instead, with a column per replicate.
#' * [design_reservoir()] returns a list when `data` is a list or a stream
#'   rather than a data frame, because there is nothing to make a data frame
#'   out of. With `batches = TRUE` it returns the rows sampled from a
#'   function's data frames as one data frame.
#' * [design_simple()] and [design_spread()] with `reserve` return a list:
#'   the main sample, then the reserves, with a leading integer `.reserve`
#'   column that is 0 for the main sample and numbers the reserves 1, 2, ...
#'   in the order they are used. With `weights = TRUE` it follows `.row`,
#'   and `.prob` is the main sample's probability for every row.
#'   [use_reserves()] turns the list into the sample to estimate from.
#'
#' Row names of the frame's own, such as the car names of `mtcars`, stay with
#' their rows. A row drawn more than once, by a design with `replace = TRUE`
#' or by [design_bootstrap()], gets names made unique by [make.unique()]:
#' `"Mazda RX4"`, then `"Mazda RX4.1"`. The automatic row names 1, 2, ...
#' number the sample's rows afresh, and `weights = TRUE` records each row's
#' place in `data` in `.row`. A tibble has no row names to keep.
#'
#' @section Design weights:
#' `weights = TRUE` prepends three columns: `.prob`, the probability that the
#' row was included in the sample; `.weight`, its reciprocal, the number of
#' population units the row stands for; and `.row`, the row of `data` it came
#' from. `.prob` and `.weight` are what an unbiased estimate needs: a
#' Horvitz-Thompson total is `sum(y * .weight)`.
#'
#' They come from the design applied to the population, not from the drawn
#' sample, and some designs have no closed form for them. [inclusion_prob()]
#' documents which, why, and what to do instead.
#'
#' Two kinds of design are weighted conditional on the clusters the draw
#' selected: [design_cluster()] with `balanced = TRUE`, and
#' [design_multistage()] without `size` when the rows taken from a cluster
#' depend on which other clusters were selected. With `a` of `A` clusters
#' selected and `m_i` of the `N_i` rows of cluster `i` taken, `.prob` is
#' `(a / A) * m_i / N_i`, the chance of the cluster being selected times the
#' chance of the row being drawn given these clusters. The total estimated
#' from it is unbiased given the clusters selected, and so over every
#' selection. [sample_summary()] says when a sample is weighted this way.
#'
#' Keep `.prob` as it is: [ht_total()], [ht_mean()] and [as_svydesign()] refuse
#' a sample whose `.prob` was changed. Weights adjusted after the draw, for
#' nonresponse say, go in `.weight`, and [as_svydesign()] is the way to
#' estimate from them; [ht_total()] and [ht_mean()] keep to `.prob` and warn
#' when `.weight` no longer matches it.
#'
#' The sample also carries its design, the frame rows it came from and a copy
#' of the frame as attributes, which is what lets [ht_total()], [ht_mean()],
#' [sample_summary()] and [as_svydesign()] work from the sample alone.
#'
#' * The frame is copied when the sample is drawn, every column of it. Later
#'   edits to `data`, including data.table's in-place `:=`, `setorder()` and
#'   `setkey()`, do not reach the sample. Each sample holds its own copy, in
#'   memory and in every file it is saved to with [saveRDS()], so a list of
#'   samples from a large frame holds the frame once per sample, and a
#'   sample sent to someone sends them the frame, with any personal data in
#'   it. [compact_sample()] keeps only the columns the design reads.
#' * Reordering the sample by `[`, `dplyr::arrange()` or
#'   `data.table::setorder()` is safe, because each row finds its frame row
#'   through `.row`.
#' * Estimate from the whole sample and use `by` for groups. Row subsets
#'   (`head()`, `s[cond, ]`, `dplyr::filter()`), [rbind()] of samples and joins
#'   that repeat rows keep the attributes but no longer match the rows that
#'   were drawn, and the functions above refuse them. [merge()], column
#'   selection (`s[, cols]`), [subset()] and [transform()] drop the attributes
#'   altogether. Add fieldwork results with [add_fieldwork()], which matches
#'   them on an id, or with `s$col <- value`.
#' * A CSV or spreadsheet file and a database table keep the rows and the
#'   `.row` column but not the attributes. [as_drawn()] rebuilds the sample
#'   from them, the design and the frame.
#'
#' `weights = TRUE` is not available for a design that samples **with
#' replacement**, or a [design_certainty()] whose `rest` does. `.prob` there
#' is the probability of being selected *at least once*, but the sample holds
#' duplicates, so `sum(y * .weight)` over it counts each duplicate at the
#' distinct-unit weight. The total comes out high by a factor of
#' `(n / N) / (1 - (1 - 1 / N)^n)`, which is about 5% at `n = N / 10`, 15% at
#' `n = 0.3 * N`, 58% at `n = N` and more beyond, applied within each stratum
#' of a stratified design. Use `replace = FALSE`, or the Hansen-Hurwitz form
#' `N / n * sum(y)`, summed over strata as `N_h / n_h * sum(y_h)` for a
#' stratified design.
#'
#' @section Seeding:
#' With a `seed`, `draw()` saves the caller's random number state, seeds,
#' draws, and restores the state on exit, so drawing a sample inside a
#' simulation does not shift the simulation's own random number stream. With
#' `seed = NULL` it draws from the caller's stream and advances it, as any
#' function that draws random numbers does.
#'
#' A seed is applied with R's default generator,
#' `set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
#' sample.kind = "Rejection")`, whatever generator the session has selected,
#' and the session's generator is restored on exit with its state. A seed
#' therefore draws the same sample after `RNGkind("L'Ecuyer-CMRG")`, which
#' parallel code often selects, or under `sample.kind = "Rounding"` as it
#' does in a session left at R's defaults. It must be a whole number from
#' `-.Machine$integer.max` to `.Machine$integer.max`. The same holds for the
#' `seed` of [inclusion_prob()], [joint_prob()] and [plot()][plot.drawn_design].
#' A design selects rows by their position in `data`, so a seed draws the
#' same rows again only from a frame with the same rows in the same order.
#'
#' With `weights = TRUE` the sample records how it was drawn in its
#' `drawn_meta` attribute: `seed` (`NULL` without one), `rng_kind` (the
#' generator used, as [RNGkind()] gives it), `version` (of drawn), `time`,
#' `frame_rows` (the rows of `data`) and `made_by` (`"draw()"`, or
#' `"as_drawn()"` for a sample rebuilt from a file). [sample_summary()]
#' prints it.
#'
#' @param data A data frame. [design_reservoir()] also accepts a list, a
#'   connection, or a zero-argument generator function.
#' @param design A design built by one of the [design_simple()] family.
#' @param seed Optional seed, applied only for this draw. See "Seeding".
#' @param weights Attach `.prob`, `.weight` and `.row` columns. See "Design
#'   weights".
#'
#' @return The sampled rows. See "What you get back".
#'
#' @examples
#' df <- data.frame(id = 1:100, site = rep(letters[1:4], each = 25))
#'
#' draw(df, design_simple(n = 10), seed = 1)
#'
#' # The same design, reused
#' by_site <- design_stratified(strata = "site", n = 12)
#' table(draw(df, by_site, seed = 1)$site)
#' table(draw(df[1:60, ], by_site, seed = 1)$site)
#'
#' # With design weights, ready for estimation
#' df$spend <- seq_len(100)
#' s <- draw(df, by_site, seed = 1, weights = TRUE)
#' sum(s$spend * s$.weight)   # estimates sum(df$spend) = 5050
#'
#' @seealso [designs] for the full list of constructors, [inclusion_prob()] for
#'   the probabilities themselves, and [ht_total()] to estimate a population
#'   total with a standard error.
#' @export
draw <- function(data, design, seed = NULL, weights = FALSE) {
  check_design(design, data, "draw")
  check_seed(seed)
  check_flag(weights, "weights")
  warn_region(design, data)

  if (!isTRUE(weights)) {
    out <- with_seed(seed, draw_design(design, data))
    return(table_ready(drop_selection(out)))
  }

  if (!is.data.frame(data)) {
    stop("`weights = TRUE` needs a data frame, not ", a_class(data), ".",
         call. = FALSE)
  }
  added <- c(".prob", ".weight", ".row")
  clash <- intersect(added, names(data))
  if (length(clash) > 0L) {
    shown <- paste0("`", clash, "`")
    stop("`data` already has ",
         if (length(shown) > 1L) {
           paste(paste(shown[-length(shown)], collapse = ", "), "and",
                 shown[length(shown)])
         } else shown,
         "; weights = TRUE would overwrite ",
         if (length(clash) == 1L) "it" else "them", ".", call. = FALSE)
  }
  if (anyDuplicated(names(data))) {
    dup <- unique(names(data)[duplicated(names(data))])
    stop("`data` has duplicated column name(s): ",
         format_bad(dup), ". Rename them before drawing with weights, which ",
         "cannot tell them apart.", call. = FALSE)
  }
  if (samples_with_replacement(design)) {
    stop("`weights = TRUE` is not available for ", replacement_design(design),
         ".\n", replacement_problem(design, "would be"), call. = FALSE)
  }

  # Probabilities first, so a design with no closed form fails before any
  # sampling happens rather than after. A design weighted conditional on its
  # clusters has probabilities only once they are selected, so it gets its
  # checks here and its probabilities after the draw.
  conditional <- conditional_design(design, data)
  if (conditional) {
    conditional_setup(design, data)
  } else {
    p_all <- exact_inclusion(design, data)
    if (inherits(design, "drawn_design_certainty")) {
      check_rest_reach(design, data, p_all)
    }
  }

  # Tag the population, draw from the tagged copy, then read the source row off
  # each drawn row. Matching on values afterwards would misattribute duplicate
  # rows and every with-replacement draw.
  key <- ".drawn_row_id"
  if (key %in% names(data)) {
    stop("`data` already has a column called `", key,
         "`, which weights = TRUE needs. Rename it.", call. = FALSE)
  }
  tagged <- data
  tagged[[key]] <- seq_len(nrow(data))

  # Read before the draw, since with_seed() changes the generator while the
  # draw runs.
  meta <- draw_meta(seed, nrow(data))
  out <- with_seed(seed, draw_design(design, tagged))
  selection <- attr(out, "drawn_selection", exact = TRUE)
  out <- drop_selection(out)
  if (conditional) p_all <- conditional_inclusion(design, data, selection)
  pos <- as.integer(out[[key]])
  out[[key]] <- NULL

  p <- p_all[pos]
  w <- rep(NA_real_, length(p))
  ok <- !is.na(p) & p > 0
  w[ok] <- 1 / p[ok]

  # cbind() on a data.frame drops a tibble's class, and the contract is that
  # what comes back matches what went in. Build the columns in place instead.
  out <- reindex(out, seq_len(nrow(out)))
  out[[".prob"]] <- p
  out[[".weight"]] <- w
  out[[".row"]] <- pos
  out <- out[, c(added, setdiff(names(out), added)), drop = FALSE]

  # Carry the design, the source rows, the frame and its probabilities so
  # ht_total() can recover the joint inclusion probabilities without being
  # handed the population again. The frame is a copy, because a data.table
  # edited in place after the draw would otherwise change the sample's
  # estimates.
  attr(out, "drawn_design") <- design
  # `.row` holds `pos` itself, so the attribute gets its own copy; an
  # in-place data.table::setorder() of the sample would reorder both.
  attr(out, "drawn_rows") <- copy_column(pos)
  attr(out, "drawn_population") <- snapshot_frame(data)
  attr(out, "drawn_pi") <- p_all
  # The clusters selected, in the order selected, from which drawn_state()
  # works out the conditional probabilities again.
  if (conditional) attr(out, "drawn_selection") <- selection
  attr(out, "drawn_meta") <- meta
  table_ready(out)
}

#' How a sample was drawn, for the record it carries
#'
#' The seed (`NULL` without one), the generator the draw used (the pinned
#' one with a seed, else the session's), the version of drawn, the time, the
#' number of frame rows, and the function that made the sample. A sample
#' rebuilt by as_drawn() records neither seed nor generator, which the file
#' it was read from does not hold.
#'
#' @noRd
draw_meta <- function(seed, frame_rows, made_by = "draw()") {
  drawn_here <- identical(made_by, "draw()")
  list(
    seed = seed,
    rng_kind = if (!drawn_here) {
      NULL
    } else if (is.null(seed)) {
      RNGkind()
    } else {
      seed_kind()
    },
    version = as.character(utils::packageVersion("drawn")),
    time = Sys.time(),
    frame_rows = as.integer(frame_rows),
    made_by = made_by
  )
}

#' Remove the clusters a cluster or multistage draw records
#'
#' Their draw_design() methods attach the codes of the clusters they selected
#' as `drawn_selection`. Only a sample weighted conditional on its clusters
#' keeps them, set by draw() itself.
#'
#' @noRd
drop_selection <- function(out) {
  if (!is.null(attr(out, "drawn_selection", exact = TRUE))) {
    attr(out, "drawn_selection") <- NULL
  }
  out
}

#' Does a design put the same row in a sample more than once?
#'
#' A certainty design has no `replace` of its own; its `rest` decides.
#'
#' @noRd
samples_with_replacement <- function(design) {
  isTRUE(design$replace) ||
    (inherits(design, "drawn_design_certainty") &&
       samples_with_replacement(design$rest))
}

#' A with-replacement design, named for a message
#' @noRd
replacement_design <- function(design) {
  if (inherits(design, "drawn_design_certainty")) {
    "a certainty design over a with-replacement `rest`"
  } else {
    "a with-replacement design"
  }
}

#' Why a with-replacement design gets no weights, and what to do instead
#'
#' draw() gives this reason when it refuses `weights = TRUE`, with `tense`
#' "would be", and drawn_state() when it refuses a sample that carries such
#' a design, with "is". A certainty design takes `replace` in its `rest`.
#'
#' @noRd
replacement_problem <- function(design, tense = "is") {
  paste0(
    "`.prob` ", tense, " the chance of being selected *at least once*, but ",
    "the sample holds\nduplicates, so summing `y * .weight` over it ",
    "overstates the total by a factor of\n(n / N) / (1 - (1 - 1 / N)^n), ",
    "about 5% at n = N / 10, 15% at n = 0.3 * N, 58% at\nn = N and more ",
    "beyond, within each stratum of a stratified design. Use\n",
    if (inherits(design, "drawn_design_certainty")) {
      paste0("replace = FALSE in `rest`, or estimate with the Hansen-Hurwitz ",
             "form,\nN / n * sum(y), summed over strata as ",
             "N_h / n_h * sum(y_h).")
    } else {
      paste0("replace = FALSE, or estimate with the Hansen-Hurwitz form, ",
             "N / n * sum(y), summed\nover strata as N_h / n_h * sum(y_h).")
    })
}

#' Let `:=` work on a data.table built by base subsetting
#'
#' Base `[` and `[[<-` return a data.table without the spare column slots
#' data.table allocates, so the first `:=` on it warns about an invalid
#' self-reference and copies. setalloccol() restores the slots. It runs last,
#' because `attr<-` would copy the table again.
#'
#' @noRd
table_ready <- function(out) {
  if (inherits(out, "data.table") &&
      requireNamespace("data.table", quietly = TRUE)) {
    out <- data.table::setalloccol(out)
  }
  out
}

#' Draw method for a design class
#'
#' The generic [draw()] calls to make a sample from a design. Each design in
#' the package has a method, and a design of your own can have one too.
#'
#' @section Writing a design of your own:
#' A design is a list with a class vector ending in `"drawn_design"`, so that
#' [is_design()] accepts it, such as
#' `structure(list(n = 5), class = c("first_rows", "drawn_design"))`. Its
#' method takes `design` and `data` and returns rows of `data`, with every
#' column of `data`, as `data[idx, , drop = FALSE]` gives them. It draws its
#' random numbers from the current stream; [draw()] applies the seed around
#' it. A package registers the method with `@exportS3Method
#' drawn::draw_design`, and a script with [.S3method()].
#'
#' That is enough for [draw()] without weights, and for [inclusion_prob()]
#' and [joint_prob()] with `simulate = TRUE`, which count how often each row
#' is drawn and find the rows through a column they add to `data`. The exact
#' probabilities, `draw(weights = TRUE)` and the estimators need the
#' design's inclusion probabilities, which the package computes for its own
#' designs only, so they refuse a design of your own.
#'
#' @param design A design object.
#' @param data The data to sample from.
#' @return The sampled rows. The cluster and multistage methods attach the
#'   codes of the clusters they selected as the attribute `drawn_selection`,
#'   which [draw()] removes or, for a sample weighted conditional on its
#'   clusters, keeps.
#' @examples
#' first_rows <- structure(list(n = 3), class = c("first_rows", "drawn_design"))
#' .S3method("draw_design", "first_rows", function(design, data) {
#'   data[seq_len(design$n), , drop = FALSE]
#' })
#'
#' df <- data.frame(id = 1:10)
#' draw(df, first_rows)
#' inclusion_prob(df, first_rows, simulate = TRUE, R = 20, seed = 1)
#' @export
draw_design <- function(design, data) {
  UseMethod("draw_design")
}

#' @export
draw_design.default <- function(design, data) {
  stop("No draw method for design type \"", design_type(design), "\".",
       call. = FALSE)
}
