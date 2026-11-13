#' Draw a sample
#'
#' Applies a sampling design to data.
#'
#' @section What you get back:
#' A data frame with the same class and the same columns, in the same order, as
#' `data`. `weights = TRUE` prepends three columns (see "Design weights"), and
#' two designs qualify it further:
#'
#' * [design_bootstrap()] prepends a `.replicate` column identifying which
#'   replicate each row belongs to, so all replicates come back in one frame.
#'   Split them with `split(out, out$.replicate)` if you need a list.
#' * [design_reservoir()] returns a list when `data` is a stream rather than a
#'   data frame, because there is nothing to make a data frame out of.
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
#' * The frame is copied when the sample is drawn. Later edits to `data`,
#'   including data.table's in-place `:=`, `setorder()` and `setkey()`, do not
#'   reach the sample, and saving the sample with [saveRDS()] saves the copy
#'   with it.
#' * Reordering the sample by `[`, `dplyr::arrange()` or
#'   `data.table::setorder()` is safe, because each row finds its frame row
#'   through `.row`.
#' * Estimate from the whole sample and use `by` for groups. Row subsets
#'   (`head()`, `s[cond, ]`, `dplyr::filter()`), [rbind()] of samples and joins
#'   that repeat rows keep the attributes but no longer match the rows that
#'   were drawn, and the functions above refuse them. [merge()], column
#'   selection (`s[, cols]`), [subset()] and [transform()] drop the attributes
#'   altogether. Add fieldwork results with `s$col <- value`.
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
#' `seed` is saved, applied, and unwound: `.Random.seed` is restored on exit, so
#' drawing a sample inside a simulation does not shift the simulation's own
#' random number stream. Leave it `NULL` to draw from the current stream.
#'
#' @param data A data frame. [design_reservoir()] also accepts a list, a
#'   connection, or a zero-argument generator function.
#' @param design A design built by one of the [design_simple()] family.
#' @param seed Optional seed, applied only for this draw.
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
  if (!is_design(design)) {
    stop("`design` must come from one of the design_*() constructors, not a ",
         class(design)[1], ". See ?designs.", call. = FALSE)
  }
  check_flag(weights, "weights")

  if (!isTRUE(weights)) {
    return(table_ready(with_seed(seed, draw_design(design, data))))
  }

  if (!is.data.frame(data)) {
    stop("`weights = TRUE` needs a data frame, not a ", class(data)[1], ".",
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
    stop("`weights = TRUE` is not available for ",
         if (inherits(design, "drawn_design_certainty")) {
           "a certainty design over a with-replacement `rest`"
         } else {
           "a with-replacement design"
         }, ".\n",
         "`.prob` would be the chance of being selected *at least once*, but ",
         "the sample holds\nduplicates, so summing `y * .weight` over it ",
         "overstates the total by a factor of\n(n / N) / (1 - (1 - 1 / N)^n), ",
         "about 5% at n = N / 10, 15% at n = 0.3 * N, 58% at\nn = N and more ",
         "beyond, within each stratum of a stratified design. Use\n",
         "replace = FALSE, or estimate with the Hansen-Hurwitz form, ",
         "N / n * sum(y), summed\nover strata as N_h / n_h * sum(y_h).",
         call. = FALSE)
  }

  # Probabilities first, so a design with no closed form fails before any
  # sampling happens rather than after.
  p_all <- exact_inclusion(design, data)

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

  out <- with_seed(seed, draw_design(design, tagged))
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
  table_ready(out)
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

#' Design-specific draw method
#'
#' @param design A design object.
#' @param data The data to sample from.
#' @return The sampled rows.
#' @keywords internal
#' @export
draw_design <- function(design, data) {
  UseMethod("draw_design")
}

#' @export
draw_design.default <- function(design, data) {
  stop("No draw method for design type \"", design_type(design), "\".",
       call. = FALSE)
}
