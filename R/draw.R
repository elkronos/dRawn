#' Draw a sample
#'
#' Applies a sampling design to data.
#'
#' @section What you get back:
#' A data frame with the same class and the same columns, in the same order, as
#' `data`. Two designs qualify that:
#'
#' * [design_bootstrap()] prepends a `.replicate` column identifying which
#'   replicate each row belongs to, so all replicates come back in one frame.
#'   Split them with `split(out, out$.replicate)` if you need a list.
#' * [design_reservoir()] returns a list when `data` is a stream rather than a
#'   data frame, because there is nothing to make a data frame out of.
#'
#' @section Design weights:
#' `weights = TRUE` prepends two columns: `.prob`, the probability that the row
#' was included in the sample, and `.weight`, its reciprocal — the number of
#' population units the row stands for. Those are what an unbiased estimate
#' needs: a Horvitz-Thompson total is `sum(y * .weight)`.
#'
#' They come from the design applied to the population, not from the drawn
#' sample, and some designs have no closed form for them. [inclusion_prob()]
#' documents which, why, and what to do instead.
#'
#' `weights = TRUE` is not available for a design that samples **with
#' replacement**. `.prob` there is the probability of being selected *at least
#' once*, but the sample holds duplicates, so `sum(y * .weight)` over it counts
#' each duplicate at the distinct-unit weight and comes out around 15% high.
#' Use `replace = FALSE`, or the Hansen-Hurwitz form `N / n * sum(y)`.
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
#' @param weights Attach `.prob` and `.weight` columns. See "Design weights".
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
    return(with_seed(seed, draw_design(design, data)))
  }

  if (!is.data.frame(data)) {
    stop("`weights = TRUE` needs a data frame, not a ", class(data)[1], ".",
         call. = FALSE)
  }
  clash <- intersect(c(".prob", ".weight"), names(data))
  if (length(clash) > 0L) {
    stop("`data` already has ", paste0("`", clash, "`", collapse = " and "),
         "; weights = TRUE would overwrite ",
         if (length(clash) == 1L) "it" else "them", ".", call. = FALSE)
  }
  if (anyDuplicated(names(data))) {
    dup <- unique(names(data)[duplicated(names(data))])
    stop("`data` has duplicated column name(s): ",
         format_bad(dup), ". Rename them before drawing with weights, which ",
         "cannot tell them apart.", call. = FALSE)
  }
  if (isTRUE(design$replace)) {
    stop("`weights = TRUE` is not available for a with-replacement design.\n",
         "`.prob` would be the chance of being selected *at least once*, but ",
         "the sample holds\nduplicates, so summing `y * .weight` over it ",
         "double-counts and comes out roughly 15%\nhigh. Use replace = FALSE, ",
         "or estimate with the Hansen-Hurwitz form, N / n * sum(y).",
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
  pos <- out[[key]]
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
  out <- out[, c(".prob", ".weight", setdiff(names(out), c(".prob", ".weight"))),
             drop = FALSE]

  # Carry the design and the source rows so ht_total() can recover the joint
  # inclusion probabilities without being handed the population again.
  attr(out, "drawn_design") <- design
  attr(out, "drawn_rows") <- pos
  attr(out, "drawn_population") <- data
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
