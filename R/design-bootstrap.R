#' Bootstrap sampling
#'
#' Generates bootstrap replicates. `"simple"` resamples rows independently;
#' `"block"` is a moving-block bootstrap, which preserves the serial dependence
#' in ordered data by concatenating `ceiling(n / block_length)` independently
#' chosen blocks per replicate.
#'
#' [draw()] returns all replicates in one data frame with a leading `.replicate`
#' column. Split them with `split(out, out$.replicate)`.
#'
#' @param n_replicates Number of replicates to generate.
#' @param n Rows per replicate. `NULL` uses `nrow(data)`, the standard
#'   nonparametric bootstrap.
#' @param method `"simple"` or `"block"`.
#' @param block_length Block length for `method = "block"`. `NULL` uses
#'   `floor(nrow(data)^(1/3))`, a common rule of thumb.
#'
#' @return A design object, for use with [draw()].
#'
#' @examples
#' df <- data.frame(id = 1:100, value = (1:100) / 10)
#'
#' reps <- draw(df, design_bootstrap(n_replicates = 5, n = 20), seed = 1)
#' table(reps$.replicate)
#'
#' # A statistic per replicate
#' vapply(split(reps, reps$.replicate), function(r) mean(r$value), numeric(1))
#'
#' @family designs
#' @seealso [draw()]
#' @export
design_bootstrap <- function(n_replicates = 1000L, n = NULL,
                             method = c("simple", "block"),
                             block_length = NULL) {
  new_design("bootstrap", list(
    n_replicates = check_count(n_replicates, "n_replicates", allow_zero = FALSE),
    n = if (is.null(n)) NULL else check_count(n, "n"),
    method = match.arg(method),
    block_length = if (is.null(block_length)) NULL else
      check_count(block_length, "block_length", allow_zero = FALSE)
  ))
}

#' @export
draw_design.drawn_design_bootstrap <- function(design, data) {
  validate_data(data)

  n_rows <- nrow(data)
  n <- design$n %||% n_rows

  if (design$method == "block") {
    block_length <- design$block_length %||% max(1L, floor(n_rows^(1 / 3)))
    if (block_length > n_rows) {
      stop("`block_length` (", block_length, ") exceeds the number of rows (",
           n_rows, ").", call. = FALSE)
    }
  }

  parts <- lapply(seq_len(design$n_replicates), function(i) {
    idx <- if (design$method == "simple") {
      sample.int(n_rows, size = n, replace = TRUE)
    } else {
      n_blocks <- ceiling(n / block_length)
      starts <- sample.int(n_rows - block_length + 1L, n_blocks, replace = TRUE)
      # Several blocks per replicate, so it is not confined to one run.
      unlist(lapply(starts, function(s) s:(s + block_length - 1L)),
             use.names = FALSE)[seq_len(n)]
    }
    out <- reindex(data, idx)
    cbind(.replicate = rep(i, n), out)
  })

  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_bootstrap <- function(design, data) {
  stop("A bootstrap resamples the sample; it is not a probability sample of a ",
       "finite population, so it has no inclusion probability.\n",
       "If you want the expected number of times each row appears, that is ",
       "n_replicates * n / nrow(data).", call. = FALSE)
}
