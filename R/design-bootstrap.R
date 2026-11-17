#' Bootstrap sampling
#'
#' Generates bootstrap replicates. `"simple"` resamples rows independently.
#' `"block"` is the moving-block bootstrap (Künsch 1989), which keeps the
#' serial dependence in ordered data by concatenating `ceiling(n /
#' block_length)` runs of `block_length` consecutive rows per replicate, each
#' starting at a row chosen at random. A run must fit inside the frame, so
#' fewer runs cover the first and last `block_length - 1` rows, and the first
#' row and the last turn up only `1 / block_length` as often as a row in the
#' middle. `"circular"` (Politis and Romano 1992) joins the last row to the
#' first, so a run can wrap around the end and every row turns up equally
#' often on average.
#'
#' By default [draw()] returns all replicates in one data frame with a leading
#' `.replicate` column. Split them with `split(out, out$.replicate)`. That
#' frame holds `n_replicates * n` rows, a copy of the data for each replicate,
#' and [draw()] warns, with class `drawn_bootstrap_size`, when it would exceed
#' 10 million rows. `output = "indices"` and `output = "weights"` return the
#' same replicates as a matrix instead.
#'
#' @param n_replicates Number of replicates to generate.
#' @param n Rows per replicate. `NULL` uses `nrow(data)`, the standard
#'   nonparametric bootstrap. `n = 0` gives replicates with no rows. The
#'   frame then has no rows but keeps its `.replicate` column, so `split()`
#'   returns no replicates at all, while
#'   `split(out, factor(out$.replicate, seq_len(n_replicates)))` returns
#'   `n_replicates` empty ones. A plan from [plan_size()] gives its `n`.
#' @param method `"simple"`, `"block"` or `"circular"`.
#' @param block_length Run length for `method = "block"` and `"circular"`, and
#'   an error for `method = "simple"`, where it would have no effect. `NULL`
#'   uses `floor(nrow(data)^(1/3))`, a common rule of thumb.
#' @param output What [draw()] returns.
#'   * `"rows"`, the default, returns the resampled rows of every replicate in
#'     one data frame, with a leading `.replicate` column.
#'   * `"indices"` returns an integer matrix with `n` rows and a column per
#'     replicate, holding the rows of `data` that each replicate drew, in
#'     order. Row `i` of replicate `r` is `data[out[i, r], ]`.
#'   * `"weights"` returns an integer matrix with a row for each row of `data`
#'     and a column per replicate, holding the number of times each replicate
#'     drew that row. With the default `n`, these are replicate weights for
#'     `survey::svrepdesign(type = "bootstrap", combined.weights = TRUE)`.
#'
#'   For one seed, all three hold the same replicates.
#'
#' @section Not a survey design:
#' The bootstrap resamples the data you already have to approximate the
#' sampling distribution of a statistic (Efron 1979). It is not a probability
#' sample of a finite population, so it has no inclusion probabilities and
#' cannot be passed to [ht_total()]. For the design-based variance of a total
#' or mean, use the estimators there.
#'
#' Resampling rows treats them as independent draws from one distribution.
#' For such data and a statistic those estimators do not cover, such as a
#' median, the bootstrap is the standard tool. A sample drawn with strata or
#' clusters is not such data, and resampling its rows one by one ignores the
#' strata and treats the rows of a cluster as independent. Resample the
#' clusters within strata instead, for example with
#' `survey::as.svrepdesign(as_svydesign(s), type = "subbootstrap")`. Runs of
#' rows keep short-range dependence in ordered data (Künsch 1989).
#'
#' @return A design object, for use with [draw()].
#'
#' @references
#' Efron, B. (1979). Bootstrap methods: another look at the jackknife.
#' *Annals of Statistics*, 7, 1–26.
#'
#' Künsch, H. R. (1989). The jackknife and the bootstrap for general
#' stationary observations. *Annals of Statistics*, 17, 1217–1241.
#'
#' Politis, D. N. and Romano, J. P. (1992). A circular block-resampling
#' procedure for stationary data. In R. LePage and L. Billard (eds.),
#' *Exploring the Limits of Bootstrap*, 263–270. Wiley.
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
#' # The same replicates as row numbers, without copying the data
#' idx <- draw(df, design_bootstrap(n_replicates = 5, n = 20,
#'                                  output = "indices"), seed = 1)
#' apply(idx, 2, function(i) mean(df$value[i]))
#'
#' # Runs of 10 rows that may wrap from the last row to the first
#' circ <- draw(df, design_bootstrap(n_replicates = 200, method = "circular",
#'                                   block_length = 10, output = "weights"),
#'              seed = 1)
#' range(rowMeans(circ))
#'
#' @family designs
#' @seealso [draw()]
#' @export
design_bootstrap <- function(n_replicates = 1000L, n = NULL,
                             method = c("simple", "block", "circular"),
                             block_length = NULL,
                             output = c("rows", "indices", "weights")) {
  method <- check_choice(method, c("simple", "block", "circular"), "method")
  output <- check_choice(output, c("rows", "indices", "weights"), "output")
  if (method == "simple" && !is.null(block_length)) {
    stop("`block_length` only applies to method = \"block\" or ",
         "\"circular\". Set the method, or leave `block_length` as NULL.",
         call. = FALSE)
  }
  new_design("bootstrap", list(
    n_replicates = check_count(n_replicates, "n_replicates", allow_zero = FALSE),
    n = if (is.null(n)) NULL else check_count(plan_n(n), "n"),
    method = method,
    block_length = if (is.null(block_length)) NULL else
      check_count(block_length, "block_length", allow_zero = FALSE),
    output = output
  ))
}

#' @export
draw_design.drawn_design_bootstrap <- function(design, data) {
  validate_data(data)
  # A design saved before `output` existed returns rows, as it did then.
  output <- design$output %||% "rows"
  if (output == "rows" && ".replicate" %in% names(data)) {
    stop("`data` already has a column called `.replicate`, which labels the ",
         "bootstrap replicates. Rename it.", call. = FALSE)
  }

  n_rows <- nrow(data)
  n <- design$n %||% n_rows
  if (output == "rows" && as.numeric(n) * design$n_replicates > 1e7) {
    warn_drawn(paste0(
      "`design_bootstrap()` returns ",
      fmt_n(as.numeric(n) * design$n_replicates), " rows here, ",
      fmt_n(design$n_replicates), " replicates of ", fmt_n(n), ". ",
      "output = \"indices\" or \"weights\" gives the same replicates as a ",
      "matrix of row numbers or counts."
    ), "drawn_bootstrap_size")
  }

  idx <- bootstrap_index(design, n_rows, n)
  switch(output,
    indices = idx,
    weights = {
      counts <- vapply(seq_len(ncol(idx)),
                       function(r) tabulate(idx[, r], nbins = n_rows),
                       integer(n_rows))
      dim(counts) <- c(n_rows, ncol(idx))
      counts
    },
    rows = {
      out <- reindex(data, as.vector(idx))
      k <- ncol(out)
      # cbind() on a data.frame demotes a tibble, and every other design keeps
      # the input's class, so the column is built in place. It is moved to the
      # front by position, since selecting by name drops a duplicated column
      # name from a tibble. A data frame's `[[<-` and `[` make duplicated
      # names unique, so the names are set last.
      out[[".replicate"]] <- rep(seq_len(ncol(idx)), each = n)
      out <- out[, c(k + 1L, seq_len(k)), drop = FALSE]
      names(out) <- c(".replicate", names(data))
      out
    }
  )
}

#' The rows each bootstrap replicate draws
#'
#' An integer matrix with `n` rows and one column per replicate. Replicates
#' are drawn one after another, with the same calls to the random number
#' generator as when each replicate was built and subset on its own.
#'
#' @noRd
bootstrap_index <- function(design, n_rows, n) {
  method <- design$method
  if (method != "simple") {
    block_length <- bootstrap_block_length(design, n_rows)
    if (block_length > n_rows) {
      stop("`block_length` (", block_length, ") exceeds the number of rows (",
           n_rows, ").", call. = FALSE)
    }
    n_blocks <- ceiling(n / block_length)
    offset <- rep(seq_len(block_length) - 1L, times = n_blocks)
    keep <- seq_len(n)
  }
  one <- switch(method,
    simple = function() sample.int(n_rows, size = n, replace = TRUE),
    block = function() {
      starts <- sample.int(n_rows - block_length + 1L, n_blocks, replace = TRUE)
      (rep(starts, each = block_length) + offset)[keep]
    },
    circular = function() {
      starts <- sample.int(n_rows, n_blocks, replace = TRUE)
      # Past the last row, a run carries on from the first. In doubles, so a
      # frame of more than 2^30 rows cannot overflow.
      at <- (rep(as.numeric(starts), each = block_length) - 1 + offset) %%
        n_rows + 1
      as.integer(at[keep])
    }
  )
  idx <- vapply(seq_len(design$n_replicates), function(r) one(), integer(n))
  dim(idx) <- c(n, design$n_replicates)
  idx
}

#' The run length of a block bootstrap: as given, else floor(N^(1/3))
#' @noRd
bootstrap_block_length <- function(design, n_rows) {
  # 64^(1/3) is 3.9999999999999996 in floating point, so the small offset
  # keeps a perfect cube from rounding down.
  design$block_length %||% max(1L, as.integer(floor(n_rows^(1 / 3) + 1e-9)))
}

#' What a refusal says about how often each row appears
#' @noRd
bootstrap_counts <- function(design) {
  if (identical(design$method, "block")) {
    return(paste0(
      "Under method = \"block\" rows do not all appear equally often. Runs ",
      "must fit inside the frame, so the first and last block_length - 1 ",
      "rows appear less often than the rest; under method = \"circular\" ",
      "every row appears n_replicates * n / nrow(data) times on average."
    ))
  }
  paste0("If you want the expected number of times each row appears, that is ",
         "n_replicates * n / nrow(data).")
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_bootstrap <- function(design, data) {
  stop("A bootstrap resamples the sample; it is not a probability sample of a ",
       "finite population, so it has no inclusion probability.\n",
       bootstrap_counts(design), call. = FALSE)
}
