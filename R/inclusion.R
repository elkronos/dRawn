#' Inclusion probabilities and design weights
#'
#' `inclusion_prob()` returns the first-order inclusion probability of every row
#' of `data` under `design` — the probability that the row lands in a sample —
#' without drawing one. `sampling_weight()` returns `1 / inclusion_prob()`, the
#' number of population units each sampled row stands for.
#'
#' These are what make a sample usable for estimation. A Horvitz-Thompson total
#' is `sum(y / pi)` over the sampled rows, and its unbiasedness rests entirely on
#' `pi` being the design's real inclusion probability rather than a plausible
#' guess.
#'
#' @section Which designs have a closed form:
#' Most do, and those are computed exactly:
#'
#' \tabular{ll}{
#'   [design_simple()]      \tab `n / N` \cr
#'   [design_stratified()]  \tab `n_h / N_h` within each stratum \cr
#'   [design_systematic()]  \tab `1 / interval`, for every row, exactly \cr
#'   [design_cluster()]     \tab `n_clusters / N_clusters` \cr
#'   [design_multistage()]  \tab `(n_clusters / N_clusters) * (n_h / N_h)`, equal allocation only \cr
#'   [design_weighted()]    \tab `n * p_i` for the `"systematic"` and `"poisson"` methods \cr
#'   [design_reservoir()]   \tab `n / N` \cr
#'   [design_temporal()]    \tab `per_interval / N_bucket` within each interval \cr
#'   [design_spatial()]     \tab `n / N_in_region` \cr
#' }
#'
#' Four cases have no closed form, and the package refuses to invent one:
#'
#' * `design_cluster(balanced = TRUE)` — the per-cluster take is the smallest
#'   *selected* cluster's size, which is itself random. Simulation on clusters of
#'   2/4/6/8 gives 0.50, 0.41, 0.34, 0.25, against the 0.50 a naive
#'   `n_clusters / N_clusters` would claim for every row.
#' * `design_multistage(allocation = "proportional")` — the stage-two allocation
#'   depends on which clusters were selected. Simulation on clusters of 3/5/7/9
#'   gives 0.17, 0.20, 0.17, 0.15, against a naive 0.33, 0.20, 0.14, 0.11.
#' * `design_weighted(method = "successive")` — the default. Successive sampling
#'   has no closed-form inclusion probability; this is the whole reason the other
#'   two methods exist. Note that its realised probabilities are not merely
#'   unknown but genuinely different from proportional-to-size: for weights
#'   `1:10` at `n = 5`, simulation puts them up to 35% away.
#' * [design_bootstrap()] — resampling with replacement from the sample is not a
#'   probability sample of a finite population, so there is no `pi` to report.
#'
#' For the first three, pass `simulate = TRUE` to estimate them by Monte Carlo
#' instead. The estimate carries `R`-sized error, which is fine for checking a
#' design and not fine for publishing a variance.
#'
#' @param data A data frame.
#' @param design A design object.
#' @param simulate Estimate the probabilities by repeated draws rather than in
#'   closed form. Required for the designs listed above; allowed for any design,
#'   which is a convenient way to check the exact formulas.
#' @param R Number of simulated draws when `simulate = TRUE`.
#' @param seed Optional seed for the simulation.
#'
#' @return A numeric vector with one element per row of `data`. Rows that the
#'   design can never select — outside the region, outside the time window —
#'   get `0`.
#'
#' @examples
#' df <- data.frame(id = 1:20, site = rep(c("a", "b"), times = c(15, 5)))
#'
#' # 8 rows allocated proportionally: 6 of 15, then 2 of 5
#' inclusion_prob(df, design_stratified("site", n = 8))
#'
#' sampling_weight(df, design_stratified("site", n = 8))
#'
#' # The default weighted design has no closed form; ask for a simulation
#' w <- data.frame(id = 1:5, w = c(1, 1, 1, 1, 16))
#' inclusion_prob(w, design_weighted("w", n = 2), simulate = TRUE, R = 2000,
#'                seed = 1)
#'
#' # Systematic PPS does have one, and it is exactly proportional to size
#' inclusion_prob(w, design_weighted("w", n = 2, method = "systematic"))
#'
#' @seealso [draw()], which attaches these to a sample when `weights = TRUE`.
#' @export
inclusion_prob <- function(data, design, simulate = FALSE, R = 5000,
                           seed = NULL) {
  if (!is_design(design)) {
    stop("`design` must come from one of the design_*() constructors, not a ",
         class(design)[1], ". See ?designs.", call. = FALSE)
  }
  validate_data(data)
  check_flag(simulate, "simulate")

  if (isTRUE(simulate)) {
    return(simulate_inclusion(data, design, R = check_count(R, "R",
                                                            allow_zero = FALSE),
                              seed = seed))
  }
  exact_inclusion(design, data)
}

#' @rdname inclusion_prob
#' @export
sampling_weight <- function(data, design, simulate = FALSE, R = 5000,
                          seed = NULL) {
  p <- inclusion_prob(data, design, simulate = simulate, R = R, seed = seed)
  out <- rep(NA_real_, length(p))
  out[p > 0] <- 1 / p[p > 0]
  out
}

#' Monte Carlo inclusion probabilities
#'
#' Draws `R` samples and counts how often each row appears. Works for every
#' design, including the four with no closed form.
#'
#' @noRd
simulate_inclusion <- function(data, design, R, seed) {
  n_rows <- nrow(data)
  key <- ".drawn_row_id"
  if (key %in% names(data)) {
    stop("`data` already has a column called `", key,
         "`, which the simulation needs. Rename it.", call. = FALSE)
  }
  tagged <- data
  tagged[[key]] <- seq_len(n_rows)

  hits <- integer(n_rows)
  with_seed(seed, {
    for (i in seq_len(R)) {
      s <- draw_design(design, tagged)
      ids <- s[[key]]
      # A row drawn twice in one replicate is still one inclusion.
      hits[unique(ids)] <- hits[unique(ids)] + 1L
    }
  })
  hits / R
}

#' @noRd
no_closed_form <- function(what, alternative) {
  stop(what, " has no closed-form inclusion probability.\n",
       alternative, "\n",
       "Or pass simulate = TRUE to estimate it by Monte Carlo.", call. = FALSE)
}

#' @noRd
exact_inclusion <- function(design, data) {
  UseMethod("exact_inclusion")
}

#' @noRd
exact_inclusion.default <- function(design, data) {
  stop("No inclusion probability method for design type \"",
       design_type(design), "\".", call. = FALSE)
}

