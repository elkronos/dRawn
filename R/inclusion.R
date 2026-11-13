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
#'   [design_systematic()]  \tab `1 / interval` with a random `start`, or `n / N` with `n`, for every row walked, exactly \cr
#'   [design_cluster()]     \tab `n_clusters / N_clusters`; with `balanced = TRUE`, that times the cluster's average take over its size, worked out exactly (see [design_cluster()]) \cr
#'   [design_multistage()]  \tab `(n_clusters / N_clusters) * (n_h / N_h)`, when the per-cluster take `n_h` is constant; with replacement, `(n_clusters / N_clusters) * (1 - (1 - 1 / N_h)^n_h)`; with `size`, the cluster's `pi_h = n_clusters * M_h / sum(M)` in place of `n_clusters / N_clusters` \cr
#'   [design_weighted()]    \tab `n * p_i` for the `"systematic"` and `"poisson"` methods; with `strata`, `n_h * p_i` within each stratum \cr
#'   [design_certainty()]   \tab `1` above the threshold, `rest`'s own below it \cr
#'   [design_reservoir()]   \tab `n / min(N, max_items)`, and `0` past `max_items` \cr
#'   [design_temporal()]    \tab `per_interval / N_bucket` within each interval \cr
#'   [design_spatial()]     \tab `n / N_in_region` \cr
#'   [design_spread()]      \tab `n / N`, or `n * p_i` with a `size` column; with `strata`, the same within each stratum \cr
#' }
#'
#' A design with `reserve` gives the probabilities of its main sample, which
#' the reserves do not change. The design of a sample from [use_reserves()]
#' gives `n / (N - u)`, and `0` for the `u` rows found unavailable.
#'
#' Four cases have no closed form, and the package refuses to invent one:
#'
#' * `design_multistage(allocation = "proportional")` — the stage-two allocation
#'   depends on which clusters were selected. Simulation on clusters of 3/5/7/9
#'   gives 0.17, 0.20, 0.17, 0.15, against a naive 0.33, 0.20, 0.14, 0.11.
#' * `design_multistage()` where the per-cluster take is not constant — `n` not
#'   divisible by `n_clusters`, or, without replacement, a cluster smaller than
#'   `n / n_clusters`. The allocation runs over the *selected* clusters, so the
#'   remainder goes to the largest of them and an undersized cluster is capped
#'   with its shortfall dealt to whichever clusters came with it. A row's
#'   probability then depends on the company it keeps. On clusters of
#'   20/10/8/6 with `n = 7` over 2 clusters, the average take is out by up to
#'   17%, and a cluster smaller than the take can push the figure above 1.
#'
#'   For these two, the exact probability is the probability given the
#'   clusters selected, averaged over every possible selection, and that
#'   average has no formula. `draw(weights = TRUE)` weights the sample
#'   conditional on the clusters it selected instead, which gives an
#'   unbiased total without the average; see "Design weights" in [draw()].
#' * `design_weighted(method = "successive")` — the default. Successive sampling
#'   has no closed-form inclusion probability; this is the whole reason the other
#'   two methods exist. Note that its realised probabilities are not merely
#'   unknown but genuinely different from proportional-to-size: for weights
#'   `1:10` at `n = 5`, simulation puts them up to 35% away.
#' * [design_bootstrap()] — resampling with replacement from the sample is not a
#'   probability sample of a finite population, so there is no `pi` to report.
#'
#' A [design_systematic()] with a fixed `start` is refused too, for another
#' reason. It is deterministic: each row has probability 0 or 1, and
#' `simulate = TRUE` reports exactly that. Without `simulate` those values
#' are refused, because `draw(weights = TRUE)` would give every sampled row a
#' weight of 1 and [ht_total()] would report the sum over the sample as the
#' total of the frame. So are [sampling_weight()], [joint_prob()],
#' `draw(weights = TRUE)` and the estimators. Leave `start` as `NULL` for a
#' probability sample.
#'
#' For all but the bootstrap, pass `simulate = TRUE` to estimate them by Monte
#' Carlo instead. The estimate carries `R`-sized error, which is fine for
#' checking a design and not fine for publishing a variance. Simulating a
#' bootstrap is refused rather than answered: every row turns up in some
#' replicate, so the count converges to 1 for all of them.
#'
#' @section What a simulation costs:
#' `simulate = TRUE` makes `R` draws, 5,000 by default, so it takes `R`
#' times as long as one [draw()]. Most designs draw from a frame of 10,000
#' rows in well under a second. The local pivotal method of [design_spread()]
#' does work in proportion to `N^2` on every draw, so at 10,000 rows each
#' draw takes seconds and the default `R` takes hours. When the first draws
#' project a run of more than 30 seconds, a message says how long it will
#' take, and an interactive session shows a progress bar. Lower `R` for a
#' quicker answer with more Monte Carlo error: its standard error is
#' `sqrt(p * (1 - p) / R)`, at most 0.007 for `R = 5000`.
#'
#' A warning that every draw raises, such as a stratum allocated no rows or
#' `max_items` cutting a frame short, is given once at the end, with the
#' number of draws that raised it.
#'
#' @param data A data frame.
#' @param design A design object.
#' @param simulate Estimate the probabilities by repeated draws rather than in
#'   closed form. Required for the designs listed above; allowed for any
#'   probability design, which is a convenient way to check the exact formulas.
#' @param R Number of simulated draws when `simulate = TRUE`. See "What a
#'   simulation costs".
#' @param seed Optional seed for the simulation, applied as in [draw()] (see
#'   "Seeding" there). Checked whether or not `simulate` is `TRUE`.
#'
#' @return A numeric vector with one element per row of `data`.
#'
#'   Rows the design can never select — outside the region, outside the time
#'   window, past `max_items`, or dropped under `na_rm = TRUE` — get `0` from
#'   `inclusion_prob()` and `NA` from `sampling_weight()`, since no finite
#'   number of population rows is represented by a row that cannot be drawn.
#'   [sample_summary()] counts them. A weight of 0 is not one of these. The
#'   designs that select in proportion to size refuse it, so drop such rows
#'   from `data` first.
#'
#'   A design that [draw()] refuses, such as an `n` larger than the frame or
#'   than the rows inside a region, is refused here too, with the same
#'   message, rather than given probabilities above 1 or a census.
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
  check_design(design, data, "inclusion_prob")
  validate_data(data)
  check_flag(simulate, "simulate")
  check_seed(seed)
  warn_region(design, data)

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

#' Rows a size-proportional design takes with certainty
#'
#' A design whose inclusion probabilities are proportional to size gives a
#' row whose share of `n` reaches 1 a probability of exactly 1, and spreads
#' what is left of `n` over the other rows, repeating until no share passes
#' 1. `certainty_set()` says which rows of `data` those are, and where the
#' cut-off falls, before any sample is drawn.
#'
#' @param data A data frame.
#' @param design [design_weighted()] with `method = "systematic"` or
#'   `"poisson"`, or [design_spread()] with a `size` column. Other designs are
#'   refused, since none of them has inclusion probabilities proportional to
#'   size.
#'
#' @return An integer vector of the rows of `data` whose inclusion probability
#'   is 1, in frame order, the same rows as
#'   `which(inclusion_prob(data, design) == 1)`. Its attribute `threshold` is
#'   the smallest size among them, or `NA` when there are none. Of the rows
#'   the design can select, those whose size is at least `threshold` are the
#'   set. With `strata`, [design_weighted()] and [design_spread()] put the
#'   cut-off in each stratum at its own size, so there `threshold` is `NA`.
#'   For [design_weighted()] without `strata`, [design_certainty()] with
#'   that `threshold` on the weights column therefore takes the same rows
#'   with certainty, and a `rest` of the same
#'   method with `n` reduced by the number of rows in the set gives the other
#'   rows the same probabilities too.
#'
#' @examples
#' df <- data.frame(id = 1:10, w = c(1, 2, 3, 4, 5, 6, 7, 8, 60, 90))
#' d <- design_weighted("w", n = 4, method = "systematic")
#' cs <- certainty_set(df, d)
#' df[cs, ]
#' attr(cs, "threshold")
#'
#' # The same probabilities, with the cut-off made explicit
#' same <- design_certainty("w", attr(cs, "threshold"),
#'                          rest = design_weighted("w", n = 4 - length(cs),
#'                                                 method = "systematic"))
#' all.equal(inclusion_prob(df, same), inclusion_prob(df, d))
#'
#' @seealso [inclusion_prob()], [design_certainty()]
#' @export
certainty_set <- function(data, design) {
  if (!is_design(design)) {
    stop("`design` must come from one of the design_*() constructors, not ",
         a_class(design), ". See ?designs.", call. = FALSE)
  }
  size <- if (inherits(design, "drawn_design_weighted") &&
              design$method %in% c("systematic", "poisson")) {
    design$weights
  } else if (inherits(design, "drawn_design_spread")) {
    design$size
  }
  if (is.null(size)) {
    stop("`design` is ", design_label(design), ", whose inclusion ",
         "probabilities are not proportional to size. certainty_set() takes ",
         "design_weighted() with method = \"systematic\" or \"poisson\", or ",
         "design_spread() with `size`.", call. = FALSE)
  }
  validate_data(data)
  rows <- which(is_certain(exact_inclusion(design, data)))
  # With strata, each stratum has a cut-off of its own.
  structure(rows, threshold = if (length(rows) && is.null(design$strata)) {
    min(data[[size]][rows])
  } else {
    NA_real_
  })
}

#' Monte Carlo inclusion probabilities
#'
#' Draws `R` samples and counts how often each row appears. Works for every
#' design, including those with no closed form.
#'
#' @noRd
simulate_inclusion <- function(data, design, R, seed) {
  refuse_simulation(design)
  # The probabilities of a design with `reserve` are its main sample's.
  design <- main_design(design)
  n_rows <- nrow(data)
  key <- ".drawn_row_id"
  if (key %in% names(data)) {
    stop("`data` already has a column called `", key,
         "`, which the simulation needs. Rename it.", call. = FALSE)
  }
  tagged <- data
  tagged[[key]] <- seq_len(n_rows)

  hits <- integer(n_rows)
  simulate_draws(R, seed, function() {
    # A row drawn twice in one replicate is still one inclusion.
    ids <- unique(drawn_ids(draw_design(design, tagged), key, design))
    hits[ids] <<- hits[ids] + 1L
  })
  hits / R
}

#' The frame rows of one simulated draw
#'
#' Read from the column the simulators add to `data`. A draw method that
#' drops it, as one written outside the package might, is refused rather
#' than counted as drawing nothing.
#'
#' @noRd
drawn_ids <- function(s, key, design) {
  ids <- if (is.data.frame(s)) s[[key]] else unlist(s)
  if (is.null(ids)) {
    stop("The draw_design() method for design type \"", design_type(design),
         "\" returned rows without the `", key, "` column the simulation ",
         "adds to `data`. A method must return rows of `data` with all its ",
         "columns; see ?draw_design.", call. = FALSE)
  }
  ids
}

#' Run `R` simulated draws
#'
#' `one()` makes one draw and records it, and the draws run under `seed`.
#' A warning of drawn's own class is muffled and counted, and each kind is
#' given once at the end with the number of draws that raised it, rather than
#' once per draw. A run is timed over its first five draws, or over the draws
#' of its first two seconds when they are slower, and one projected to take
#' more than 30 seconds says so with a message and, in an interactive
#' session, shows a progress bar.
#'
#' @noRd
simulate_draws <- function(R, seed, one) {
  i <- 0L
  kinds <- character(0)
  first <- character(0)
  draws <- integer(0)
  last <- integer(0)
  record <- function(w) {
    k <- match(class(w)[1L], kinds)
    if (is.na(k)) {
      kinds <<- c(kinds, class(w)[1L])
      first <<- c(first, conditionMessage(w))
      draws <<- c(draws, 0L)
      last <<- c(last, 0L)
      k <- length(kinds)
    }
    if (last[k] != i) {
      draws[k] <<- draws[k] + 1L
      last[k] <<- i
    }
    invokeRestart("muffleWarning")
  }

  bar <- NULL
  on.exit(if (!is.null(bar)) close(bar), add = TRUE)
  start <- elapsed_time()
  timed <- FALSE
  with_seed(seed, withCallingHandlers({
    for (i in seq_len(R)) {
      one()
      if (!timed) {
        spent <- elapsed_time() - start
        if (i >= min(5L, R) || spent >= 2) {
          timed <- TRUE
          bar <- announce_run(R, i, spent)
        }
      }
      if (!is.null(bar)) utils::setTxtProgressBar(bar, i)
    }
  }, drawn_warning = record))

  for (k in seq_along(kinds)) {
    warn_drawn(paste0(first[k], " This happened in ", fmt_n(draws[k]),
                      " of the ", fmt_n(R), " simulated draws."), kinds[k])
  }
  invisible(NULL)
}

#' Say how long a simulation will take, when that is more than 30 seconds
#'
#' `done` draws took `spent` seconds. Returns a progress bar, already at
#' `done`, in an interactive session, else `NULL`.
#'
#' @noRd
announce_run <- function(R, done, spent) {
  total <- spent / done * R
  if (done >= R || total <= 30) return(NULL)
  message("Simulating ", fmt_n(R), " draws will take about ",
          format_duration(total), ", at ", format_duration(spent / done),
          " a draw. A smaller `R` takes less time, with more Monte Carlo ",
          "error.")
  if (!is_interactive()) return(NULL)
  utils::txtProgressBar(min = 0, max = R, initial = done, style = 3)
}

#' @noRd
no_closed_form <- function(what, alternative) {
  stop_declined(paste0(
    what, " has no closed-form inclusion probability.\n", alternative, "\n",
    "Or estimate it by Monte Carlo with ",
    "inclusion_prob(data, design, simulate = TRUE)."))
}

#' @noRd
exact_inclusion <- function(design, data) {
  UseMethod("exact_inclusion")
}

#' @noRd
exact_inclusion.default <- function(design, data) {
  stop_declined(paste0(
    "Design type \"", design_type(design), "\" is not one of drawn's own, ",
    "so drawn has no inclusion probabilities for it, and `weights = TRUE`, ",
    "the exact inclusion_prob() and joint_prob() and the estimators are not ",
    "available. inclusion_prob(data, design, simulate = TRUE) estimates the ",
    "probabilities by Monte Carlo; see ?draw_design."))
}

