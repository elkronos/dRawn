#' Spatially balanced sampling
#'
#' Draws a sample that is spread evenly across the space of one or more
#' columns — map coordinates, or any numeric variables known for the whole
#' frame — using the local pivotal method (Grafström, Lundström and Schelin
#' 2012). Neighbouring rows compete for selection, so a sampled row makes its
#' nearest neighbours less likely to be sampled too, and the sample avoids the
#' clumps and gaps a simple random sample leaves behind.
#'
#' The inclusion probabilities are exactly the ones asked for — `n / N` for
#' every row, or proportional to `size` — so the sample is estimable in the
#' usual way. What changes is *which combinations* of rows can occur: samples
#' that bunch up in one corner become rare. When what you measure varies
#' smoothly across the space, as most environmental, geographic and
#' operational quantities do, that is worth a large reduction in variance at no
#' extra cost.
#'
#' @section Why this, and not a region:
#' [design_spatial()] restricts a simple random sample to a polygon; it decides
#' *where* sampling happens, not *how evenly*. `design_spread()` decides how
#' evenly, and needs no polygon or `sf` — only the columns to spread across. The
#' two compose naturally: filter the frame to your region first, then spread.
#'
#' @section Spreading over more than space:
#' Nothing requires the columns to be coordinates. Spreading across auxiliary
#' variables that predict the outcome — last year's value, size, age — balances
#' the sample on them the way stratification would, without choosing strata
#' boundaries (Grafström and Schelin 2014). With `scale = TRUE`, the default,
#' each column is standardised first so that one measured in larger units does
#' not decide the distances alone; set `scale = FALSE` for coordinates that
#' already share a unit, such as projected metres.
#'
#' A Date column in `across` counts in days and a date-time (POSIXct) column
#' in seconds since 1970, so time can be spread over alongside space. Keep
#' `scale = TRUE` when mixing them with coordinates, or the column in the
#' larger unit decides the distances. `size` has to be numeric.
#'
#' A POINT geometry column, such as the geometry of an `sf` frame, counts as
#' its two coordinates, read with [sf::st_coordinates()]: `across =
#' "geometry"` spreads over the map. Distances are measured in the
#' geometry's own units, so longitude and latitude are treated as plain
#' numbers; a projected CRS gives distances in metres.
#'
#' @section Variance:
#' The joint inclusion probabilities of a pivotal design have no closed form,
#' and treating the sample as simple random would give away exactly the gain
#' the spreading bought. [ht_total()] and [ht_mean()] instead use the local
#' mean variance estimator (Stevens and Olsen 2003; Grafström and Schelin
#' 2014), which compares each sampled row with the mean of itself and its three
#' nearest sampled neighbours. It is an approximation, reported as
#' `method = "local mean"`, and usually overstates the variance. The table
#' gives the ratio of its average to the true variance of the total, over
#' 1,000 draws from one frame per cell, by `n` and the sampling fraction
#' `f = n / N`. The frames are uniform on the unit square, spread with
#' `scale = FALSE`; the smooth surface is the one in the examples, and the
#' other variable is normal noise unrelated to position.
#'
#' | Variable | `n` | `f` = 0.05 | 0.2 | 0.5 | 0.8 |
#' |:--|--:|--:|--:|--:|--:|
#' | smooth surface | 8 | 2.2 | 2.5 | 2.3 | 3.6 |
#' | | 20 | 1.6 | 1.5 | 2.3 | 1.4 |
#' | | 40 | 1.1 | 1.3 | 1.9 | 2.1 |
#' | unrelated to position | 8 | 1.0 | 1.4 | 1.8 | 4.7 |
#' | | 20 | 1.0 | 1.2 | 1.8 | 4.5 |
#' | | 40 | 1.0 | 1.3 | 1.9 | 4.9 |
#'
#' For a variable unrelated to position the ratio is close to `1 / (1 - f)`,
#' since the estimator then reduces to the with-replacement variance of a
#' simple random sample. It also varies from frame to frame: with `n = 20`
#' and `f = 0.05`, four frames of 400 rows gave 0.96 to 1.15 for that
#' variable and 1.36 to 1.52 for the smooth surface, so a ratio slightly
#' below 1 occurs too. Below about 10 sampled rows the ratio climbs: `n = 4`
#' gave 2.4 for the smooth surface at `f = 0.05`, and 18.9 on a frame of 8
#' rows.
#'
#' [deff()] divides this variance by a simple random sample's, so it carries
#' the same overstatement: at a ratio of 2, a spread sample with a true design
#' effect of 0.8 reports 1.6, worse than a simple random sample when it is
#' better. Rows taken with certainty, because their `size` puts them past it
#' or because `n` is the whole frame, are in every sample and add no variance,
#' so the estimator leaves them out.
#'
#' @section Strata:
#' With `strata`, `n` is split over the strata as [design_stratified()]
#' splits it, by `allocation`, and each stratum is then spread on its own, as
#' if it were a frame of its own: `scale` standardises the columns within the
#' stratum, `size` sets probabilities proportional to size within it, and a
#' row competes only with rows of its own stratum. [inclusion_prob()] gives
#' `n_h / N_h` within stratum `h`, or `n_h * p_i` with `size`, exactly.
#'
#' The strata are drawn independently, so [ht_total()] and [ht_mean()] add up
#' the local mean variance of each stratum, and the degrees of freedom are
#' the sampled rows below certainty less the strata that hold them. A stratum
#' with one sampled row below certainty has no neighbour to compare it with,
#' and the variance is refused with the stratum named, as a stratified
#' sample's is; draw at least two rows in every stratum. [as_svydesign()]
#' hands the strata to `survey`, [allocation_table()] shows the split before
#' anything is drawn, and [sample_summary()] lists the rows drawn from each
#' stratum.
#'
#' Drawing each stratum with a design of its own and binding the samples with
#' [rbind()] does not give a sample to estimate from. The bound sample keeps
#' the first draw's frame and design, so [ht_total()] refuses it; use
#' `strata` instead.
#'
#' @section Reserves:
#' `reserve` lists that many more rows after the main sample, to visit in
#' place of rows found unavailable in the field, as [design_simple()] does.
#' The reserves are a second spread sample, of `reserve` rows at equal
#' probability from the rows the main sample left, so they cover the space
#' as evenly as the main sample does, and they come in random order.
#' [draw()] returns the main sample in frame order followed by the reserves,
#' with an integer `.reserve` column: 0 for the main sample, and 1, 2, ...
#' for the reserves, in the order they are used. The main sample is the one
#' the design without `reserve` draws with the same seed.
#'
#' [use_reserves()] takes the list and what the visits found, and returns
#' the sample to estimate from. It weights the rows used as a spread sample
#' of `n` from the frame without the `u` rows found unavailable, at
#' probability `n / (N - u)`. That is the average of the rows' probabilities,
#' not each row's own, so the weights are an approximation. The variance is
#' the local mean estimator, as for any spread sample. Only an unstratified
#' design without `size`, whose rows all have probability `n / N`, takes
#' `reserve`.
#'
#' @section Cost:
#' Each step finds a nearest neighbour among the rows still undecided, so the
#' draw takes time proportional to `N^2`: under a second for 5,000 rows, and
#' several seconds at 20,000. With `strata`, the walk runs within each
#' stratum and the time grows with the square of each stratum's size: 20,000
#' rows in four strata of 5,000 took under a third of the time they took as
#' one frame. For much larger frames, use strata, or thin the frame first.
#'
#' @param across Column names to spread the sample across: two coordinates for
#'   a map, or any numeric auxiliary variables. At least one. A Date or
#'   date-time column, and a POINT geometry column, are read as numbers; see
#'   "Spreading over more than space".
#' @param n Number of rows to draw, the total across `strata` when given, or
#'   a plan from [plan_size()].
#'   Fixed, not expected. With counts in `allocation` it can be left out and
#'   is their sum; when given, it must equal that sum.
#' @param size Optional column of positive size measures. When given, rows are
#'   included with probability proportional to size — `n * size / sum(size)`,
#'   with rows past certainty taken whole, exactly as
#'   [design_weighted()] — and spread at the same time.
#' @param scale Standardise each column of `across` to unit standard deviation
#'   before measuring distance. See "Spreading over more than space".
#' @param na_rm Drop rows with a missing value in `across`, `size` or
#'   `strata` instead of raising an error.
#' @param strata Optional column name(s) defining strata, to spread the
#'   sample within each one. Several columns are cross-classified. See
#'   "Strata".
#' @param allocation How `n` is split across `strata`: `"proportional"` to
#'   each stratum's size, or `"equal"`, by the largest-remainder method of
#'   [design_stratified()], with a stratum that would get more rows than it
#'   holds taken whole; or the count for each stratum, as a named vector such
#'   as `c(north = 10, south = 5)` when `strata` is one column, or a data
#'   frame of the `strata` columns and `n`. Counts follow the rules of
#'   [design_stratified()]: every stratum of the frame needs one, a count for
#'   a stratum the frame does not have is refused, and no count may exceed
#'   the rows in its stratum. Only with `strata`.
#' @param reserve Number of reserve rows to list after the main sample, to
#'   replace rows found unavailable in the field. Not available with `size`
#'   or `strata`. See "Reserves".
#'
#' @return A design object, for use with [draw()].
#'
#' @references
#' Deville, J.-C. and Tillé, Y. (1998). Unequal probability sampling without
#' replacement through a splitting method. *Biometrika*, 85, 89–101.
#'
#' Grafström, A., Lundström, N. L. P. and Schelin, L. (2012). Spatially
#' balanced sampling through the pivotal method. *Biometrics*, 68, 514–520.
#'
#' Grafström, A. and Schelin, L. (2014). How to select representative samples.
#' *Scandinavian Journal of Statistics*, 41, 277–290.
#'
#' Stevens, D. L. and Olsen, A. R. (2003). Variance estimation for spatially
#' balanced samples of environmental resources. *Environmetrics*, 14,
#' 593–610.
#'
#' Stevens, D. L. and Olsen, A. R. (2004). Spatially balanced sampling of
#' natural resources. *Journal of the American Statistical Association*, 99,
#' 262–278.
#'
#' @examples
#' set.seed(1)
#' plots <- data.frame(x = stats::runif(400), y = stats::runif(400))
#' # Something that varies smoothly over the map
#' plots$biomass <- 100 + 80 * plots$x + 60 * sin(4 * plots$y) +
#'   stats::rnorm(400, 0, 10)
#'
#' spread <- design_spread(c("x", "y"), n = 40, scale = FALSE)
#' s <- draw(plots, spread, seed = 1, weights = TRUE)
#' ht_total(s, "biomass")
#' sum(plots$biomass)   # the truth
#'
#' # The same total from a simple random sample: a wider interval
#' ht_total(draw(plots, design_simple(n = 40), seed = 1, weights = TRUE),
#'          "biomass")
#'
#' # See the difference: spread points avoid each other
#' op <- par(mfrow = c(1, 2))
#' plot(design_simple(n = 40), plots, type = "map", coords = c("x", "y"),
#'      seed = 1)
#' plot(spread, plots, type = "map", seed = 1)
#' par(op)
#'
#' # Spread within two regions, each given rows in proportion to its size
#' plots$region <- ifelse(plots$x < 0.5, "west", "east")
#' by_region <- design_spread(c("x", "y"), n = 40, scale = FALSE,
#'                            strata = "region")
#' allocation_table(plots, by_region)
#' ht_total(draw(plots, by_region, seed = 1, weights = TRUE), "biomass")
#'
#' @family designs
#' @seealso [draw()], [design_spatial()] to restrict sampling to a region,
#'   [design_weighted()] for size-proportional selection without spreading,
#'   [design_stratified()] for the allocation rules, [allocation_table()],
#'   [use_reserves()].
#' @export
design_spread <- function(across, n, size = NULL, scale = TRUE,
                          na_rm = FALSE, strata = NULL,
                          allocation = c("proportional", "equal"),
                          reserve = 0) {
  if (!missing(n)) n <- plan_n(n)
  across <- check_columns(across, "across")
  if (is.null(strata)) {
    if (!missing(allocation)) {
      stop("`allocation` splits `n` over strata, and `strata` is NULL. Name ",
           "the strata columns in `strata`, or leave `allocation` out.",
           call. = FALSE)
    }
    allocation <- NULL
    n <- check_count(n, "n")
  } else {
    strata <- check_columns(strata, "strata")
    if (is.numeric(allocation) || is.data.frame(allocation)) {
      allocation <- check_stratum_counts(allocation, strata)
      n <- counts_n(if (!missing(n)) n, allocation)
    } else {
      if (!is.character(allocation)) {
        stop("`allocation` must be \"proportional\", \"equal\" or the count ",
             "for each stratum, not ", format_bad(allocation), ".",
             call. = FALSE)
      }
      allocation <- check_choice(allocation, c("proportional", "equal"),
                                 "allocation")
      if (missing(n) || is.null(n)) {
        stop("`n` is missing. Give the total number of rows to draw, or the ",
             "count for each stratum in `allocation`.", call. = FALSE)
      }
      n <- check_count(n, "n")
    }
  }
  size <- if (is.null(size)) NULL else check_columns(size, "size", 1L)
  scale <- check_flag(scale, "scale")
  na_rm <- check_flag(na_rm, "na_rm")
  reserve <- check_count(reserve, "reserve")
  if (reserve > 0L && (!is.null(size) || !is.null(strata))) {
    by <- if (!is.null(size)) "size" else "strata"
    stop("`reserve` needs every row at the same probability, and `", by,
         "` ", if (by == "size") {
           "makes it proportional to size"
         } else {
           "splits `n` over strata"
         }, ". Leave out `", by, "` or `reserve`.", call. = FALSE)
  }
  new_design("spread", list(
    across = across,
    strata = strata,
    n = n,
    allocation = allocation,
    size = size,
    scale = scale,
    na_rm = na_rm,
    reserve = reserve
  ))
}

#' Rows a spread design can use, their target probabilities and positions
#'
#' Shared by the draw and the probability and variance methods, so all three
#' agree on which rows are dropped, how the columns are scaled, and what each
#' row's probability is.
#'
#' `keep` holds the frame rows the design can use, and `x` and `pi` the
#' scaled positions and target probabilities of those rows. The rest is by
#' position in `keep`: `idx` lists the rows of each stratum, in key order,
#' `code` gives each row's stratum as a position in `labels`, `sizes` and
#' `n_alloc`, and `n_alloc` is the rows drawn from each stratum. A design
#' without `strata` has one stratum, with no `code` or `labels`. Each stratum
#' is scaled and sized on its own, as a frame of its own would be, and its
#' split of `n` comes from `stratum_allocation()`, as a stratified design's
#' does.
#'
#' @noRd
spread_frame <- function(design, data) {
  cols <- c(design$across, design$size)
  strata <- design$strata
  validate_data(data, required_columns = c(cols, strata))
  geom <- vapply(cols, function(nm) inherits(data[[nm]], "sfc"), logical(1))
  check_key_columns(data, cols[!geom], "across")
  if (!is.null(strata)) check_key_columns(data, strata, "strata")
  for (nm in design$across) {
    col <- data[[nm]]
    if (!is.numeric(col) && !inherits(col, c("Date", "POSIXt", "sfc"))) {
      stop("`", nm, "` must be numeric to spread across, not ",
           class(col)[1], ". A Date, date-time or POINT geometry column ",
           "works too; convert anything else to numbers first, for example ",
           "with as.numeric().", call. = FALSE)
    }
  }
  if (!is.null(design$size) && !is.numeric(data[[design$size]])) {
    stop("`", design$size, "` must be numeric to size by, not ",
         class(data[[design$size]])[1], ".", call. = FALSE)
  }
  x_cols <- coordinate_columns(data, design$across, "across", time = TRUE)
  bad <- spread_missing(x_cols,
                        if (!is.null(design$size)) data[[design$size]])
  if (!is.null(strata)) bad <- bad | is.na(group_key(data, strata))
  check_na_policy(bad, design$na_rm,
                  paste0("a missing value in ",
                         paste0("`", c(cols, strata), "`", collapse = ", ")))
  keep <- which(!bad)
  if (!length(keep)) stop("Every row has a missing value.", call. = FALSE)
  # Rows use_reserves() found unavailable are outside the design's frame.
  if (!is.null(design[["unavailable"]])) {
    keep <- intersect(keep, available_rows(design, nrow(data)))
  }

  x <- matrix(vapply(x_cols, function(v) as.numeric(v[keep]),
                     numeric(length(keep))),
              nrow = length(keep), dimnames = list(NULL, names(x_cols)))
  if (any(!is.finite(x))) {
    stop("`across` contains non-finite values.", call. = FALSE)
  }

  if (is.null(strata)) {
    if (design$n > length(keep)) {
      stop("`n` (", design$n, ") cannot exceed the number of rows (",
           length(keep), ").", call. = FALSE)
    }
    check_reserve_size(design, length(keep))
    al <- list(idx = list(seq_along(keep)), sizes = length(keep),
               n_alloc = design$n)
  } else {
    al <- stratum_allocation(design, table_frame(
      lapply(strata, function(col) data[[col]][keep]), strata, list()))
  }
  w <- NULL
  if (!is.null(design$size)) {
    w <- data[[design$size]][keep]
    check_weights(w, design$size)
  }

  pi <- numeric(length(keep))
  for (h in seq_along(al$idx)) {
    i <- al$idx[[h]]
    if (design$scale) {
      s <- apply(x[i, , drop = FALSE], 2, stats::sd)
      s[!is.finite(s) | s == 0] <- 1
      x[i, ] <- sweep(x[i, , drop = FALSE], 2, s, "/")
    }
    pi[i] <- if (is.null(w)) {
      al$n_alloc[h] / length(i)
    } else {
      pps_pi(w[i], al$n_alloc[h])
    }
  }
  list(keep = keep, x = x, pi = pi, idx = al$idx, code = al$code,
       labels = al$labels, sizes = al$sizes, n_alloc = al$n_alloc)
}

#' Rows with a missing value in any column a spread design reads
#'
#' `x_cols` is `coordinate_columns()` of `across`, so an empty point counts
#' as missing; `size` is that column, or `NULL`.
#'
#' @noRd
spread_missing <- function(x_cols, size = NULL) {
  Reduce(`|`, lapply(c(x_cols, if (!is.null(size)) list(size)), is.na))
}

#' @export
draw_design.drawn_design_spread <- function(design, data) {
  fr <- spread_frame(design, data)
  if (!is.null(design$strata)) {
    warn_empty_groups(fr$n_alloc, fr$sizes, fr$labels, fix = paste0(
      if (has_counts(design)) "Give " else "Raise `n`, or give ",
      "each stratum a count of at least 1 in `allocation` (or 2, so each ",
      "stratum also supports a variance)."))
  }
  # Strata are drawn one after another in key order, each by its own
  # pivotal walk.
  picked <- unlist(lapply(fr$idx, function(i) {
    i[local_pivotal(fr$pi[i], fr$x[i, , drop = FALSE])]
  }), use.names = FALSE)
  r <- reserve_size(design)
  if (r == 0L) return(reindex(data, fr$keep[picked], sort = TRUE))

  # The reserves are a second pivotal sample, of `reserve` rows at equal
  # probability from the rows the main sample left, so they are spread over
  # the space too. They are used in random order.
  left <- setdiff(seq_along(fr$keep), picked)
  spare <- left[local_pivotal(rep(r / length(left), length(left)),
                              fr$x[left, , drop = FALSE])]
  spare <- spare[sample.int(length(spare))]
  out <- reindex(data, fr$keep[c(sort(picked), spare)])
  with_reserve(out, design$n, r)
}

#' The local pivotal method, LPM2
#'
#' Repeatedly picks an undecided row at random and lets it compete with its
#' nearest undecided neighbour: the pair's probabilities are pooled and one of
#' them is pushed to 0 or 1, in a way that leaves each one's expected value
#' unchanged (Deville and Tillé 1998). Every row's inclusion probability is
#' therefore exactly its target, while neighbours rarely end up in the sample
#' together. Ties for nearest neighbour are broken at random.
#'
#' Every row with a target above 0 competes, however small the target. Each
#' step settles one of the pair at exactly 0 or 1, so the loop ends after at
#' most one step per row. A pooled probability less than `eps` below 1 is
#' the rounding error of a sum that is 1 in exact arithmetic and is settled
#' at 1, and a remainder less than `eps` above 1 is settled at 0 for the same
#' reason. A small pooled probability is kept, so two rows of `1e-13` pool to
#' `2e-13` and not to 0.
#'
#' @noRd
local_pivotal <- function(pi, x, eps = 1e-12) {
  p <- pi
  idx <- which(p > 0 & p < 1)
  while (length(idx) > 1L) {
    i <- idx[sample.int(length(idx), 1L)]
    others <- idx[idx != i]
    d <- colSums((t(x[others, , drop = FALSE]) - x[i, ])^2)
    nearest <- others[d == min(d)]
    j <- if (length(nearest) > 1L) nearest[sample.int(length(nearest), 1L)]
         else nearest
    s <- p[i] + p[j]
    if (s < 1) {
      pooled <- if (s > 1 - eps) 1 else s
      if (stats::runif(1) < p[j] / s) {
        p[j] <- pooled; p[i] <- 0
      } else {
        p[i] <- pooled; p[j] <- 0
      }
    } else {
      left <- if (s - 1 < eps) 0 else s - 1
      if (stats::runif(1) < (1 - p[j]) / (2 - s)) {
        p[i] <- 1; p[j] <- left
      } else {
        p[j] <- 1; p[i] <- left
      }
    }
    idx <- idx[p[idx] > 0 & p[idx] < 1]
  }
  # A single row left over only when the probabilities did not sum to a whole
  # number; settle it with its own remaining probability.
  if (length(idx) == 1L) {
    p[idx] <- as.numeric(stats::runif(1) < p[idx])
  }
  which(p == 1)
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_spread <- function(design, data) {
  fr <- spread_frame(design, data)
  out <- numeric(nrow(data))
  out[fr$keep] <- fr$pi
  out
}

#' @noRd
joint_inclusion.drawn_design_spread <- function(design, data, rows) {
  no_joint_form(
    "`design_spread()`",
    paste0("ht_total() uses the local mean variance estimator, which needs ",
           "no joint probabilities.")
  )
}

#' Local mean variance estimator for a spatially balanced sample
#'
#' Each sampled row is compared with the mean of `y / pi` over its
#' neighbourhood -- itself and its three nearest sampled neighbours, in the
#' same scaled space the design spread across -- and the squared deviations are
#' summed with a `k / (k - 1)` correction for the neighbourhood's own size
#' (Stevens and Olsen 2003; Grafström and Schelin 2014).
#'
#' Where the variable is unrelated to position this reduces to the
#' with-replacement variance of a simple random sample, so it errs
#' conservative as the sampling fraction grows. A per-row `(1 - pi)` correction
#' was tried and rejected: it removed that excess for unstructured variables
#' but left smooth ones understated by about a quarter at a 20-40% sampling
#' fraction, and understating is the worse failure.
#'
#' `nb` is `local_neighbours(x, k)` when the caller has it already, as the
#' domains of one sample do.
#'
#' @noRd
local_mean_variance <- function(x, y, pi_i, k = 4L, nb = NULL) {
  n <- length(y)
  if (n < 2L) {
    return(list(variance = NA_real_, method = "local mean",
                note = "A variance needs at least two sampled rows."))
  }
  u <- y / pi_i
  nb <- nb %||% local_neighbours(x, k)
  kk <- ncol(nb)
  local <- rowMeans(matrix(u[nb], nrow = n))
  dev <- kk / (kk - 1) * (u - local)^2
  list(variance = sum(dev), method = "local mean",
       note = paste0("Local mean approximation for a spatially balanced ",
                     "sample: the local pivotal method has no closed-form ",
                     "joint inclusion probabilities."))
}

#' Each sampled row's neighbourhood for the local mean
#'
#' Row `i` of the result holds the row itself and its `k - 1` nearest rows of
#' `x`, by squared distance, `k` at most the number of rows. The
#' neighbourhoods depend on where the rows lie and not on `y`.
#'
#' @noRd
local_neighbours <- function(x, k = 4L) {
  x <- as.matrix(x)
  n <- nrow(x)
  k <- min(k, n)
  nb <- vapply(seq_len(n), function(i) {
    d <- colSums((t(x) - x[i, ])^2)
    d[i] <- -Inf                             # itself first
    order(d)[seq_len(k)]
  }, integer(k))
  matrix(nb, nrow = n, ncol = k, byrow = TRUE)
}

#' Variance of a spread sample
#'
#' Rows taken with certainty, in a census or because a `size` column puts them
#' past it, are in every possible sample and add no variance. Kept in, each
#' one's distance from its neighbours' mean counted as sampling variation. A
#' census then reported a standard error of about 5% of its exact total, and
#' on a frame with three dominant sizes the estimate averaged 1,270 times the
#' true variance. They are left out, both as rows and as neighbours, and the
#' estimator runs over the rest.
#'
#' With `strata`, each stratum is a spread sample of its own, drawn
#' independently of the others, so the variance is the sum of the strata's
#' local mean variances, with neighbours found within each stratum.
#'
#' @noRd
ht_variance.drawn_design_spread <- function(design, data, rows, y, pi_i,
                                            domain = NULL, ctx = NULL) {
  fr <- from_context(ctx, "spread", spread_frame(design, data))
  pos <- match(rows, fr$keep)
  if (anyNA(pos)) {
    return(list(variance = NA_real_, note = paste0(
      "Some sampled rows were dropped by this design's `na_rm`, so they have ",
      "no place in the space the variance is computed over.")))
  }
  free <- !is_certain(pi_i)
  if (!any(free)) return(all_certain())
  if (!is.null(design$strata)) {
    return(spread_strata_variance(design, data, fr, rows, pos, y, pi_i, free,
                                  domain, ctx))
  }
  if (sum(free) < 2L) {
    return(list(variance = NA_real_, method = "local mean", note = paste0(
      "Only one sampled row was not taken with certainty, so the local mean ",
      "estimator has no neighbour to compare it with. Draw a larger `n`.")))
  }
  x <- fr$x[pos[free], , drop = FALSE]
  local_mean_variance(x, y[free], pi_i[free],
                      nb = from_context(ctx, "neighbours", local_neighbours(x)))
}

#' The local mean variance of a spread sample with strata
#'
#' The sum over strata of `local_mean_variance()` on each stratum's sampled
#' rows below certainty. A stratum with one such row has no neighbour to
#' compare it with, and is refused by `lonely_note()`, unless `domain` shows
#' it holds no unit of the domain, when it adds 0. `fr` is `spread_frame()`
#' of the frame `data`, `pos` the sampled rows' positions in `fr$keep`, and
#' `free` marks the sampled rows below certainty.
#'
#' @noRd
spread_strata_variance <- function(design, data, fr, rows, pos, y, pi_i,
                                   free, domain = NULL, ctx = NULL) {
  k <- length(fr$sizes)
  # Each frame row's stratum, counting only the rows below certainty.
  code <- rep(NA_integer_, nrow(data))
  open <- !is_certain(fr$pi)
  code[fr$keep[open]] <- fr$code[open]
  group <- structure(code, levels = fr$labels, class = "factor")
  note <- lonely_note(group, rows, c("stratum", "strata"), lonely_fix(design),
                      domain)
  if (!is.null(note)) {
    return(list(variance = NA_real_, method = "local mean", note = note))
  }
  h <- fr$code[pos]
  by_stratum <- lapply(seq_len(k), function(g) which(free & h == g))
  nb <- from_context(ctx, "neighbours", lapply(by_stratum, function(i) {
    if (length(i) >= 2L) local_neighbours(fr$x[pos[i], , drop = FALSE])
  }))
  v <- 0
  for (g in seq_len(k)) {
    i <- by_stratum[[g]]
    if (length(i) < 2L) next
    v <- v + local_mean_variance(fr$x[pos[i], , drop = FALSE], y[i],
                                 pi_i[i], nb = nb[[g]])$variance
  }
  list(variance = v, method = "local mean",
       note = paste0("Local mean approximation for a spatially balanced ",
                     "sample, summed over strata: the local pivotal method ",
                     "has no closed-form joint inclusion probabilities."))
}
