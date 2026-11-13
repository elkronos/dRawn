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
#' @section Cost:
#' Each step finds a nearest neighbour among the rows still undecided, so the
#' draw takes time proportional to `N^2`: under a second for 5,000 rows, and
#' several seconds at 20,000. For much larger frames, spread within strata or
#' thin the frame first.
#'
#' @param across Column names to spread the sample across: two coordinates for
#'   a map, or any numeric auxiliary variables. At least one. A Date or
#'   date-time column, and a POINT geometry column, are read as numbers; see
#'   "Spreading over more than space".
#' @param n Number of rows to draw. Fixed, not expected.
#' @param size Optional column of positive size measures. When given, rows are
#'   included with probability proportional to size — `n * size / sum(size)`,
#'   with rows past certainty taken whole, exactly as
#'   [design_weighted()] — and spread at the same time.
#' @param scale Standardise each column of `across` to unit standard deviation
#'   before measuring distance. See "Spreading over more than space".
#' @param na_rm Drop rows with a missing value in `across` or `size` instead of
#'   raising an error.
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
#' @family designs
#' @seealso [draw()], [design_spatial()] to restrict sampling to a region,
#'   [design_weighted()] for size-proportional selection without spreading.
#' @export
design_spread <- function(across, n, size = NULL, scale = TRUE,
                          na_rm = FALSE) {
  new_design("spread", list(
    across = check_columns(across, "across"),
    n = check_count(n, "n"),
    size = if (is.null(size)) NULL else check_columns(size, "size", 1L),
    scale = check_flag(scale, "scale"),
    na_rm = check_flag(na_rm, "na_rm")
  ))
}

#' Rows a spread design can use, their target probabilities and positions
#'
#' Shared by the draw and the probability and variance methods, so all three
#' agree on which rows are dropped, how the columns are scaled, and what each
#' row's probability is.
#'
#' @noRd
spread_frame <- function(design, data) {
  cols <- c(design$across, design$size)
  validate_data(data, required_columns = cols)
  geom <- vapply(cols, function(nm) inherits(data[[nm]], "sfc"), logical(1))
  check_key_columns(data, cols[!geom], "across")
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
  check_na_policy(bad, design$na_rm,
                  paste0("a missing value in ",
                         paste0("`", cols, "`", collapse = ", ")))
  keep <- which(!bad)
  if (!length(keep)) stop("Every row has a missing value.", call. = FALSE)

  x <- matrix(vapply(x_cols, function(v) as.numeric(v[keep]),
                     numeric(length(keep))),
              nrow = length(keep), dimnames = list(NULL, names(x_cols)))
  if (any(!is.finite(x))) {
    stop("`across` contains non-finite values.", call. = FALSE)
  }
  if (design$scale) {
    s <- apply(x, 2, stats::sd)
    s[!is.finite(s) | s == 0] <- 1
    x <- sweep(x, 2, s, "/")
  }

  if (design$n > length(keep)) {
    stop("`n` (", design$n, ") cannot exceed the number of rows (",
         length(keep), ").", call. = FALSE)
  }
  pi <- if (is.null(design$size)) {
    rep(design$n / length(keep), length(keep))
  } else {
    w <- data[[design$size]][keep]
    check_weights(w, design$size)
    pps_pi(w, design$n)
  }
  list(keep = keep, x = x, pi = pi)
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
  picked <- local_pivotal(fr$pi, fr$x)
  reindex(data, fr$keep[picked], sort = TRUE)
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
#' @noRd
local_pivotal <- function(pi, x, eps = 1e-9) {
  p <- pi
  p[p < eps] <- 0
  p[p > 1 - eps] <- 1
  open <- p > 0 & p < 1
  idx <- which(open)
  while (length(idx) > 1L) {
    i <- idx[sample.int(length(idx), 1L)]
    others <- idx[idx != i]
    d <- colSums((t(x[others, , drop = FALSE]) - x[i, ])^2)
    nearest <- others[d == min(d)]
    j <- if (length(nearest) > 1L) nearest[sample.int(length(nearest), 1L)]
         else nearest
    s <- p[i] + p[j]
    if (s < 1) {
      if (stats::runif(1) < p[j] / s) {
        p[j] <- s; p[i] <- 0
      } else {
        p[i] <- s; p[j] <- 0
      }
    } else {
      if (stats::runif(1) < (1 - p[j]) / (2 - s)) {
        p[i] <- 1; p[j] <- s - 1
      } else {
        p[j] <- 1; p[i] <- s - 1
      }
    }
    for (k in c(i, j)) {
      if (p[k] < eps) p[k] <- 0
      if (p[k] > 1 - eps) p[k] <- 1
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
    paste0("The local pivotal method has no closed-form joint inclusion ",
           "probability. Pass\nsimulate = TRUE to estimate them; ht_total() ",
           "uses the local mean variance\nestimator instead.")
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
  if (sum(free) < 2L) {
    return(list(variance = NA_real_, method = "local mean", note = paste0(
      "Only one sampled row was not taken with certainty, so the local mean ",
      "estimator has no neighbour to compare it with. Draw a larger `n`.")))
  }
  x <- fr$x[pos[free], , drop = FALSE]
  local_mean_variance(x, y[free], pi_i[free],
                      nb = from_context(ctx, "neighbours", local_neighbours(x)))
}
