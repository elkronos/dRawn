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
#' @section Variance:
#' The joint inclusion probabilities of a pivotal design have no closed form,
#' and treating the sample as simple random would give away exactly the gain
#' the spreading bought. [ht_total()] and [ht_mean()] instead use the local
#' mean variance estimator (Stevens and Olsen 2003; Grafström and Schelin
#' 2014), which compares each sampled row with the mean of itself and its three
#' nearest sampled neighbours. It is an approximation, reported as
#' `method = "local mean"`, and errs on the conservative side: in the
#' package's simulations it runs from about 1 to 1.6 times the true variance,
#' higher for small samples and large sampling fractions.
#'
#' @section Cost:
#' Each step finds a nearest neighbour among the rows still undecided, so the
#' draw takes time proportional to `N^2`: under a second for 5,000 rows, and
#' several seconds at 20,000. For much larger frames, spread within strata or
#' thin the frame first.
#'
#' @param across Column names to spread the sample across: two coordinates for
#'   a map, or any numeric auxiliary variables. At least one.
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
  check_key_columns(data, cols, "across")
  for (nm in cols) {
    if (!is.numeric(data[[nm]])) {
      stop("`", nm, "` must be numeric to spread across or size by, not ",
           class(data[[nm]])[1], ".", call. = FALSE)
    }
  }
  bad <- Reduce(`|`, lapply(data[cols], is.na))
  check_na_policy(bad, design$na_rm,
                  paste0("a missing value in ",
                         paste0("`", cols, "`", collapse = ", ")))
  keep <- which(!bad)
  if (!length(keep)) stop("Every row has a missing value.", call. = FALSE)

  x <- as.matrix(data[keep, design$across, drop = FALSE])
  storage.mode(x) <- "double"
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
#' @noRd
local_mean_variance <- function(x, y, pi_i, k = 4L) {
  n <- length(y)
  if (n < 2L) {
    return(list(variance = NA_real_, method = "local mean",
                note = "A variance needs at least two sampled rows."))
  }
  u <- y / pi_i
  x <- as.matrix(x)
  k <- min(k, n)
  dev <- vapply(seq_len(n), function(i) {
    d <- colSums((t(x) - x[i, ])^2)
    d[i] <- -Inf                             # itself first
    nb <- order(d)[seq_len(k)]
    kk <- length(nb)
    kk / (kk - 1) * (u[i] - mean(u[nb]))^2
  }, numeric(1))
  list(variance = sum(dev), method = "local mean",
       note = paste0("Local mean approximation for a spatially balanced ",
                     "sample: the local pivotal method has no closed-form ",
                     "joint inclusion probabilities."))
}

#' @noRd
ht_variance.drawn_design_spread <- function(design, data, rows, y, pi_i) {
  fr <- spread_frame(design, data)
  pos <- match(rows, fr$keep)
  if (anyNA(pos)) {
    return(list(variance = NA_real_, note = paste0(
      "Some sampled rows were dropped by this design's `na_rm`, so they have ",
      "no place in the space the variance is computed over.")))
  }
  local_mean_variance(fr$x[pos, , drop = FALSE], y, pi_i)
}
