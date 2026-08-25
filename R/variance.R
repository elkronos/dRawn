#' Joint inclusion probabilities
#'
#' The probability that a *pair* of rows both land in the sample. First-order
#' probabilities from [inclusion_prob()] give you an unbiased total;
#' second-order probabilities are what let you put a standard error on it.
#'
#' @section Which designs have them:
#' Closed forms exist, and are used, for the designs whose selection is either
#' independent across groups or a simple random sample within them:
#'
#' \tabular{ll}{
#'   [design_simple()]      \tab `n(n-1) / (N(N-1))` for a pair, without replacement \cr
#'   [design_stratified()]  \tab within a stratum as above; across strata, independent \cr
#'   [design_cluster()]     \tab same cluster: `a/A`; different clusters: `a(a-1)/(A(A-1))` \cr
#'   [design_multistage()]  \tab the two stages multiplied, equal allocation only \cr
#'   [design_reservoir()]   \tab identical to simple random sampling \cr
#'   [design_temporal()]    \tab within an interval as above; across intervals, independent \cr
#'   [design_spatial()]     \tab simple random sampling inside the region \cr
#'   [design_weighted()]    \tab `"poisson"` only, where rows are independent: `pi_i * pi_j` \cr
#'   [design_systematic()]  \tab `1/interval` for rows sharing a residue class, otherwise **zero** \cr
#' }
#'
#' Systematic sampling is the awkward one. Most pairs can never co-occur, so
#' their joint probability is genuinely 0 and no design-unbiased variance
#' estimator exists. [ht_total()] says so rather than returning a number.
#'
#' `design_weighted(method = "systematic")` has joint probabilities, but they
#' depend on the order units are visited and need a dedicated algorithm. Use
#' `sampling::UPsystematicpi2()` for those.
#'
#' @param data A data frame.
#' @param design A design object.
#' @param rows Optional row indices. Supply these — usually the rows you drew —
#'   to get the submatrix for them instead of the full `nrow(data)` square,
#'   which is what makes this usable on a large population.
#'
#' @return A square matrix with one row and column per element of `rows`
#'   (or per row of `data`). The diagonal holds first-order probabilities.
#'
#' @examples
#' df <- data.frame(id = 1:10)
#' round(joint_prob(df, design_simple(n = 4)), 3)
#'
#' # Only for the rows you drew
#' joint_prob(df, design_simple(n = 4), rows = c(2, 5, 7))
#'
#' @seealso [inclusion_prob()], [ht_total()]
#' @export
joint_prob <- function(data, design, rows = NULL) {
  if (!is_design(design)) {
    stop("`design` must come from one of the design_*() constructors, not a ",
         class(design)[1], ". See ?designs.", call. = FALSE)
  }
  validate_data(data)
  rows <- rows %||% seq_len(nrow(data))
  if (!is.numeric(rows) || anyNA(rows) || any(rows < 1) ||
      any(rows > nrow(data))) {
    stop("`rows` must be valid row indices into `data`.", call. = FALSE)
  }
  joint_inclusion(design, data, as.integer(rows))
}

#' @noRd
no_joint_form <- function(what, alternative) {
  stop(what, " has no closed-form joint inclusion probability.\n", alternative,
       call. = FALSE)
}

#' @noRd
joint_inclusion <- function(design, data, rows) UseMethod("joint_inclusion")

#' @noRd
joint_inclusion.default <- function(design, data, rows) {
  no_joint_form(
    paste0("`", design_type(design), "`"),
    "Only the designs listed in ?joint_prob have one."
  )
}

#' Pairwise matrix from a group id and a per-group (n, N)
#'
#' Rows in the same group are a simple random sample of `n_g` from `N_g`; rows
#' in different groups are selected independently.
#'
#' @noRd
joint_from_groups <- function(group, n_g, size_g, pi_i) {
  k <- length(group)
  out <- outer(pi_i, pi_i)                     # independent case
  same <- outer(group, group, "==")
  if (any(same)) {
    n <- n_g[as.character(group)]
    N <- size_g[as.character(group)]
    within <- outer(seq_len(k), seq_len(k), function(a, b) {
      nn <- n[a]; NN <- N[a]
      ifelse(NN > 1, nn * (nn - 1) / (NN * (NN - 1)), 0)
    })
    out[same] <- within[same]
  }
  diag(out) <- pi_i
  dimnames(out) <- NULL
  out
}

#' @noRd
joint_inclusion.drawn_design_simple <- function(design, data, rows) {
  if (design$replace) {
    no_joint_form("`design_simple(replace = TRUE)`",
                  "Use replace = FALSE, or estimate variance by resampling.")
  }
  N <- nrow(data)
  n <- design$n
  pi_i <- rep(n / N, length(rows))
  out <- matrix(if (N > 1) n * (n - 1) / (N * (N - 1)) else 0,
                length(rows), length(rows))
  diag(out) <- pi_i
  out
}

#' @noRd
joint_inclusion.drawn_design_reservoir <- function(design, data, rows) {
  joint_inclusion.drawn_design_simple(
    new_design("simple", list(n = min(design$n, nrow(data)), replace = FALSE)),
    data, rows
  )
}

#' @noRd
joint_inclusion.drawn_design_stratified <- function(design, data, rows) {
  if (design$replace) {
    no_joint_form("`design_stratified(replace = TRUE)`", "Use replace = FALSE.")
  }
  validate_data(data, required_columns = design$strata)
  check_key_columns(data, design$strata, "strata")

  keys <- data[design$strata]
  bad <- Reduce(`|`, lapply(keys, is.na))
  group <- interaction(keys, drop = TRUE, lex.order = TRUE)
  idx <- split(seq_len(nrow(data))[!bad], group[!bad])
  sizes <- lengths(idx)
  n_alloc <- allocate(design$n, sizes, design$allocation,
                      design$min_per_stratum, cap = TRUE,
                      spread = stratum_spread(design, data, idx))

  pi_all <- exact_inclusion(design, data)
  joint_from_groups(as.character(group[rows]), n_alloc, sizes, pi_all[rows])
}

#' @noRd
joint_inclusion.drawn_design_temporal <- function(design, data, rows) {
  validate_data(data, required_columns = design$time)
  tz <- design$tz %||% "UTC"
  tv <- parse_time(data[[design$time]], tz, design$time)
  breaks <- make_breaks(parse_time(design$from, tz, "from"),
                        parse_time(design$to, tz, "to"),
                        design$interval, design$unit)
  pi_all <- exact_inclusion(design, data)

  bucket <- rep(NA_character_, nrow(data))
  if (length(breaks) >= 2L) {
    inw <- !is.na(tv) & tv >= breaks[1] & tv < breaks[length(breaks)]
    bucket[inw] <- as.character(findInterval(tv[inw], breaks))
  }
  sizes <- table(bucket[!is.na(bucket)])
  n_g <- pmin(sizes, design$per_interval)

  g <- bucket[rows]
  g[is.na(g)] <- paste0("__out__", seq_len(sum(is.na(g))))  # never co-occur
  joint_from_groups(g, stats::setNames(as.integer(n_g), names(n_g)),
                    stats::setNames(as.integer(sizes), names(sizes)),
                    pi_all[rows])
}

#' @noRd
joint_inclusion.drawn_design_spatial <- function(design, data, rows) {
  require_suggested("sf", "spatial sampling")
  inside <- spatial_inside(design, data)
  N <- sum(inside)
  n <- min(design$n, N)
  pi_all <- exact_inclusion(design, data)
  out <- matrix(if (N > 1) n * (n - 1) / (N * (N - 1)) else 0,
                length(rows), length(rows))
  out[!inside[rows], ] <- 0
  out[, !inside[rows]] <- 0
  diag(out) <- pi_all[rows]
  out
}

#' @noRd
joint_inclusion.drawn_design_cluster <- function(design, data, rows) {
  if (design$balanced) {
    no_joint_form("`design_cluster(balanced = TRUE)`",
                  "Use balanced = FALSE.")
  }
  cl <- count_clusters(design, data)
  A <- cl$total
  a <- design$n_clusters
  lab <- as.character(cl$labels[rows])

  same <- outer(lab, lab, "==")
  same[is.na(same)] <- FALSE
  out <- matrix(if (A > 1) a * (a - 1) / (A * (A - 1)) else 0,
                length(rows), length(rows))
  out[same] <- a / A
  pi_all <- exact_inclusion(design, data)
  diag(out) <- pi_all[rows]
  out
}

#' @noRd
joint_inclusion.drawn_design_multistage <- function(design, data, rows) {
  if (design$allocation == "proportional") {
    no_joint_form("`design_multistage(allocation = \"proportional\")`",
                  "Use allocation = \"equal\".")
  }
  if (design$replace) {
    no_joint_form("`design_multistage(replace = TRUE)`", "Use replace = FALSE.")
  }
  if (design$n %% design$n_clusters != 0L) {
    no_joint_form(
      paste0("`design_multistage()` with n = ", design$n, " over ",
             design$n_clusters, " clusters"),
      paste0("The per-cluster take is not constant, so pairs in different\n",
             "clusters have no single joint probability. Choose an `n` that\n",
             "divides evenly by `n_clusters` (here, a multiple of ",
             design$n_clusters, ").")
    )
  }
  cl <- count_clusters(design, data)
  A <- cl$total
  a <- design$n_clusters
  m <- design$n %/% design$n_clusters
  sizes <- table(cl$labels[cl$present])

  lab <- as.character(cl$labels[rows])
  Nh <- as.numeric(sizes[lab])
  same <- outer(lab, lab, "==")
  same[is.na(same)] <- FALSE

  # Different clusters: both clusters selected, then each row within its own.
  out <- outer(m / Nh, m / Nh) * (if (A > 1) a * (a - 1) / (A * (A - 1)) else 0)
  # Same cluster: that cluster selected, then both rows drawn from it.
  win <- outer(seq_along(rows), seq_along(rows), function(i, j) {
    NN <- Nh[i]
    ifelse(NN > 1, (a / A) * (m * (m - 1)) / (NN * (NN - 1)), 0)
  })
  out[same] <- win[same]
  pi_all <- exact_inclusion(design, data)
  diag(out) <- pi_all[rows]
  out
}

#' @noRd
joint_inclusion.drawn_design_weighted <- function(design, data, rows) {
  if (design$method == "poisson") {
    pi_all <- exact_inclusion(design, data)
    out <- outer(pi_all[rows], pi_all[rows])   # independent
    diag(out) <- pi_all[rows]
    return(out)
  }
  if (design$method == "systematic") {
    no_joint_form(
      "`design_weighted(method = \"systematic\")`",
      paste0("Its joint probabilities depend on the order units are visited ",
             "and need a\ndedicated algorithm. `sampling::UPsystematicpi2()` ",
             "computes them.")
    )
  }
  no_joint_form("`design_weighted(method = \"successive\")`",
                "It has no closed-form first-order probability either.")
}

#' @noRd
joint_inclusion.drawn_design_systematic <- function(design, data, rows) {
  if (!is.null(design$start)) {
    stop("A systematic design with a fixed `start` is not a probability ",
         "sample.", call. = FALSE)
  }
  k <- design$interval
  res <- (rows - 1L) %% k
  out <- matrix(0, length(rows), length(rows))
  out[outer(res, res, "==")] <- 1 / k
  diag(out) <- 1 / k
  out
}

#' @noRd
joint_inclusion.drawn_design_bootstrap <- function(design, data, rows) {
  no_joint_form("`design_bootstrap()`",
                "A bootstrap is not a probability sample of a finite population.")
}


# ---- Horvitz-Thompson estimation -------------------------------------------

#' Estimate a population total from a sample
#'
#' Forms the Horvitz-Thompson total `sum(y / pi)` and, where the design allows
#' it, a design-unbiased variance and confidence interval.
#'
#' @section Variance:
#' For fixed-size designs the Sen-Yates-Grundy estimator is used, which is
#' non-negative more often than the general Horvitz-Thompson form and is the
#' usual choice. Poisson sampling has a random size, so the independent-units
#' form `sum((1 - pi) / pi^2 * y^2)` is used instead.
#'
#' A variance needs joint inclusion probabilities, and not every design has
#' them — see [joint_prob()]. Where they are unavailable the estimate is still
#' returned, with `variance` as `NA` and a note saying why. Systematic sampling
#' is the notable case: most pairs of rows can never co-occur, so no
#' design-unbiased variance estimator exists at all.
#'
#' @param sample A data frame returned by [draw()] with `weights = TRUE`.
#' @param y The variable to total: a column name, or a numeric vector as long as
#'   `sample`.
#' @param variance How to compute it. `"auto"` uses the analytic estimator when
#'   the design has one and falls back to the jackknife when it does not;
#'   `"analytic"` insists on the analytic form and errors otherwise;
#'   `"jackknife"` always resamples; `"none"` skips it. The result reports which
#'   was used.
#' @param level Confidence level for the interval.
#'
#' @return A list with `total`, `variance`, `se`, `ci`, `n`, `method` and
#'   `note`, with a `print()` method.
#'
#' @examples
#' set.seed(1)
#' pop <- data.frame(
#'   id = 1:200,
#'   site = rep(c("a", "b"), times = c(150, 50)),
#'   spend = round(stats::runif(200, 10, 500))
#' )
#'
#' s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
#' ht_total(s, "spend")
#'
#' sum(pop$spend)   # the truth
#'
#' @seealso [joint_prob()], [inclusion_prob()]
#' @export
ht_total <- function(sample, y, variance = c("auto", "analytic", "jackknife",
                                            "none"), level = 0.95) {
  variance <- match.arg(variance)
  if (!is.data.frame(sample)) {
    stop("`sample` must be a data frame returned by draw().", call. = FALSE)
  }
  design <- attr(sample, "drawn_design")
  rows <- attr(sample, "drawn_rows")
  pop <- attr(sample, "drawn_population")
  if (is.null(design) || is.null(rows) || is.null(pop)) {
    stop("`sample` does not carry its design. Draw it with ",
         "draw(..., weights = TRUE), and call ht_total() on that result ",
         "before subsetting it.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }

  yv <- if (is.character(y) && length(y) == 1L) {
    if (!y %in% names(sample)) {
      stop("`sample` has no column `", y, "`.", call. = FALSE)
    }
    sample[[y]]
  } else {
    y
  }
  if (!is.numeric(yv) || length(yv) != nrow(sample)) {
    stop("`y` must be a numeric column name or a numeric vector as long as ",
         "`sample`.", call. = FALSE)
  }
  if (anyNA(yv)) {
    stop(sum(is.na(yv)), " value(s) of `y` are missing.", call. = FALSE)
  }

  pi_i <- sample$.prob
  total <- sum(yv / pi_i)

  var_out <- switch(variance,
    none = list(variance = NA_real_, method = "none",
                note = "Variance not requested."),
    jackknife = jackknife_variance(design, sample, pop, yv, pi_i),
    analytic = tryCatch(
      c(ht_variance(design, pop, rows, yv, pi_i), list(method = "analytic")),
      error = function(e) list(variance = NA_real_, method = "analytic",
                               note = conditionMessage(e))),
    tryCatch(
      c(ht_variance(design, pop, rows, yv, pi_i), list(method = "analytic")),
      error = function(e) {
        jk <- jackknife_variance(design, sample, pop, yv, pi_i)
        jk$note <- paste0("No analytic variance for this design, so the ",
                          "jackknife was used instead. ",
                          sub("\n.*", "", conditionMessage(e)))
        jk
      })
  )

  se <- if (is.na(var_out$variance) || var_out$variance < 0) {
    NA_real_
  } else {
    sqrt(var_out$variance)
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  ci <- if (is.na(se)) c(NA_real_, NA_real_) else total + c(-1, 1) * z * se

  structure(
    list(total = total, variance = var_out$variance, se = se, ci = ci,
         level = level, n = nrow(sample), design = design_type(design),
         method = var_out$method %||% "analytic", note = var_out$note),
    class = "drawn_ht"
  )
}

#' Delete-a-group jackknife variance
#'
#' Works from the sample alone, which is what makes it available where the
#' analytic form is not. Groups are the primary sampling units: whole clusters
#' where the design has them, otherwise individual rows. Deleting a group,
#' inflating the surviving weights to compensate, and looking at how far the
#' estimate moves is a direct measure of how much the estimate depended on
#' which groups were drawn.
#'
#' @noRd
jackknife_variance <- function(design, sample, pop, y, pi_i) {
  if (inherits(design, "drawn_design_systematic")) {
    return(list(variance = NA_real_, method = "jackknife",
                note = paste0("A systematic sample has one primary sampling ",
                              "unit -- the random start -- so deleting rows ",
                              "does not reflect the design's randomness and ",
                              "the jackknife is not valid here.")))
  }
  grp <- jackknife_groups(design, sample)
  m <- length(unique(grp))
  if (m < 2L) {
    return(list(variance = NA_real_, method = "jackknife",
                note = paste0("The jackknife needs at least two primary ",
                              "sampling units; this sample has ", m, ".")))
  }
  base <- y / pi_i
  ug <- unique(grp)
  reps <- vapply(ug, function(g) sum(base[grp != g]) * m / (m - 1), numeric(1))

  # Without a finite population correction the jackknife treats the frame as
  # infinite and overstates the variance -- by a factor of four when a quarter
  # of the clusters were taken.
  fpc <- jackknife_fpc(design, pop, m)
  v <- fpc * ((m - 1) / m) * sum((reps - mean(reps))^2)
  list(variance = v, method = "jackknife",
       note = paste0("Jackknife over ", m, " primary sampling unit(s)",
                     if (fpc < 1) paste0(", with a finite population ",
                                         "correction of ", signif(fpc, 3))
                     else "", "."))
}

#' The share of primary sampling units NOT taken
#'
#' Applied for single-stage designs, where deleting a sampling unit accounts for
#' all the randomness there is. Deliberately not applied to multi-stage designs:
#' the delete-a-cluster jackknife sees only between-cluster variation, so the
#' second stage is already missing, and correcting for the first stage on top of
#' that understates the total. Leaving it out is the ultimate-cluster
#' approximation, which errs conservative -- the usual choice in survey practice.
#'
#' @noRd
jackknife_fpc <- function(design, pop, m) {
  if (inherits(design, "drawn_design_multistage")) return(1)
  total <- tryCatch({
    if (inherits(design, "drawn_design_cluster")) {
      count_clusters(design, pop)$total
    } else {
      nrow(pop)
    }
  }, error = function(e) NA_real_)
  if (!is.finite(total) || total <= m) return(1)
  1 - m / total
}

#' @noRd
jackknife_groups <- function(design, sample) {
  col <- if (inherits(design, c("drawn_design_cluster",
                                "drawn_design_multistage"))) {
    design$clusters
  } else if (inherits(design, "drawn_design_stratified")) {
    NULL   # rows within strata are the sampling units
  } else {
    NULL
  }
  if (!is.null(col) && col %in% names(sample)) {
    as.character(sample[[col]])
  } else {
    as.character(seq_len(nrow(sample)))
  }
}

#' @noRd
ht_variance <- function(design, data, rows, y, pi_i) {
  if (inherits(design, "drawn_design_systematic")) {
    stop("Systematic sampling has no design-unbiased variance estimator: most ",
         "pairs of rows can never appear together, so their joint inclusion ",
         "probability is 0. Repeat the draw with different starts, or treat ",
         "the sample as simple random, which is conservative.", call. = FALSE)
  }

  # Poisson: units are independent, so no joint matrix is needed.
  if (inherits(design, "drawn_design_weighted") && design$method == "poisson") {
    return(list(variance = sum((1 - pi_i) / pi_i^2 * y^2), note = NULL))
  }

  pij <- joint_inclusion(design, data, rows)
  k <- length(rows)
  if (k < 2L) {
    return(list(variance = NA_real_,
                note = "A variance needs at least two sampled rows."))
  }

  # Sen-Yates-Grundy, for fixed-size designs.
  yk <- y / pi_i
  d <- outer(yk, yk, "-")^2
  num <- outer(pi_i, pi_i) - pij
  w <- num / pij
  w[!is.finite(w)] <- 0                 # a zero joint probability contributes nothing
  diag(w) <- 0
  list(variance = 0.5 * sum(w * d), note = NULL)
}

#' @export
print.drawn_ht <- function(x, ...) {
  cat("Horvitz-Thompson total  (", x$design, " design, n = ", x$n, ")\n",
      sep = "")
  cat("  total    ", format(x$total, big.mark = ","), "\n", sep = "")
  if (is.na(x$se)) {
    cat("  se       NA\n")
    if (!is.null(x$note)) {
      cat("\n", strwrap(x$note, prefix = "  "), sep = "\n")
      cat("\n")
    }
  } else {
    cat("  se       ", format(x$se, big.mark = ","),
        "  (", x$method, ")\n", sep = "")
    cat("  ", format(100 * x$level), "% CI  ",
        format(x$ci[1], big.mark = ","), " to ",
        format(x$ci[2], big.mark = ","), "\n", sep = "")
    if (identical(x$method, "jackknife") && !is.null(x$note)) {
      cat("\n", strwrap(x$note, prefix = "  "), sep = "\n"); cat("\n")
    }
  }
  invisible(x)
}
