est_pop <- function() {
  set.seed(404)
  N <- 300
  d <- data.frame(
    id   = 1:N,
    site = rep(c("a", "b", "c"), times = c(150, 100, 50)),
    cl   = rep(paste0("c", 1:30), each = 10),
    y    = round(stats::runif(N, 20, 400))
  )
  d$w <- d$y + 100
  d
}

test_that("ht_mean recovers the population mean under a range of designs", {
  d <- est_pop()
  truth <- mean(d$y)

  for (des in list(
    design_simple(n = 80),
    design_stratified("site", n = 80),
    design_weighted("w", n = 80, method = "poisson")
  )) {
    est <- replicate(400, {
      s <- draw(d, des, weights = TRUE)
      ht_mean(s, "y", variance = "none")$mean
    })
    se <- stats::sd(est) / sqrt(length(est))
    expect_lt(abs(mean(est) - truth), 4 * se)
  }
})

test_that("the two estimators differ, and Hajek is the steadier one", {
  d <- est_pop()
  # Unequal probabilities are what separate them; with equal ones they agree.
  eq <- draw(d, design_simple(n = 60), seed = 3, weights = TRUE)
  expect_equal(ht_mean(eq, "y", "hajek", variance = "none")$mean,
               ht_mean(eq, "y", "ht", variance = "none")$mean)

  des <- design_weighted("w", n = 60, method = "poisson")
  pair <- vapply(1:200, function(i) {
    s <- draw(d, des, seed = i, weights = TRUE)
    c(hajek = ht_mean(s, "y", "hajek", variance = "none")$mean,
      ht    = ht_mean(s, "y", "ht", variance = "none")$mean)
  }, numeric(2))
  expect_lt(stats::sd(pair["hajek", ]), stats::sd(pair["ht", ]))
})

test_that("the Hajek variance tracks the actual sampling variance", {
  d <- est_pop()
  des <- design_stratified("site", n = 90)
  out <- vapply(1:400, function(i) {
    s <- draw(d, des, seed = i, weights = TRUE)
    e <- ht_mean(s, "y")
    c(e$mean, e$variance)
  }, numeric(2))
  empirical <- stats::var(out[1, ])
  expect_gt(mean(out[2, ]) / empirical, 0.75)
  expect_lt(mean(out[2, ]) / empirical, 1.35)
})

test_that("the interval widens with the level and is centred on the estimate", {
  d <- est_pop()
  s <- draw(d, design_stratified("site", n = 60), seed = 11, weights = TRUE)
  a <- ht_mean(s, "y", level = 0.90)
  b <- ht_mean(s, "y", level = 0.99)
  expect_equal(a$mean, b$mean)
  expect_lt(diff(a$ci), diff(b$ci))
  expect_equal(mean(a$ci), a$mean)
})

test_that("ht_mean accepts a vector as well as a column name", {
  d <- est_pop()
  s <- draw(d, design_simple(n = 50), seed = 2, weights = TRUE)
  expect_equal(ht_mean(s, "y", variance = "none")$mean,
               ht_mean(s, s$y, variance = "none")$mean)
})

test_that("ht_mean reports why it has no variance rather than inventing one", {
  d <- est_pop()
  # One cluster: nothing to compare it with, analytically or by jackknife
  s <- draw(d, design_cluster("cl", n_clusters = 1), seed = 1, weights = TRUE)
  e <- ht_mean(s, "y")
  expect_true(is.finite(e$mean))
  expect_true(is.na(e$se))
  expect_true(nzchar(e$note))
  expect_output(print(e), "se\\s+NA")
})

test_that("ht_mean refuses samples it cannot estimate from", {
  d <- est_pop()
  bare <- draw(d, design_simple(n = 20), seed = 1)
  expect_error(ht_mean(bare, "y"), "weights = TRUE")

  s <- draw(d, design_simple(n = 20), seed = 1, weights = TRUE)
  expect_error(ht_mean(s, "nope"), "no column")
  expect_error(ht_mean(s, "site"), "numeric")
  expect_error(ht_mean(s, "y", level = 1), "between 0 and 1")

  na_y <- s
  na_y$y[1] <- NA
  expect_error(ht_mean(na_y, "y"), "missing")
})

test_that("deff falls below 1 when stratification helps and rises when it hurts", {
  set.seed(88)
  N <- 400
  d <- data.frame(
    id   = 1:N,
    site = rep(c("a", "b", "c", "d"), each = 100),
    cl   = rep(paste0("c", 1:40), each = 10)
  )
  # The grouping variable is the whole story of y, so stratifying on it removes
  # nearly all the variance and clustering on it keeps nearly all of it.
  d$y <- rep(c(20, 80, 160, 260), each = 100) + round(stats::rnorm(N, 0, 6))

  st <- deff(ht_total(draw(d, design_stratified("site", n = 40), seed = 1,
                           weights = TRUE), "y"))
  cl <- deff(ht_total(draw(d, design_cluster("cl", n_clusters = 4), seed = 1,
                           weights = TRUE), "y"))
  expect_lt(st, 0.5)
  expect_gt(cl, 2)
})

test_that("deff is available from both estimators and NA without a variance", {
  d <- est_pop()
  s <- draw(d, design_stratified("site", n = 60), seed = 5, weights = TRUE)
  expect_true(is.finite(deff(ht_total(s, "y"))))
  expect_true(is.finite(deff(ht_mean(s, "y"))))
  expect_true(is.na(deff(ht_mean(s, "y", variance = "none"))))
  expect_error(deff(1:10), "ht_total")
})

test_that("the jackknife deletes clusters, not rows", {
  # The primary sampling unit of a cluster design is the cluster. Deleting one
  # row at a time instead treats 80 correlated rows as 80 independent ones and
  # inflates the variance; the note names the count, so it is checkable.
  d <- est_pop()
  s <- draw(d, design_cluster("cl", n_clusters = 8), seed = 1, weights = TRUE)
  expect_equal(nrow(s), 80)

  for (r in list(ht_total(s, "y", variance = "jackknife"),
                 ht_mean(s, "y", variance = "jackknife"))) {
    expect_match(r$note, "Jackknife over 8 primary")
    expect_false(grepl("over 80 primary", r$note))
  }

  # And it lands on the analytic answer, which is the point of grouping right.
  expect_equal(ht_total(s, "y", variance = "jackknife")$variance,
               ht_total(s, "y", variance = "analytic")$variance,
               tolerance = 1e-6)
})

test_that("ht_mean falls back to the jackknife and says it did", {
  d <- est_pop()
  # One row per cluster: the exact two-stage form cannot see within-cluster
  # variation, so it declines and the delete-a-cluster jackknife steps in.
  s <- draw(d, design_multistage("cl", n_clusters = 10, n = 10), seed = 1,
            weights = TRUE)
  e <- ht_mean(s, "y")
  expect_equal(e$method, "jackknife")
  expect_true(is.finite(e$se))
  expect_match(e$note, "jackknife was used")
})

test_that("printing an estimate shows the pieces someone would report", {
  d <- est_pop()
  s <- draw(d, design_stratified("site", n = 60), seed = 5, weights = TRUE)
  out <- utils::capture.output(print(ht_mean(s, "y")))
  expect_true(any(grepl("mean", out)))
  expect_true(any(grepl("estimate", out)))
  expect_true(any(grepl("95% CI", out, fixed = TRUE)))
  expect_true(any(grepl("deff", out)))
})

test_that("deff is exactly 1 for a simple random sample", {
  # The reference variance divided an n-term sum by N - 1, leaving it short by
  # n(N-1)/(N(n-1)) -- so deff came out at 1 + 1/n and a plain SRS of 10
  # printed as "worse than simple random sampling".
  set.seed(31)
  pop <- data.frame(y = stats::rnorm(5000, 100, 20))
  for (n in c(10, 20, 50, 200)) {
    got <- vapply(1:200, function(i) {
      deff(ht_total(draw(pop, design_simple(n = n), seed = i, weights = TRUE),
                    "y"))
    }, numeric(1))
    expect_equal(mean(got), 1, tolerance = 1e-9,
                 label = paste("deff at n =", n))
  }
})

test_that("deff is NA rather than negative when there is no variance", {
  set.seed(5)
  pop <- data.frame(y = stats::rnorm(60, 10, 3), cl = rep(1:6, each = 10))
  r <- ht_total(draw(pop, design_cluster("cl", n_clusters = 1), seed = 1,
                     weights = TRUE), "y")
  expect_true(is.na(deff(r)))
  expect_true(is.na(drawn:::deff_value(-5, 1:10, rep(0.5, 10), 100, 10,
                                       "total")))
})

test_that("an empty sample is refused with a reason, not a crash", {
  # A start past the end of a short frame selects nothing
  d <- data.frame(y = 1:10)
  s <- suppressWarnings(draw(d, design_systematic(interval = 20), seed = 2,
                             weights = TRUE))
  expect_equal(nrow(s), 0L)
  expect_error(ht_mean(s, "y"), "no rows")
  expect_error(ht_total(s, "y"), "no rows")
})

# ---- a sample's rows are matched to the frame through `.row` ----------------

reorder_designs <- function() {
  list(
    stratified = design_stratified("site", n = 30),
    systematic = design_systematic(interval = 10),
    certainty  = design_certainty("w", 440, design_stratified("site", n = 20))
  )
}

estimate_summary <- function(s) {
  t <- ht_total(s, "y")
  m <- ht_mean(s, "y")
  b <- ht_total(s, "y", by = "site")
  list(total = t$total, se = t$se, df = t$df, deff = t$deff, method = t$method,
       mean_se = m$se, by_se = b$se[order(b$site)])
}

test_that("reordering a weighted sample leaves every estimate unchanged", {
  # The design's rows were recovered by position, so sorting a sample moved
  # each row onto another row's stratum, cluster or systematic position:
  # a stratified standard error fell from 3014 to 1332 on a sort by `y`.
  d <- est_pop()
  for (des in reorder_designs()) {
    s <- draw(d, des, seed = 3, weights = TRUE)
    want <- estimate_summary(s)
    expect_equal(estimate_summary(s[order(s$y), ]), want)
    expect_equal(estimate_summary(s[rev(seq_len(nrow(s))), ]), want)
  }
})

test_that("an in-place data.table sort of a sample changes nothing", {
  skip_if_not_installed("data.table")
  d <- est_pop()
  for (des in reorder_designs()) {
    for (s in list(draw(d, des, seed = 3, weights = TRUE),
                   draw(data.table::as.data.table(d), des, seed = 3,
                        weights = TRUE))) {
      want <- estimate_summary(s)
      drawn_order <- attr(s, "drawn_rows")
      data.table::setorderv(s, "y")
      expect_false(is.unsorted(s$y))
      expect_identical(attr(s, "drawn_rows"), drawn_order)
      expect_equal(estimate_summary(s), want)
    }
  }
})

test_that("an in-place data.table edit of the frame does not reach the sample", {
  # The sample held the caller's frame by reference: setkey() on it after
  # drawing doubled a cluster design's standard error.
  skip_if_not_installed("data.table")
  dt <- data.table::as.data.table(est_pop())
  for (des in list(design_cluster("cl", n_clusters = 6),
                   design_certainty("w", 440, design_simple(n = 20)))) {
    s <- draw(dt, des, seed = 3, weights = TRUE)
    want <- estimate_summary(s)
    data.table::setkeyv(dt, "w")
    data.table::set(dt, j = "site", value = "z")
    expect_equal(estimate_summary(s), want)
    data.table::setorderv(dt, "id")
  }
})

test_that("a subset, a bound or a merged sample is refused, not misread", {
  d <- est_pop()
  s <- draw(d, design_stratified("site", n = 30), seed = 3, weights = TRUE)
  expect_error(ht_total(head(s), "y"), "`sample` has 6 of the 30 rows",
               fixed = TRUE)
  expect_error(ht_total(s[s$site == "a", ], "y"), "use `by` for groups",
               fixed = TRUE)
  expect_error(ht_mean(rbind(s, s), "y"), "more than once")
  other <- draw(d, design_stratified("site", n = 30), seed = 4, weights = TRUE)
  other <- other[!other$.row %in% s$.row, ]
  expect_error(ht_total(rbind(s, other), "y"), "that were not in the sample")

  merged <- merge(s, data.frame(id = s$id, resp = 1), by = "id")
  expect_error(ht_total(merged, "y"), "does not carry its design")
  expect_error(ht_total(merged, "y"), "merge()", fixed = TRUE)
  expect_error(ht_total(s[, c("y", "site")], "y"), "column selection")
})

test_that("a damaged `.row` column is refused", {
  d <- est_pop()
  s <- draw(d, design_simple(n = 20), seed = 3, weights = TRUE)
  gone <- s
  gone$.row <- NULL
  expect_error(ht_total(gone, "y"), "has no `.row` column", fixed = TRUE)
  far <- s
  far$.row[1] <- 1000L
  expect_error(ht_total(far, "y"), "not rows of the 300-row frame")
  half <- s
  half$.row <- half$.row + 0.5
  expect_error(ht_total(half, "y"), "not rows of the 300-row frame")
  chr <- s
  chr$.row <- as.character(chr$.row)
  expect_error(ht_total(chr, "y"), "character `.row` column", fixed = TRUE)
})

test_that("a sample saved by drawn 0.2.0 still estimates, by position", {
  d <- est_pop()
  s <- draw(d, design_stratified("site", n = 30), seed = 3, weights = TRUE)
  old <- s
  old$.row <- NULL
  attr(old, "drawn_pi") <- NULL
  expect_equal(ht_total(old, "y")$se, ht_total(s, "y")$se)
  expect_error(ht_total(old[1:10, ], "y"), "has 10 rows but was drawn with 30")
})

test_that("a `.prob` changed or removed after the draw is refused", {
  d <- est_pop()
  s <- draw(d, design_simple(n = 20), seed = 3, weights = TRUE)

  gone <- s
  gone$.prob <- NULL
  expect_error(ht_total(gone, "y"), "has no `.prob` column", fixed = TRUE)
  missing_p <- s
  missing_p$.prob[1] <- NA
  expect_error(ht_total(missing_p, "y"),
               "not probabilities in (0, 1] in 1 row(s): NA", fixed = TRUE)
  zero <- s
  zero$.prob[2] <- 0
  expect_error(ht_mean(zero, "y"), "in 1 row(s): 0.", fixed = TRUE)
  chr <- s
  chr$.prob <- as.character(chr$.prob)
  expect_error(ht_total(chr, "y"), "character `.prob` column", fixed = TRUE)

  # Within (0, 1], but not what was drawn: the total would move while the
  # variance still followed the design.
  edited <- s
  edited$.prob[1:3] <- edited$.prob[1:3] * 1.1
  expect_error(ht_total(edited, "y"), "in 3 row(s) that differ", fixed = TRUE)
  expect_error(ht_mean(edited, "y"), "estimate from them with as_svydesign()",
               fixed = TRUE)
  expect_error(sample_summary(edited), "`.prob` was changed after draw()",
               fixed = TRUE)

  # A sample saved by drawn 0.2.0 records no probabilities to compare with,
  # but its `.prob` must still be a probability.
  old <- s
  old$.row <- NULL
  attr(old, "drawn_pi") <- NULL
  old$.prob[1] <- NA
  expect_error(ht_total(old, "y"), "in 1 row(s): NA", fixed = TRUE)
})

test_that("a `y` that is not finite is refused with a count", {
  d <- est_pop()
  s <- draw(d, design_simple(n = 20), seed = 3, weights = TRUE)
  s$y[c(2, 5)] <- c(Inf, NaN)
  expect_error(ht_total(s, "y"),
               "`y` is Inf or NaN in 2 row(s) of `sample`: 2, 5.", fixed = TRUE)
  expect_error(ht_mean(s, "y", by = "site"), "finite `y`")
  s$y[c(2, 5)] <- c(-Inf, 1)
  expect_error(ht_mean(s, "y"), "`y` is -Inf in 1 row(s)", fixed = TRUE)
})

test_that("an edited `.weight` warns, and the estimate keeps to `.prob`", {
  d <- est_pop()
  s <- draw(d, design_stratified("site", n = 30), seed = 1, weights = TRUE)
  want <- ht_total(s, "y")
  w <- s
  w$.weight <- w$.weight * 1.1
  expect_warning(got <- ht_total(w, "y"), class = "drawn_weight_edited")
  expect_equal(got$total, want$total)
  expect_equal(got$se, want$se)
  expect_warning(ht_mean(w, "y"), "differ from 1 / `.prob` in 30 row(s)",
                 fixed = TRUE)
  expect_warning(ht_total(w, "y", by = "site"), class = "drawn_weight_edited")
  expect_warning(ht_total(w, "y"), class = "drawn_warning")
  expect_no_warning(ht_total(s, "y"))
})

test_that("a design that gives the frame other probabilities now is refused", {
  # Editing the stored probabilities of rows that were not drawn stands in for
  # a sample drawn by another version of drawn, or in a session whose settings
  # change what the design computes.
  d <- est_pop()
  s <- draw(d, design_simple(n = 20), seed = 3, weights = TRUE)
  p <- attr(s, "drawn_pi")
  p[setdiff(seq_len(nrow(d)), s$.row)[1:2]] <- 0.5
  attr(s, "drawn_pi") <- p
  expect_error(ht_total(s, "y"), paste0(
    "gives 2 of its 300 frame rows a different inclusion probability in this ",
    "session or this version of drawn"), fixed = TRUE)
  expect_error(ht_mean(s, "y"), "draw the sample again")
})

test_that("dplyr groups are reported, not silently ignored", {
  d <- est_pop()
  s <- draw(d, design_stratified("site", n = 30), seed = 1, weights = TRUE)
  # The class and `groups` attribute dplyr::group_by() gives, without dplyr.
  g <- s
  groups <- data.frame(site = c("a", "b", "c"))
  groups$.rows <- unname(split(seq_len(nrow(s)), s$site))
  attr(g, "groups") <- groups
  class(g) <- c("grouped_df", "data.frame")
  expect_warning(m <- ht_mean(g, "y"), "Use `by = \"site\"`", fixed = TRUE)
  expect_equal(m$mean, ht_mean(s, "y")$mean)
  expect_warning(ht_total(g, "y"), class = "drawn_grouped_df")
  expect_no_warning(ht_total(g, "y", by = "site"))
})
