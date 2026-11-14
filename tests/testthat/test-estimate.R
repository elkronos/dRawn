est_pop <- function() {
  local_seed(404)
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
    est <- vapply(1:400, function(i) {
      s <- draw(d, des, seed = i, weights = TRUE)
      ht_mean(s, "y", variance = "none")$mean
    }, numeric(1))
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

test_that("equal probabilities over the rows reached part Hajek from HT", {
  # A reservoir capped at max_items is a simple random sample of the rows it
  # reaches, so its weights sum to 100, not 200, and the HT mean, which
  # divides by all 200, is half the Hajek mean.
  d <- data.frame(y = 1:200)
  expect_warning(
    s <- draw(d, design_reservoir(n = 20, max_items = 100), seed = 1,
              weights = TRUE),
    "max_items"
  )
  expect_equal(sum(s$.weight), 100)
  hajek <- ht_mean(s, "y", variance = "none")$mean
  expect_equal(ht_mean(s, "y", "ht", variance = "none")$mean, hajek / 2)
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
  s$when <- as.Date("2024-01-01") + seq_len(nrow(s))
  expect_error(ht_mean(s, "when"), "numeric")
  expect_error(ht_mean(s, "y", level = 1), "between 0 and 1")

  na_y <- s
  na_y$y[1] <- NA
  expect_error(ht_mean(na_y, "y"), "missing")
})

test_that("deff falls below 1 when stratification helps and rises when it hurts", {
  local_seed(88)
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
  # Anything else goes to survey's deff(), which has no design effect for it.
  skip_if_not_installed("survey")
  expect_warning(expect_null(deff(1:10)), "no design effect")
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
  local_seed(31)
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
  local_seed(5)
  pop <- data.frame(y = stats::rnorm(60, 10, 3), cl = rep(1:6, each = 10))
  r <- ht_total(draw(pop, design_cluster("cl", n_clusters = 1), seed = 1,
                     weights = TRUE), "y")
  expect_true(is.na(deff(r)))
  expect_true(is.na(drawn:::deff_value(-5, 1:10, rep(0.5, 10), 100, 10,
                                       "total")))
})

test_that("an empty sample gives a total of 0 and no Hajek mean, not a crash", {
  # A start past the end of a short frame selects nothing. The total of that
  # draw is 0, which keeps the estimator unbiased; the Hajek mean is 0 / 0.
  d <- data.frame(y = 1:10)
  s <- suppressWarnings(draw(d, design_systematic(interval = 20), seed = 2,
                             weights = TRUE))
  expect_equal(nrow(s), 0L)
  expect_error(ht_mean(s, "y"), "no rows")
  r <- ht_total(s, "y")
  expect_equal(r$total, 0)
  expect_true(is.na(r$se))
  expect_match(r$note, "The sample is empty")
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
  gap <- s
  gap$.row[2] <- NA
  expect_error(ht_total(gap, "y"), "frame it was drawn from: NA.",
               fixed = TRUE)
  gap$.row[5] <- 1000L
  expect_error(ht_mean(gap, "y"), "frame it was drawn from: 1000, NA.",
               fixed = TRUE)
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

test_that("a sample saved by drawn 0.2.0 has its design brought up to date", {
  # A systematic design saved then has no `n` field, and `design$n` matched
  # its `na_rm` by partial name, so the design read as `n = FALSE`.
  d <- est_pop()
  s <- draw(d, design_systematic(interval = 10, order_by = "w"), seed = 2,
            weights = TRUE)
  old <- s
  old$.row <- NULL
  attr(old, "drawn_pi") <- NULL
  des <- attr(s, "drawn_design")
  attr(old, "drawn_design") <- structure(
    unclass(des)[c("interval", "start", "order_by", "na_rm")],
    class = class(des))
  want <- ht_total(s, "y")
  got <- ht_total(old, "y")
  expect_equal(got$se, want$se)
  expect_equal(got$df, want$df)
  expect_equal(sample_summary(old)$N_reachable, nrow(d))
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
  many <- s
  many$.prob[1:7] <- -(1:7)
  expect_error(ht_total(many, "y"),
               "in 7 row(s): -1, -2, -3, -4, -5 and 2 more. `.prob` was",
               fixed = TRUE)
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

  # A sample saved by drawn 0.2.0 records no probabilities for the frame, so
  # its `.prob` must be a probability and the one the design gives its row.
  old <- s
  old$.row <- NULL
  attr(old, "drawn_pi") <- NULL
  old$.prob[1] <- NA
  expect_error(ht_total(old, "y"), "in 1 row(s): NA", fixed = TRUE)
  expect_error(sample_summary(old), "in 1 row(s): NA", fixed = TRUE)
  old$.prob <- NULL
  expect_error(sample_summary(old), "has no `.prob` column", fixed = TRUE)
  old <- edited
  old$.row <- NULL
  attr(old, "drawn_pi") <- NULL
  expect_error(ht_total(old, "y"), paste0(
    "gives 3 of its 20 rows a different inclusion probability in this ",
    "session or this version of drawn than its `.prob` column records"),
    fixed = TRUE)
  expect_error(ht_mean(old, "y"), "`.prob` changed after draw()",
               fixed = TRUE)
  expect_error(sample_summary(old), "matched to its frame rows by position",
               fixed = TRUE)
  skip_if_not_installed("survey")
  expect_error(as_svydesign(old), "than its `.prob` column records",
               fixed = TRUE)
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

# ---- proportions, one-sided bounds, design effects and reachable rows -------

test_that("whether y is a proportion is read from the frame, not the sample", {
  # The frame's column holds 3s as well, so this is the mean of a count. A
  # logit interval, which cannot pass 1, would never cover a mean above 1.
  pop <- data.frame(id = 1:100, kids = c(rep(0, 60), rep(1, 36), rep(3, 4)))
  pop$any <- as.numeric(pop$kids > 0)
  s <- draw(pop, design_simple(n = 10), seed = 1, weights = TRUE)
  expect_true(all(s$kids %in% 0:1))
  expect_equal(ht_mean(s, "kids")$ci_scale, "linear")
  expect_equal(ht_mean(s, "kids", proportion = TRUE)$ci_scale, "logit")

  # A frame column of 0s and 1s, and a logical vector, are proportions
  expect_equal(ht_mean(s, "any")$ci_scale, "logit")
  expect_equal(ht_mean(s, s$kids > 0)$ci_scale, "logit")
  expect_equal(ht_mean(s, "any", proportion = FALSE)$ci_scale, "linear")

  # A 0/1 column added after the draw, or a numeric vector, has no frame
  # column to read, so its sampled values decide and the note says so
  s$found <- as.numeric(s$kids > 0)
  expect_equal(ht_mean(s, "found")$ci_scale, "logit")
  expect_equal(ht_mean(s, s$found)$ci_scale, "logit")
  expect_equal(ht_mean(s, "found")$ci, ht_mean(s, "any")$ci)
  expect_equal(ht_mean(s, "found", proportion = TRUE)$ci,
               ht_mean(s, "any")$ci)
  expect_equal(ht_mean(s, "found", proportion = FALSE)$ci_scale, "linear")
  # A single number stands for a constant, not a proportion
  expect_equal(ht_mean(s, 1)$ci_scale, "linear")

  expect_error(ht_mean(s, "kids", proportion = NA),
               "`proportion` must be NULL, TRUE or FALSE, not NA.",
               fixed = TRUE)
  s$kids[1] <- 3
  expect_error(ht_mean(s, "kids", proportion = TRUE),
               "needs `y` to be 0 or 1 in every sampled row")
  expect_error(ht_mean(s, "kids", ci_method = "beta"),
               "Pass `proportion = TRUE`")
  expect_error(ht_mean(s, "any", ci_method = "exact"),
               "`ci_method` must be one of \"logit\" or \"beta\"")
})

test_that("a 0/1 column recorded in fieldwork keeps its logit interval", {
  # Deviations found in 2 of 60 sampled rows, recorded after the draw. With
  # no frame column to read, the sampled values decide, as they always did,
  # and the note says so. The symmetric interval ran below 0 and put the
  # upper bound at 0.078 against 0.120.
  pop <- data.frame(id = 1:600)
  s <- draw(pop, design_simple(n = 60), seed = 1, weights = TRUE)
  s$dev <- as.numeric(seq_len(60) %in% c(5, 40))
  r <- ht_mean(s, "dev")
  expect_equal(r$ci_scale, "logit")
  half <- stats::qt(0.975, 59) * r$se / (r$mean * (1 - r$mean))
  expect_equal(r$ci, stats::plogis(stats::qlogis(r$mean) + c(-1, 1) * half))
  expect_equal(round(r$ci, 4), c(0.0086, 0.1202))
  expect_match(r$note, "`dev` is not a column of the frame, so it was taken",
               fixed = TRUE)
  expect_output(print(r), "taken for a proportion")
  expect_match(ht_mean(s, s$dev)$note, "`y` is not a column of the frame",
               fixed = TRUE)
  s$half <- rep(c("a", "b"), 30)
  by_half <- ht_mean(s, "dev", by = "half")
  expect_equal(attr(by_half, "ci_scale"), "logit")
  expect_match(attr(by_half, "notes"), "`dev` is not a column of the frame",
               fixed = TRUE, all = FALSE)

  # proportion = FALSE gives the symmetric interval, with no note
  lin <- ht_mean(s, "dev", proportion = FALSE)
  expect_equal(lin$ci_scale, "linear")
  expect_lt(lin$ci[1], 0)
  expect_null(lin$note)
  # proportion = TRUE settles it, so there is nothing to note
  expect_null(ht_mean(s, "dev", proportion = TRUE)$note)

  # A fieldwork column of zeros gets the beta interval
  s$none <- 0
  z <- ht_mean(s, "none")
  expect_equal(z$ci, c(0, 1 - 0.025^(1 / 60)))
  expect_match(z$note, "not a column of the frame")
})

test_that("a proportion of 0 or 1 gets a beta interval and a note", {
  # No deviations in 60 sampled invoices: the variance is 0, and a logit
  # interval would be [0, 0]. Clopper-Pearson over 60 rows gives the bounds.
  pop <- data.frame(id = 1:240, deviation = FALSE)
  s <- draw(pop, design_simple(n = 60), seed = 1, weights = TRUE)
  r <- ht_mean(s, "deviation")
  expect_equal(r$mean, 0)
  expect_equal(r$ci_scale, "beta")
  expect_equal(r$ci, c(0, 1 - 0.025^(1 / 60)))
  expect_match(r$note, "beta interval, with the 60 sampled rows")
  expect_output(print(r), "beta, t, 59 df")
  expect_output(print(r), "would have no width")

  up <- ht_mean(s, "deviation", alternative = "less")
  expect_equal(up$ci, c(0, 1 - 0.05^(1 / 60)))
  expect_equal(round(up$ci[2], 4), 0.0487)
  expect_output(print(up), "95% upper bound  0.0487")

  s$deviation <- TRUE
  expect_equal(ht_mean(s, "deviation")$ci, c(0.025^(1 / 60), 1))

  # Rows taken with certainty carry no sampling error and are not counted
  pop$w <- c(rep(1000, 10), rep(1, 230))
  cs <- draw(pop, design_certainty("w", 1000, design_simple(n = 20)),
             seed = 1, weights = TRUE)
  rc <- ht_mean(cs, "deviation")
  expect_equal(rc$ci, c(0, 1 - 0.025^(1 / 20)))
  expect_match(rc$note, "20 sampled rows below certainty")

  # A census of zeros is exact
  cen <- draw(pop[1:20, ], design_simple(n = 20), seed = 1, weights = TRUE)
  expect_equal(ht_mean(cen, "deviation")$ci, c(0, 0))
})

test_that("a proportion of 0 at 0 degrees of freedom gets Clopper-Pearson", {
  # A Poisson sample of one row leaves 0 degrees of freedom. The beta
  # interval at an estimate of 0 needs no variance, but it was dropped, and
  # the note said both that there was no interval and that it was a beta one.
  pop <- data.frame(id = 1:50, w = 1:50, flag = FALSE)
  s <- draw(pop, design_weighted("w", n = 1, method = "poisson"), seed = 2,
            weights = TRUE)
  expect_equal(nrow(s), 1L)
  r <- ht_mean(s, "flag", estimator = "ht")
  expect_equal(r$df, 0)
  expect_equal(r$ci_scale, "beta")
  expect_equal(r$ci, c(0, 0.975))
  expect_match(r$note, "beta interval, with the 1 sampled row")
  expect_false(grepl("no t interval", r$note))
})

test_that("one-sided bounds are the ends of a two-sided interval", {
  pop <- data.frame(id = 1:200, site = rep(c("a", "b"), each = 100),
                    y = 20 + (1:200 * 37) %% 101)
  pop$flag <- pop$y > 60
  s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
  two <- ht_total(s, "y", level = 0.90)
  up <- ht_total(s, "y", alternative = "less")
  lo <- ht_total(s, "y", alternative = "greater")
  expect_equal(up$ci, c(-Inf, two$ci[2]))
  expect_equal(lo$ci, c(two$ci[1], Inf))
  expect_equal(up$alternative, "less")
  expect_output(print(up), "95% upper bound  ")
  expect_output(print(lo), "95% lower bound  ")

  # A proportion's bound has 0 or 1 at its open end, on either scale
  m2 <- ht_mean(s, "flag", level = 0.90)
  expect_equal(ht_mean(s, "flag", alternative = "less")$ci, c(0, m2$ci[2]))
  b2 <- ht_mean(s, "flag", level = 0.90, ci_method = "beta")
  expect_equal(ht_mean(s, "flag", alternative = "greater",
                       ci_method = "beta")$ci, c(b2$ci[1], 1))

  expect_error(ht_total(s, "y", alternative = "up"),
               "`alternative` must be one of")
  expect_error(ht_mean(s, "y", alternative = 1),
               "`alternative` must be one of")
})

test_that("a constant y has a variance of exactly 0 and no design effect", {
  # Unequal weights left the residuals a rounding error from 0, which gave a
  # standard error near 1e-16 and design effects from noise.
  pop <- data.frame(id = 1:240, w = 1 + (1:240 * 7) %% 19)
  s <- draw(pop, design_weighted("w", n = 30, method = "poisson"), seed = 8,
            weights = TRUE)
  s$k <- 5
  m <- ht_mean(s, "k")
  expect_identical(m$mean, 5)
  expect_identical(m$variance, 0)
  expect_equal(m$ci, c(5, 5))
  expect_true(is.na(m$deff))
  expect_match(m$note, "`y` is 5 in every sampled row")
  expect_output(print(m), "5 in every sampled row")
  tot <- ht_total(s, "k")
  expect_true(tot$se > 0)
  expect_true(is.na(tot$deff))

  # Data with a large offset is not constant, and keeps its design effect
  s$big <- 1e6 + (seq_len(nrow(s)) %% 3) / 2
  expect_true(is.finite(ht_mean(s, "big")$deff))
  expect_true(is.finite(ht_total(s, "big")$deff))
})

test_that("deff is measured against the rows the design can reach", {
  # An equal-probability sample of 20 from the first 100 of 400 rows. Against
  # all 400 rows the total's design effect came out 0.05.
  pop <- data.frame(y = 50 + (1:400 * 37) %% 23)
  expect_warning(
    s <- draw(pop, design_reservoir(n = 20, max_items = 100), seed = 1,
              weights = TRUE),
    "max_items")
  expect_equal(deff(ht_total(s, "y")), 1)
  expect_equal(deff(ht_mean(s, "y")), 1)
  # The Horvitz-Thompson mean is the total over all 400 rows
  expect_equal(deff(ht_mean(s, "y", estimator = "ht")),
               deff(ht_total(s, "y")))
})

test_that("an estimate says how many frame rows it cannot cover", {
  skew <- data.frame(id = 1:1000, g = rep(c("big", "rare"), c(990, 10)),
                     y = rep(c(1, 1000), c(990, 10)))
  expect_warning(
    s <- draw(skew, design_stratified("g", n = 20), seed = 1, weights = TRUE),
    class = "drawn_empty_stratum")
  r <- ht_total(s, "y")
  expect_equal(r$unreachable, 10)
  expect_output(print(r), "unreachable  10 frame rows at probability 0")
  expect_output(print(ht_total(s, "y", variance = "none")), "unreachable  10")
  expect_equal(ht_mean(s, "y")$unreachable, 10)

  all_rows <- draw(skew, design_simple(n = 20), seed = 1, weights = TRUE)
  expect_equal(ht_total(all_rows, "y")$unreachable, 0)
  out <- capture.output(print(ht_total(all_rows, "y")))
  expect_false(any(grepl("unreachable", out)))
})

# ---- the estimate as numbers ------------------------------------------------

methods_pop <- function() {
  id <- 1:300
  w <- 1 + (id * 13) %% 50
  w[c(5, 160, 260)] <- 500
  data.frame(id = id, site = rep(c("a", "b", "c"), times = c(150, 100, 50)),
             y = 20 + (id * 37) %% 101, hit = (id * 7) %% 5 == 0,
             none = FALSE, w = w)
}

test_that("an estimate gives coef(), vcov(), confint() and a data frame", {
  s <- draw(methods_pop(), design_stratified("site", n = 60), seed = 1,
            weights = TRUE)
  r <- ht_total(s, "y")
  expect_equal(r$variable, "y")
  expect_equal(coef(r), c(y = r$total))
  expect_equal(vcov(r), matrix(r$variance, 1, 1, dimnames = list("y", "y")))
  ci <- confint(r)
  expect_equal(unname(ci[1, ]), r$ci)
  expect_equal(colnames(ci), c("2.5 %", "97.5 %"))
  expect_equal(confint(r, "y"), ci)
  expect_equal(confint(r, 1), ci)
  expect_error(confint(r, level = 95), "strictly between 0 and 1")

  df <- as.data.frame(r)
  expect_equal(names(df), c("variable", "n", "total", "variance", "se",
                            "ci_lower", "ci_upper", "level", "df", "deff",
                            "method"))
  expect_equal(unlist(df[c("total", "se", "ci_lower", "ci_upper", "df")]),
               c(total = r$total, se = r$se, ci_lower = r$ci[1],
                 ci_upper = r$ci[2], df = r$df))

  # A vector `y` has no name, and the row says so.
  m <- ht_mean(s, s$y)
  expect_true("variable" %in% names(m))
  expect_null(m$variable)
  expect_equal(names(coef(m)), "mean")
  both <- rbind(as.data.frame(ht_mean(s, "y")), as.data.frame(m))
  expect_equal(both$variable, c("y", NA))
  expect_equal(both$mean, c(m$mean, m$mean))
})

test_that("confint() at another level is the estimate made at that level", {
  pop <- methods_pop()
  s <- draw(pop, design_stratified("site", n = 60), seed = 1, weights = TRUE)
  same <- function(fit) {
    expect_equal(unname(confint(fit(0.95), level = 0.8)[1, ]), fit(0.8)$ci)
  }
  same(function(l) ht_total(s, "y", level = l))
  same(function(l) ht_total(s, "y", level = l, alternative = "less"))
  same(function(l) ht_mean(s, "y", level = l, estimator = "ht"))
  # Proportions keep their scale: logit, beta, and beta at an estimate of 0.
  same(function(l) ht_mean(s, "hit", level = l))
  same(function(l) ht_mean(s, "hit", level = l, ci_method = "beta"))
  same(function(l) ht_mean(s, "none", level = l, alternative = "less"))
  expect_equal(colnames(confint(ht_total(s, "y", alternative = "less"),
                                level = 0.9)), c("0 %", "90 %"))

  # The beta interval's sample size leaves out the rows taken with certainty.
  cs <- draw(pop, design_certainty("w", 500, design_simple(n = 40)), seed = 2,
             weights = TRUE)
  b <- ht_mean(cs, "hit", ci_method = "beta")
  expect_equal(b$n_random, 40L)
  expect_equal(b$n, 43L)
  same(function(l) ht_mean(cs, "hit", level = l, ci_method = "beta"))
})

test_that("survey's SE() and cv() read an estimate", {
  skip_if_not_installed("survey")
  s <- draw(methods_pop(), design_stratified("site", n = 60), seed = 1,
            weights = TRUE)
  r <- ht_total(s, "y")
  expect_equal(survey::SE(r), c(y = r$se))
  expect_equal(survey::cv(r), c(y = r$se / r$total))
  m <- ht_mean(s, "y")
  expect_equal(survey::SE(m), c(y = m$se))
})

test_that("tidy() gives an estimate in broom's columns", {
  skip_if_not_installed("generics")
  s <- draw(methods_pop(), design_stratified("site", n = 60), seed = 1,
            weights = TRUE)
  r <- ht_total(s, "y")
  t <- generics::tidy(r)
  expect_equal(t, data.frame(term = "y", estimate = r$total,
                             std.error = r$se))
  t <- generics::tidy(r, conf.int = TRUE)
  expect_equal(c(t$conf.low, t$conf.high), r$ci)
  t <- generics::tidy(ht_mean(s, "hit"), conf.int = TRUE, conf.level = 0.8)
  expect_equal(c(t$conf.low, t$conf.high), ht_mean(s, "hit", level = 0.8)$ci)
})

test_that("deff() answers for survey's objects, and survey's deff() for these", {
  skip_if_not_installed("survey")
  s <- draw(methods_pop(), design_stratified("site", n = 60), seed = 1,
            weights = TRUE)
  des <- as_svydesign(s)
  # drawn attached after survey: this package's deff() is found first and
  # hands survey's objects to survey.
  m <- survey::svymean(~y, des, deff = TRUE)
  expect_equal(deff(m), survey::deff(m))
  by_site <- survey::svyby(~y, ~site, des, survey::svytotal, deff = TRUE)
  expect_equal(deff(by_site), survey::deff(by_site))
  expect_equal(deff(m, quietly = TRUE), survey::deff(m, quietly = TRUE))
  # survey attached after drawn: survey's deff() is found first and finds the
  # methods for the estimates made here.
  r <- ht_total(s, "y")
  expect_equal(survey::deff(r), r$deff)
  expect_equal(survey::deff(ht_mean(s, "y")), ht_mean(s, "y")$deff)
  expect_true(is.na(survey::deff(ht_total(s, "y", by = "site",
                                          contrast = c(a = 1, b = -1)))))
  expect_error(survey::deff(ht_total(s, "y", by = "site")),
               "no design effect")
})

test_that("a printed estimate and a by table name the method alike", {
  d <- data.frame(id = 1:200, g = rep(c("a", "b"), 100),
                  y = 20 + (1:200 * 37) %% 101)
  s <- draw(d, design_systematic(interval = 8), seed = 1, weights = TRUE)
  r <- ht_total(s, "y")
  b <- ht_total(s, "y", by = "g")
  expect_equal(r$method, "successive difference")
  expect_equal(unique(b$method), r$method)
  expect_output(print(r), "(successive difference, an approximation)",
                fixed = TRUE)
  expect_output(print(b), "successive difference")
  srs <- draw(d, design_simple(n = 20), seed = 1, weights = TRUE)
  expect_output(print(ht_total(srs, "y")), "(analytic)", fixed = TRUE)
})

test_that("ht_total() refuses a missing y and says why a respondent total is low", {
  d <- est_pop()
  s <- draw(d, design_stratified("site", n = 40), seed = 1, weights = TRUE)
  s$y[c(2, 7, 25)] <- NA
  msg <- tryCatch(ht_total(s, "y"), error = conditionMessage)
  expect_match(msg, "`y` is missing in 3 of the 40 sampled rows", fixed = TRUE)
  share <- sum(s$.weight[is.na(s$y)]) / sum(s$.weight)
  expect_match(msg, paste0(format(signif(100 * share, 2)), "% here"),
               fixed = TRUE)
  expect_match(msg, "ht_mean(na_rm = TRUE)", fixed = TRUE)
  expect_match(msg, "as_svydesign()", fixed = TRUE)
  expect_error(ht_mean(s, "y"), "Set `na_rm = TRUE`")
  expect_error(ht_mean(s, "y", na_rm = NA), "`na_rm` must be TRUE or FALSE")
})

test_that("ht_mean(na_rm = TRUE) is the mean over the rows with y present", {
  d <- est_pop()
  for (des in list(design_stratified("site", n = 40),
                   design_cluster("cl", n_clusters = 6))) {
    s <- draw(d, des, seed = 1, weights = TRUE)
    s$y[c(2, 7, 25)] <- NA
    m <- ht_mean(s, "y", na_rm = TRUE)
    r <- !is.na(s$y)
    expect_equal(m$mean, sum(s$y[r] * s$.weight[r]) / sum(s$.weight[r]))
    expect_identical(m$n, sum(r))
    expect_identical(m$n_missing, 3L)
    expect_true(is.na(m$deff))
    expect_match(m$note, "missing completely at random")
    expect_output(print(m), "3 sampled rows left out")
    # The same as the domain of the rows with y present
    s0 <- s
    s0$y0 <- ifelse(r, s$y, 0)
    s0$r <- r
    dom <- ht_mean(s0, "y0", by = "r")
    expect_equal(m$se, dom$se[dom$r])
    skip_if_not_installed("survey")
    sv <- survey::svymean(~y, as_svydesign(s), na.rm = TRUE)
    expect_equal(m$mean, unname(coef(sv)))
    expect_equal(m$se, as.numeric(survey::SE(sv)))
    by <- ht_mean(s, "y", by = "site", na_rm = TRUE)
    sb <- survey::svyby(~y, ~site, as_svydesign(s), survey::svymean,
                        na.rm = TRUE)
    expect_equal(by$se, unname(sb$se))
    expect_identical(by$n, as.integer(tapply(r, s$site, sum)))
  }
})

test_that("ht_mean(na_rm = TRUE) changes nothing when no y is missing", {
  d <- est_pop()
  s <- draw(d, design_stratified("site", n = 40), seed = 1, weights = TRUE)
  expect_identical(ht_mean(s, "y", na_rm = TRUE), ht_mean(s, "y"))
  expect_identical(ht_mean(s, "y", by = "site", na_rm = TRUE),
                   ht_mean(s, "y", by = "site"))
  expect_identical(ht_mean(s, "y")$n_missing, 0L)
})

test_that("ht_mean(na_rm = TRUE) refuses what it cannot estimate", {
  d <- est_pop()
  s <- draw(d, design_simple(n = 20), seed = 1, weights = TRUE)
  s$y[1] <- NA
  expect_error(ht_mean(s, "y", estimator = "ht", na_rm = TRUE),
               "number of frame rows the mean is over is unknown")
  s$y[] <- NA
  expect_error(ht_mean(s, "y", na_rm = TRUE), "missing in every one of the 20")
  # One row with y present has no standard error
  s$y[1] <- 5
  one <- ht_mean(s, "y", na_rm = TRUE)
  expect_true(is.na(one$se))
  expect_match(one$note, "The domain of rows with `y` present has one")
})

# ---- ratios -----------------------------------------------------------------

ratio_pop <- function() {
  id <- 1:300
  book <- 20 + (id * 37) %% 400
  data.frame(id = id, site = rep(c("a", "b", "c"), times = c(150, 100, 50)),
             cl = rep(paste0("c", 1:30), each = 10), book = book,
             audited = book * (0.9 + ((id * 7) %% 11) / 50))
}

test_that("a ratio and a ratio total match survey::svyratio()", {
  skip_if_not_installed("survey")
  d <- ratio_pop()
  for (des in list(design_stratified("site", n = 45),
                   design_cluster("cl", n_clusters = 6))) {
    s <- draw(d, des, seed = 3, weights = TRUE)
    sv <- survey::svyratio(~audited, ~book, as_svydesign(s))
    r <- ht_ratio(s, "audited", "book")
    expect_s3_class(r, "drawn_ratio")
    expect_equal(r$ratio, unname(coef(sv)))
    expect_equal(r$se, unname(survey::SE(sv)))
    expect_equal(r$ratio, r$total_y / r$total_x)
    expect_equal(r$df, ht_total(s, "audited")$df)

    t <- ht_total(s, "audited", ratio = "book")
    p <- predict(sv, total = sum(d$book))
    expect_s3_class(t, "drawn_ht")
    expect_equal(t$total, as.numeric(p$total))
    expect_equal(t$se, as.numeric(p$se))
    expect_equal(t$ratio, r$ratio)
    expect_equal(t$frame_total, sum(d$book))
  }
  # The design effect is measured as svyratio(deff = TRUE) measures it.
  s <- draw(d, design_stratified("site", n = 45), seed = 3, weights = TRUE)
  sv <- survey::svyratio(~audited, ~book, as_svydesign(s), deff = TRUE)
  expect_equal(ht_ratio(s, "audited", "book")$deff,
               as.numeric(attr(sv, "deff")))
})

test_that("the ratio total uses the frame total over the rows reached", {
  d <- ratio_pop()
  d$book[d$site == "c"] <- NA
  # Site c is outside the design, so its missing `book` plays no part.
  expect_warning(
    s <- draw(d, design_stratified("site", n = 40,
                                   allocation = c(a = 25, b = 15, c = 0)),
              seed = 3, weights = TRUE),
    class = "drawn_empty_stratum")
  t <- ht_total(s, "audited", ratio = "book")
  expect_equal(t$frame_total, sum(d$book[d$site != "c"]))
  expect_equal(t$frame_total, unname(frame_totals(s, ~book)["book"]))
  # Close to proportional to `book`, so far more precise than the HT total.
  expect_lt(t$se, ht_total(s, "audited")$se / 3)
  expect_output(print(t), "Ratio estimate of the total")
  expect_output(print(t), "`book` over the frame", fixed = TRUE)
})

test_that("ht_ratio() and ht_total(ratio =) refuse what they cannot estimate", {
  d <- ratio_pop()
  s <- draw(d, design_stratified("site", n = 45), seed = 3, weights = TRUE)
  s$gap <- s$book
  s$gap[3] <- NA
  expect_error(ht_ratio(s, "audited", "gap"), "`x` is missing in 1 of the 45")
  s$zero <- 0
  expect_error(ht_ratio(s, "audited", "zero"),
               "The estimated total of `zero`, `sum(x / pi)` over the sampled",
               fixed = TRUE)
  expect_error(ht_ratio(s, "audited", "site"),
               "`x` names `site`, a character column")
  expect_error(ht_ratio(s, "audited", c("book", "id")),
               "`x` must be one column name or a numeric vector")
  expect_error(ht_total(s, "audited", ratio = "zero"),
               "`ratio` names `zero`, which is not a column of the frame")
  expect_error(ht_total(s, "audited", ratio = "site"),
               "a character column of the frame")
  expect_error(ht_total(s, "audited", ratio = c("book", "id")),
               "`ratio` must be the name of one column")
  expect_error(ht_total(s, "audited", ratio = "book", by = "site"),
               "`ratio` gives one ratio estimator of the total")

  d$part <- d$book
  d$part[7] <- NA
  s <- draw(d, design_stratified("site", n = 45), seed = 3, weights = TRUE)
  expect_error(ht_total(s, "audited", ratio = "part"),
               "`part` is missing in 1 of the 300 frame rows")
  small <- compact_sample(s)
  expect_error(ht_total(small, "audited", ratio = "book"),
               "compact_sample() dropped", fixed = TRUE)
})

test_that("a ratio over one sampling unit has no se, an exact one a note", {
  d <- ratio_pop()
  s <- draw(d, design_cluster("cl", n_clusters = 1), seed = 1,
            weights = TRUE)
  r <- ht_ratio(s, "audited", "book")
  expect_true(is.na(r$se))
  expect_match(r$note, "A ratio's linearised residuals sum to zero")
  expect_true(is.na(ht_total(s, "audited", ratio = "book")$se))

  # `y` proportional to `x`: residuals of exactly 0 rather than rounding
  # error, and no design effect.
  s <- draw(d, design_simple(n = 30), seed = 1, weights = TRUE)
  s$part <- 0.3 * s$book
  r <- ht_ratio(s, "part", "book")
  expect_identical(r$variance, 0)
  expect_true(is.na(r$deff))
  expect_match(r$note, "`y / x` is 0.3 in every sampled row", fixed = TRUE)
})

test_that("a ratio has coef(), vcov(), confint(), as.data.frame() and tidy()", {
  s <- draw(ratio_pop(), design_stratified("site", n = 45), seed = 3,
            weights = TRUE)
  r <- ht_ratio(s, "audited", "book")
  expect_equal(coef(r), c("audited/book" = r$ratio))
  expect_equal(vcov(r), matrix(r$variance, 1, 1,
                               dimnames = list("audited/book",
                                               "audited/book")))
  expect_equal(unname(confint(r)[1, ]), r$ci)
  expect_equal(unname(confint(r, level = 0.8)[1, ]),
               ht_ratio(s, "audited", "book", level = 0.8)$ci)
  df <- as.data.frame(r)
  expect_equal(df$x, "book")
  expect_equal(df$ratio, r$ratio)
  expect_equal(df$se, r$se)
  expect_equal(deff(r), r$deff)
  expect_output(print(r), "Ratio of Horvitz-Thompson totals")
  expect_output(print(r), "`audited` / `book`", fixed = TRUE)

  # Vectors have no names.
  v <- ht_ratio(s, s$audited, s$book)
  expect_equal(names(coef(v)), "y/x")
  expect_equal(v$se, r$se)
  expect_null(v$x)

  skip_if_not_installed("generics")
  expect_equal(generics::tidy(r)$term, "audited/book")
  skip_if_not_installed("survey")
  expect_equal(survey::deff(r), r$deff)
})

# ---- several variables ------------------------------------------------------

test_that("several y give a table whose rows are their own estimates", {
  s <- draw(est_pop(), design_stratified("site", n = 40), seed = 1,
            weights = TRUE)
  s$hit <- s$y > 200
  b <- ht_total(s, c("y", "w"))
  expect_s3_class(b, "drawn_by")
  expect_equal(names(b), c("variable", "n", "total", "variance", "se",
                           "ci_lower", "ci_upper", "method"))
  expect_equal(b$variable, c("y", "w"))
  m <- ht_mean(s, c("y", "hit"))
  for (i in 1:2) {
    one <- ht_total(s, b$variable[i])
    expect_equal(b$total[i], one$total)
    expect_equal(b$se[i], one$se)
    expect_equal(c(b$ci_lower[i], b$ci_upper[i]), one$ci)
    one <- ht_mean(s, m$variable[i])
    expect_equal(m$mean[i], one$mean)
    expect_equal(c(m$ci_lower[i], m$ci_upper[i]), one$ci)
  }
  expect_equal(coef(b), c(y = b$total[1], w = b$total[2]))
  V <- vcov(b)
  expect_equal(dimnames(V), list(c("y", "w"), c("y", "w")))
  expect_equal(unname(diag(V)), b$variance)
  expect_equal(V[1, 2], V[2, 1])
  # The proportion keeps its logit interval at another level.
  expect_equal(unname(confint(m, level = 0.8)[2, ]),
               ht_mean(s, "hit", level = 0.8)$ci)
  expect_output(print(m), "Hajek mean by variable")
  expect_output(print(m), "logit for proportions")
  expect_error(deff(b), "several variables")
  skip_if_not_installed("generics")
  expect_equal(generics::tidy(b)$term, c("y", "w"))
})

test_that("the covariances of several estimates are survey's", {
  skip_if_not_installed("survey")
  d <- est_pop()
  d$z <- (d$id * 17) %% 23
  for (des in list(design_stratified("site", n = 40),
                   design_cluster("cl", n_clusters = 6))) {
    s <- draw(d, des, seed = 1, weights = TRUE)
    svd <- as_svydesign(s)
    b <- ht_total(s, c("y", "z"))
    expect_equal(unname(vcov(b)),
                 unname(matrix(vcov(survey::svytotal(~y + z, svd)), 2)))
    m <- ht_mean(s, c("y", "z"))
    expect_equal(unname(vcov(m)),
                 unname(matrix(vcov(survey::svymean(~y + z, svd)), 2)))
  }
})

test_that("several means with y missing are each over their own rows", {
  s <- draw(est_pop(), design_stratified("site", n = 40), seed = 1,
            weights = TRUE)
  s$y[c(2, 7)] <- NA
  m <- ht_mean(s, c("y", "w"), na_rm = TRUE)
  expect_equal(m$mean, c(ht_mean(s, "y", na_rm = TRUE)$mean,
                         ht_mean(s, "w")$mean))
  expect_equal(m$se, c(ht_mean(s, "y", na_rm = TRUE)$se, ht_mean(s, "w")$se))
  expect_identical(m$n, c(38L, 40L))
  expect_match(attr(m, "notes"), "^`y`: 2 of the 40 sampled rows",
               all = FALSE)
  expect_error(ht_total(s, c("w", "y")),
               "`y`, named in `y`, is missing in 2 of the 40", fixed = TRUE)

  # Each covariance is that of the two means' linearised residuals.
  skip_if_not_installed("survey")
  r <- !is.na(s$y)
  s$za <- 0
  s$za[r] <- (s$y[r] - m$mean[1]) / sum(s$.weight[r])
  s$zb <- (s$w - m$mean[2]) / sum(s$.weight)
  V <- vcov(survey::svytotal(~za + zb, as_svydesign(s)))
  expect_equal(vcov(m)[1, 2], V[1, 2])
})

test_that("a y that is not a column name or a numeric vector says so", {
  s <- draw(est_pop(), design_simple(n = 20), seed = 1, weights = TRUE)
  forms <- paste0("`y` must be one column name, several column names, or a ",
                  "numeric vector as long as `sample`.")
  expect_error(ht_total(s, letters[1:20]), forms, fixed = TRUE)
  expect_error(ht_total(s, letters[1:20]), "none of them a column")
  expect_error(ht_mean(s, c("y", NA)), "some of them missing or empty")
  expect_error(ht_total(s, 1:3),
               "an integer vector of length 3, and `sample` has 20 rows")
  expect_error(ht_total(s, list(1)), forms, fixed = TRUE)
  expect_error(ht_mean(s, c("y", "site")),
               "`y` names `site`, a character column, among several")
  expect_error(ht_total(s, c("y", "y")), "`y` names `y` more than once")
  expect_error(ht_total(s, c("y", "nope")), "`sample` has no column `nope`")
  expect_error(ht_total(s, c("y", "w"), by = "site",
                        contrast = c(a = 1, b = -1)),
               "`contrast` weighs the domains of one variable")
})

test_that("a factor y gives the share or count of each level, as survey's", {
  pop <- est_pop()
  pop$band <- factor(ifelse(pop$y < 150, "low",
                            ifelse(pop$y < 300, "mid", "high")),
                     levels = c("low", "mid", "high", "none"))
  pop$band_chr <- as.character(pop$band)
  s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
  m <- ht_mean(s, "band")
  expect_s3_class(m, "drawn_by")
  expect_equal(names(m), c("level", "n", "mean", "variance", "se",
                           "ci_lower", "ci_upper", "method"))
  expect_equal(levels(m$level), levels(pop$band))
  expect_equal(as.character(m$level), levels(pop$band))
  expect_equal(sum(m$mean), 1)
  expect_equal(attr(m, "ci_scale"), "logit")
  # Each row is the mean of the level's indicator.
  one <- ht_mean(s, s$band == "mid")
  expect_equal(m$mean[2], one$mean)
  expect_equal(c(m$ci_lower[2], m$ci_upper[2]), one$ci)
  # A level no sampled row has is a share of 0, with a beta interval.
  expect_equal(m$mean[4], 0)
  expect_gt(m$ci_upper[4], 0)
  expect_equal(ht_mean(s, s$band)$mean, m$mean)
  # A character column's levels are its values, sorted.
  t <- ht_total(s, "band_chr")
  expect_equal(t$level, c("high", "low", "mid"))
  expect_equal(sum(t$total), sum(1 / s$.prob))
  expect_output(print(m), "Hajek mean of `band` by level")
  expect_output(print(ht_total(s, "band", by = "site")),
                "total of `band` by level and domain")
  expect_error(deff(m), "levels of `y`")
  expect_error(ht_mean(s, "band", by = "site", contrast = c(a = 1, b = -1)),
               "Pass the indicator of one level as `y`")
  s$level <- s$site
  expect_error(ht_mean(s, "band", by = "level"), "uses for its own column")

  skip_if_not_installed("survey")
  svd <- as_svydesign(s)
  sm <- survey::svymean(~band, svd)
  expect_equal(unname(coef(m)), unname(coef(sm)))
  expect_equal(unname(vcov(m)), unname(matrix(vcov(sm), 4)))
  st <- survey::svytotal(~band_chr, svd)
  expect_equal(unname(coef(t)), unname(coef(st)))
  expect_equal(unname(vcov(t)), unname(matrix(vcov(st), 3)))
  # By domain, each level's share in each site.
  b <- ht_mean(s, "band", by = "site")
  sb <- survey::svyby(~band, ~site, svd, survey::svymean)
  at <- paste0(b$site, ":band", b$level)
  expect_equal(b$mean, unname(coef(sb)[at]))
  # survey's standard errors come a column per level, a row per site.
  expect_equal(b$se, c(as.matrix(survey::SE(sb))))
})

test_that("a single number as y stands for that value in every row", {
  s <- draw(est_pop(), design_weighted("w", n = 40, method = "poisson"),
            seed = 1, weights = TRUE)
  ones <- rep(1, nrow(s))
  n_hat <- ht_total(s, 1)
  expect_equal(n_hat$total, sum(1 / s$.prob))
  expect_equal(n_hat$se, ht_total(s, ones)$se)
  expect_equal(ht_total(s, 1, by = "site")$total,
               ht_total(s, ones, by = "site")$total)
  expect_equal(ht_ratio(s, "y", 1)$ratio, ht_mean(s, "y")$mean)
  # A single string is still a column name.
  expect_error(ht_total(s, "nope"), "`sample` has no column `nope`")
  expect_error(ht_total(s, letters[seq_len(nrow(s))]), "pass it as a factor")
})

test_that("ht_quantile() gives survey's quantiles and Woodruff intervals", {
  skip_if_not_installed("survey")
  pop <- est_pop()
  probs <- c(0.25, 0.5, 0.75)
  for (des in list(design_stratified("site", n = 40),
                   design_cluster("cl", n_clusters = 6),
                   design_multistage("cl", n_clusters = 8, n = 40))) {
    s <- draw(pop, des, seed = 2, weights = TRUE)
    q <- ht_quantile(s, "y", probs = probs)
    sq <- survey::svyquantile(~y, as_svydesign(s), probs, ci = TRUE,
                              interval.type = "mean")$y
    expect_equal(q$quantile, unname(sq[, 1]))
    expect_equal(cbind(q$ci_lower, q$ci_upper), unname(sq[, 2:3]))
    expect_equal(q$se, unname(sq[, 4]))
    expect_equal(attr(q, "df"), survey::degf(as_svydesign(s)))
    med <- ht_quantile(s, "y")
    expect_equal(c(med$quantile, med$ci, med$se),
                 unname(c(sq[2, 1], sq[2, 2:3], sq[2, 4])))
  }
  # By domain, with the normal interval survey also gives each domain.
  s <- draw(pop, design_stratified("site", n = 60), seed = 2, weights = TRUE)
  b <- ht_quantile(s, "y", probs = c(0.5, 0.75), by = "site", df = Inf)
  sb <- survey::svyby(~y, ~site, as_svydesign(s), survey::svyquantile,
                      quantiles = c(0.5, 0.75), ci = TRUE, df = Inf)
  expect_equal(b$quantile, c(t(as.matrix(sb[c("y.0.5", "y.0.75")]))),
               ignore_attr = TRUE)
  expect_equal(b$se, c(t(as.matrix(sb[c("se.y.0.5", "se.y.0.75")]))),
               ignore_attr = TRUE)
})

test_that("ht_quantile() follows the quantile rule and the indicator's mean", {
  pop <- est_pop()
  # With equal weights the rule is quantile(type = 1).
  sr <- draw(pop, design_simple(n = 30), seed = 3, weights = TRUE)
  p <- c(0, 0.1, 1 / 3, 0.5, 0.9, 1)
  q <- ht_quantile(sr, "y", p, variance = "none")
  expect_equal(q$quantile, unname(stats::quantile(sr$y, p, type = 1)))
  # F(q) and its standard error are the Hajek mean of y <= q, for a design
  # survey has no estimator for.
  sy <- draw(pop, design_systematic(interval = 7), seed = 1, weights = TRUE)
  med <- ht_quantile(sy, "y")
  ind <- ht_mean(sy, sy$y <= med$quantile, proportion = FALSE)
  expect_equal(med$cdf, ind$mean)
  expect_equal(med$cdf_se, ind$se)
  expect_equal(med$method, "successive difference")
  expect_gte(med$cdf, 0.5)
  expect_lte(med$ci[1], med$quantile)
  expect_gte(med$ci[2], med$quantile)
})

test_that("a quantile has coef(), vcov(), confint(), a frame and tidy()", {
  s <- draw(est_pop(), design_stratified("site", n = 60), seed = 1,
            weights = TRUE)
  q <- ht_quantile(s, "y")
  expect_s3_class(q, "drawn_quantile")
  expect_equal(coef(q), c(y = q$quantile))
  expect_equal(vcov(q), matrix(q$se^2, 1, 1, dimnames = list("y", "y")))
  expect_equal(unname(confint(q)[1, ]), q$ci)
  expect_equal(unname(confint(q, level = 0.8)[1, ]),
               ht_quantile(s, "y", level = 0.8)$ci)
  expect_equal(names(as.data.frame(q)),
               c("variable", "prob", "n", "quantile", "variance", "se",
                 "ci_lower", "ci_upper", "level", "df", "deff", "method"))
  expect_true(is.na(deff(q)))
  expect_output(print(q), "Quantile of `y`")
  expect_output(print(q), "Woodruff, t, 57 df", fixed = TRUE)

  b <- ht_quantile(s, "y", c(0.25, 0.5), by = "site")
  expect_equal(names(b), c("site", "prob", "n", "quantile", "variance", "se",
                           "ci_lower", "ci_upper", "method"))
  expect_equal(names(coef(b))[1:2], c("a . 0.25", "a . 0.5"))
  at <- ht_quantile(s, "y", c(0.25, 0.5), by = "site", level = 0.8)
  expect_equal(unname(confint(b, level = 0.8)),
               cbind(at$ci_lower, at$ci_upper))
  # Selected rows keep their own intervals at another level.
  expect_equal(unname(confint(b[3:4, ], level = 0.8)),
               cbind(at$ci_lower, at$ci_upper)[3:4, ])
  expect_equal(unname(diag(vcov(b))), b$variance)
  expect_output(print(b), "Quantiles of `y` by domain")
  expect_output(print(ht_quantile(s, "y", c(0.25, 0.75))),
                "Quantiles of `y`  (stratified", fixed = TRUE)
  expect_error(deff(b), "table of quantiles")
  skip_if_not_installed("generics")
  expect_equal(generics::tidy(q, conf.int = TRUE)$conf.low, q$ci[1])
  expect_equal(names(generics::tidy(b)),
               c("site", "prob", "term", "estimate", "std.error"))
})

test_that("ht_quantile() declines what it cannot measure", {
  pop <- est_pop()
  s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
  # At the largest value F(q) is 1, with a variance of 0 whatever y does.
  top <- ht_quantile(s, "y", 0.99)
  expect_equal(top$quantile, max(s$y))
  expect_true(is.na(top$se))
  expect_equal(top$method, "none")
  expect_match(top$note, "largest sampled value")
  b <- ht_quantile(s, "y", c(0.5, 0.99), by = "site")
  expect_equal(is.na(b$se), rep(c(FALSE, TRUE), 3))
  expect_match(attr(b, "notes"), "`a . 0.99`, `b . 0.99`", all = FALSE)
  # A census is exact there.
  cen <- draw(pop, design_simple(n = 300), seed = 1, weights = TRUE)
  ex <- ht_quantile(cen, "y", 1)
  expect_equal(ex$se, 0)
  expect_equal(ex$ci, rep(max(pop$y), 2))
  # An interval for F(q) below 0 bounds the quantile from above only.
  low <- ht_quantile(s, "y", 0.01)
  expect_true(is.na(low$ci[1]))
  expect_false(is.na(low$ci[2]))
  expect_true(is.na(low$se))
  expect_match(low$note, "bounds the quantile on one side only")
  # One sampled row in a domain.
  s$one <- c("x", rep("y", nrow(s) - 1L))
  expect_match(attr(ht_quantile(s, "y", by = "one"), "notes"),
               "A quantile's linearised residuals sum to zero",
               all = FALSE)
})

test_that("ht_quantile() refuses what it cannot read", {
  s <- draw(est_pop(), design_stratified("site", n = 40), seed = 1,
            weights = TRUE)
  s$gap <- s$y
  s$gap[3] <- NA
  expect_error(ht_quantile(s, "gap"), "`y` is missing in 1 of the 40")
  expect_error(ht_quantile(s, "site"), "a character column")
  expect_error(ht_quantile(s, c("y", "w")), "one column name or a numeric")
  expect_error(ht_quantile(s, "y", probs = 1.5),
               "`probs` must be one or more numbers from 0 to 1")
  expect_error(ht_quantile(s, "y", probs = c(0.5, 0.5)), "more than once")
  s$prob <- 1
  expect_error(ht_quantile(s, "y", by = "prob"), "uses for its own column")
  d <- data.frame(y = 1:10)
  e <- suppressWarnings(draw(d, design_systematic(interval = 20), seed = 2,
                             weights = TRUE))
  expect_error(ht_quantile(e, "y"), "no rows, so there is no quantile")
})
