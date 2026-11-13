var_pop <- function() {
  set.seed(7)
  N <- 240
  d <- data.frame(
    id   = 1:N,
    site = rep(c("a", "b", "c"), times = c(120, 80, 40)),
    cl   = rep(paste0("c", 1:24), each = 10),
    y    = round(stats::runif(N, 5, 200))
  )
  d$w <- d$y + 50
  d$ts <- rep(seq(as.POSIXct("2024-01-01", tz = "UTC"), by = "day",
                  length.out = 12), each = 20)
  d
}

test_that("joint probabilities match the textbook forms", {
  d <- data.frame(id = 1:10)
  m <- joint_prob(d, design_simple(n = 4))
  expect_equal(dim(m), c(10L, 10L))
  expect_equal(diag(m), rep(0.4, 10))
  expect_equal(m[1, 2], 4 * 3 / (10 * 9))

  # Poisson: rows are independent, so the joint is the product
  w <- data.frame(id = 1:8, w = c(1, 1, 2, 2, 3, 3, 4, 4))
  pd <- design_weighted("w", n = 3, method = "poisson")
  pj <- joint_prob(w, pd)
  pi <- inclusion_prob(w, pd)
  expect_equal(pj[1, 2], pi[1] * pi[2])
})

test_that("rows selects a submatrix rather than the whole square", {
  d <- data.frame(id = 1:1000)
  m <- joint_prob(d, design_simple(n = 50), rows = c(3, 17, 400))
  expect_equal(dim(m), c(3L, 3L))
  expect_equal(diag(m), rep(50 / 1000, 3))
})

test_that("cluster joints differ within and between clusters", {
  d <- data.frame(id = 1:20, cl = rep(letters[1:5], each = 4))
  m <- joint_prob(d, design_cluster("cl", n_clusters = 2))
  expect_equal(m[1, 2], 2 / 5)            # same cluster: that cluster selected
  expect_equal(m[1, 5], 2 * 1 / (5 * 4))  # different clusters: both selected
})

test_that("systematic joints are zero for rows that cannot co-occur", {
  d <- data.frame(id = 1:20)
  m <- joint_prob(d, design_systematic(interval = 5))
  expect_equal(m[1, 6], 1 / 5)   # same residue class
  expect_equal(m[1, 2], 0)       # different residue class: impossible together
})

test_that("designs without a closed form refuse", {
  d <- var_pop()
  expect_error(joint_prob(d, design_weighted("w", n = 10, method = "systematic")),
               "no closed-form joint")
  expect_error(joint_prob(d, design_weighted("w", n = 10)), "no closed-form")
  expect_error(joint_prob(d, design_cluster("cl", n_clusters = 3, balanced = TRUE)),
               "no closed-form")
  expect_error(joint_prob(d, design_bootstrap(n_replicates = 5)),
               "not a probability sample")
  expect_error(joint_prob(d, design_multistage("cl", n_clusters = 4, n = 10)),
               "does not divide|not constant")
})

test_that("ht_total recovers the total and needs weights = TRUE", {
  d <- var_pop()
  s <- draw(d, design_stratified("site", n = 40), seed = 1, weights = TRUE)
  r <- ht_total(s, "y")
  expect_s3_class(r, "drawn_ht")
  expect_equal(r$n, 40)
  expect_true(is.finite(r$se))
  expect_true(r$ci[1] < r$total && r$total < r$ci[2])

  plain <- draw(d, design_stratified("site", n = 40), seed = 1)
  expect_error(ht_total(plain, "y"), "does not carry its design")
})

test_that("ht_total validates y and level", {
  d <- var_pop()
  s <- draw(d, design_simple(n = 30), seed = 1, weights = TRUE)
  expect_error(ht_total(s, "nope"), "no column")
  expect_error(ht_total(s, letters[1:30]), "numeric column name")
  expect_error(ht_total(s, "y", level = 1.5), "between 0 and 1")
})

test_that("the estimated variance matches the estimator's true variance", {
  d <- var_pop()
  designs <- list(
    simple     = design_simple(n = 30),
    stratified = design_stratified("site", n = 30),
    cluster    = design_cluster("cl", n_clusters = 6),
    multistage = design_multistage("cl", n_clusters = 8, n = 32),
    poisson    = design_weighted("w", n = 30, method = "poisson")
  )
  R <- 700
  for (nm in names(designs)) {
    est <- numeric(R); vr <- numeric(R)
    for (i in seq_len(R)) {
      s <- draw(d, designs[[nm]], seed = i, weights = TRUE)
      if (nrow(s) < 2) { est[i] <- NA; vr[i] <- NA; next }
      r <- ht_total(s, "y")
      est[i] <- r$total; vr[i] <- r$variance
    }
    ok <- !is.na(est) & !is.na(vr)
    ratio <- mean(vr[ok]) / stats::var(est[ok])
    expect_gt(ratio, 0.75, label = paste(nm, "variance ratio"))
    expect_lt(ratio, 1.35, label = paste(nm, "variance ratio"))
  }
})

test_that("systematic sampling reports a successive-difference variance and says so", {
  d <- var_pop()
  s <- draw(d, design_systematic(interval = 8), seed = 1, weights = TRUE)
  r <- ht_total(s, "y")
  expect_true(is.finite(r$total))
  expect_equal(r$method, "successive difference")
  expect_true(r$variance > 0)
  expect_match(r$note, "no design-unbiased variance")
  expect_true(all(is.finite(r$ci)))
})

test_that("the print method shows the interval, or the reason there isn't one", {
  d <- var_pop()
  out <- capture.output(print(ht_total(draw(d, design_simple(n = 30), seed = 1,
                                            weights = TRUE), "y")))
  expect_true(any(grepl("Horvitz-Thompson total", out)))
  expect_true(any(grepl("95% CI", out)))

  out2 <- capture.output(print(ht_total(draw(d, design_cluster("cl",
                                                             n_clusters = 1),
                                             seed = 1, weights = TRUE), "y")))
  expect_true(any(grepl("se\\s+NA", out2)))

  # An approximation is named in the print, never passed off as exact
  out3 <- capture.output(print(ht_total(draw(d, design_systematic(interval = 8),
                                             seed = 1, weights = TRUE), "y")))
  expect_true(any(grepl("successive-difference approximation", out3)))
})

test_that("the jackknife supplies a variance where the analytic form cannot", {
  set.seed(3); N <- 240
  pop <- data.frame(id = 1:N, cl = rep(paste0("c", 1:24), each = 10),
                    y = round(stats::runif(N, 5, 200)))
  pop$w <- pop$y + 50

  # One row per cluster: the exact two-stage estimator declines
  r <- ht_total(draw(pop, design_multistage("cl", n_clusters = 6, n = 6),
                     seed = 1, weights = TRUE), "y")
  expect_equal(r$method, "jackknife")
  expect_true(is.finite(r$se))
  expect_match(r$note, "jackknife was used")
})

test_that("jackknife and analytic agree where both exist", {
  set.seed(3); N <- 240
  pop <- data.frame(id = 1:N, cl = rep(paste0("c", 1:24), each = 10),
                    y = round(stats::runif(N, 5, 200)))
  for (d in list(design_simple(n = 30), design_cluster("cl", n_clusters = 6))) {
    s <- draw(pop, d, seed = 1, weights = TRUE)
    a <- ht_total(s, "y", variance = "analytic")$variance
    j <- ht_total(s, "y", variance = "jackknife")$variance
    expect_equal(j, a, tolerance = 1e-6)
  }
})

test_that("systematic sampling refuses the jackknife too", {
  set.seed(3)
  pop <- data.frame(id = 1:240, y = round(stats::runif(240, 5, 200)))
  r <- ht_total(draw(pop, design_systematic(interval = 8), seed = 1,
                     weights = TRUE), "y", variance = "jackknife")
  expect_true(is.na(r$variance))
  expect_match(r$note, "one primary sampling unit")
})

test_that("variance = 'none' and 'analytic' behave as documented", {
  pop <- data.frame(id = 1:100, y = 1:100)
  s <- draw(pop, design_simple(n = 20), seed = 1, weights = TRUE)
  expect_true(is.na(ht_total(s, "y", variance = "none")$variance))
  expect_equal(ht_total(s, "y", variance = "analytic")$method, "analytic")
  expect_error(ht_total(s, "y", variance = "nope"))
})

# ---- simulated joint probabilities -----------------------------------------

test_that("simulated joints converge on the closed form where one exists", {
  d <- data.frame(id = 1:20, site = rep(c("a", "b"), each = 10))
  for (des in list(design_simple(n = 6),
                   design_stratified("site", n = 6),
                   design_systematic(interval = 4))) {
    exact <- joint_prob(d, des)
    sim <- joint_prob(d, des, simulate = TRUE, R = 4000, seed = 1)
    expect_equal(dim(sim), dim(exact))
    expect_lt(max(abs(sim - exact)), 0.035)
  }
})

test_that("simulation supplies joints for designs with no closed form", {
  set.seed(5)
  d <- data.frame(id = 1:30, w = stats::runif(30, 1, 10))
  des <- design_weighted("w", n = 8, method = "systematic")
  expect_error(joint_prob(d, des), "no closed-form joint")

  m <- joint_prob(d, des, simulate = TRUE, R = 3000, seed = 2)
  expect_equal(dim(m), c(30L, 30L))
  expect_equal(m, t(m))
  # The diagonal is the first-order probability, which this design does have.
  expect_lt(max(abs(diag(m) - inclusion_prob(d, des))), 0.04)
  # No pair can be likelier than either row on its own.
  expect_true(all(m <= outer(diag(m), diag(m), pmin) + 1e-9))
})

test_that("simulation respects `rows` and is reproducible", {
  d <- data.frame(id = 1:40)
  des <- design_simple(n = 10)
  a <- joint_prob(d, des, rows = c(2, 9, 30), simulate = TRUE, R = 500,
                  seed = 42)
  b <- joint_prob(d, des, rows = c(2, 9, 30), simulate = TRUE, R = 500,
                  seed = 42)
  expect_equal(dim(a), c(3L, 3L))
  expect_identical(a, b)
  expect_false(identical(
    a, joint_prob(d, des, rows = c(2, 9, 30), simulate = TRUE, R = 500,
                  seed = 43)))
})

test_that("simulation validates R and refuses a reserved column name", {
  d <- data.frame(id = 1:20)
  expect_error(joint_prob(d, design_simple(n = 5), simulate = TRUE, R = 0),
               "R")
  clash <- data.frame(id = 1:20, .drawn_row_id = 1:20)
  expect_error(joint_prob(clash, design_simple(n = 5), simulate = TRUE,
                          R = 10),
               "already has a column")
})

# ---- the sampling unit has to be the one the design actually uses ----------

test_that("cluster variance is honest when clusters differ in size", {
  # Whole clusters are taken, so the row count is random unless every cluster
  # is the same size. Sen-Yates-Grundy assumes a fixed size; applied row by row
  # here it came out at 0.39 of the truth, went negative on 57% of samples, and
  # produced zero-width intervals.
  set.seed(9)
  sizes <- c(3, 4, 5, 6, 3, 4, 5, 6, 7, 8, 4, 5)
  pop <- data.frame(cl = rep(paste0("c", seq_along(sizes)), times = sizes))
  pop$y <- round(stats::rnorm(nrow(pop), 100, 30) +
                   rep(stats::rnorm(length(sizes), 0, 40), times = sizes))
  des <- design_cluster("cl", n_clusters = 4)

  out <- vapply(1:1500, function(i) {
    e <- ht_total(draw(pop, des, seed = i, weights = TRUE), "y")
    c(e$total, e$variance)
  }, numeric(2))

  expect_true(all(out[2, ] >= 0))
  ratio <- mean(out[2, ]) / stats::var(out[1, ])
  expect_gt(ratio, 0.85)
  expect_lt(ratio, 1.2)

  # The cluster-level form is algebraically the delete-a-cluster jackknife.
  s <- draw(pop, des, seed = 1, weights = TRUE)
  expect_equal(ht_total(s, "y", variance = "analytic")$variance,
               ht_total(s, "y", variance = "jackknife")$variance,
               tolerance = 1e-8)

  # A design where every cluster total is identical has no variance at all,
  # and one where they differ must not report zero.
  flat <- data.frame(cl = rep(c("a", "b", "c", "d"), each = 5), y = 1)
  expect_equal(ht_total(draw(flat, design_cluster("cl", n_clusters = 2),
                             seed = 1, weights = TRUE), "y")$variance, 0)
})

test_that("the jackknife declines where deleting rows misrepresents the design", {
  set.seed(4)
  pop <- data.frame(w = stats::runif(60, 1, 20), y = round(stats::rnorm(60, 100, 25)))

  # Poisson: rows are independent and the exact form is right there.
  s <- draw(pop, design_weighted("w", n = 12, method = "poisson"), seed = 1,
            weights = TRUE)
  jk <- ht_total(s, "y", variance = "jackknife")
  expect_true(is.na(jk$variance))
  expect_match(jk$note, "variance = \"analytic\"")
  expect_true(is.finite(ht_total(s, "y", variance = "analytic")$variance))

  # Systematic PPS is still a fixed-size draw without replacement, so the
  # finite population correction belongs and is applied. Dropping it overstates
  # the variance badly once the sampling fraction is large -- 1.77 of the truth
  # at n/N = 0.375, against 1.11 with the correction.
  ps <- draw(pop, design_weighted("w", n = 12, method = "systematic"), seed = 1,
             weights = TRUE)
  r <- ht_total(ps, "y", variance = "jackknife")
  expect_equal(r$method, "jackknife")
  expect_match(r$note, "finite population correction")
})

test_that("the systematic-PPS jackknife tracks the real sampling variance", {
  set.seed(19)
  N <- 300
  pop <- data.frame(w = stats::rlnorm(N, 5, 0.9))
  # A response variable distinct from the size measure. Totalling the size
  # measure itself makes y/pi constant and every variance exactly zero, which
  # would let any estimator pass.
  pop$y <- round(pop$w * 0.3 + stats::rnorm(N, 0, 40))

  for (n in c(30, 90)) {
    des <- design_weighted("w", n = n, method = "systematic")
    o <- vapply(1:1200, function(i) {
      e <- ht_total(draw(pop, des, seed = i, weights = TRUE), "y")
      c(e$total, e$variance)
    }, numeric(2))
    ratio <- mean(o[2, ]) / stats::var(o[1, ])
    expect_gt(ratio, 0.75, label = paste("pps-systematic ratio at n =", n))
    expect_lt(ratio, 1.35, label = paste("pps-systematic ratio at n =", n))
  }
})

test_that("the jackknife treats certainty rows as fixed, not as sampling units", {
  # Certainty rows are in every possible sample. Deleting one and inflating the
  # rest invents variance the design does not have -- 4.3x too much.
  set.seed(12)
  pop <- data.frame(w = c(rep(100, 6), rep(1, 54)),
                    y = round(stats::rnorm(60, 100, 25), 3))
  des <- design_certainty("w", 50, design_simple(n = 10))
  s <- draw(pop, des, seed = 1, weights = TRUE)

  jk <- ht_total(s, "y", variance = "jackknife")
  expect_match(jk$note, "Jackknife over 10 primary")
  expect_false(grepl("over 16 primary", jk$note))
  expect_equal(jk$variance, ht_total(s, "y", variance = "analytic")$variance,
               tolerance = 1e-8)
})

test_that("a negative analytic variance is reported, not passed on", {
  # Sen-Yates-Grundy can go negative on an unlucky sample. sqrt() was already
  # guarded, but the reason was dropped, `auto` did not fall back, and deff()
  # divided by it and returned a negative design effect.
  fake <- list(variance = -1234, method = "analytic", note = NULL)
  expect_true(is.na(drawn:::deff_value(-1234, 1:10, rep(0.5, 10), 100, 10,
                                       "total")))

  pop <- data.frame(id = 1:60, y = 1:60)
  s <- draw(pop, design_simple(n = 20), seed = 1, weights = TRUE)
  bad <- drawn:::ht_variance_dispatch(
    structure(list(), class = c("drawn_design_nope", "drawn_design")),
    s, pop, attr(s, "drawn_rows"), s$y, s$.prob, "analytic")
  expect_true(is.na(bad$variance))
  expect_true(nzchar(bad$note))
})

test_that("a design neither estimator can handle says so once, truthfully", {
  set.seed(3)
  pop <- data.frame(y = stats::rnorm(60, 10, 3), cl = rep(1:6, each = 10))
  r <- ht_total(draw(pop, design_cluster("cl", n_clusters = 1), seed = 1,
                     weights = TRUE), "y")
  # Nothing computed a variance, so nothing may claim to have been used.
  expect_true(is.na(r$variance))
  expect_equal(r$method, "none")
  expect_match(r$note, "at least two clusters")
  expect_match(r$note, "at least two primary")
})

test_that("strata and clusters come from the frame, not the sample's columns", {
  d <- var_pop()
  designs <- list(
    site = design_stratified("site", n = 40, min_per_stratum = 2),
    cl = design_cluster("cl", n_clusters = 6),
    cl = design_multistage("cl", n_clusters = 6, n = 24)
  )
  for (i in seq_along(designs)) {
    key <- names(designs)[i]
    s <- draw(d, designs[[i]], seed = 2, weights = TRUE)
    want <- ht_total(s, "y")
    want_jk <- ht_total(s, "y", variance = "jackknife")
    renamed <- s
    names(renamed)[names(renamed) == key] <- "renamed"
    recoded <- s
    recoded[[key]] <- rev(recoded[[key]])
    for (e in list(renamed, recoded)) {
      got <- ht_total(e, "y")
      expect_equal(got$se, want$se, label = key)
      expect_equal(got$df, want$df, label = key)
      expect_equal(ht_total(e, "y", variance = "jackknife")$se, want_jk$se,
                   label = key)
    }
  }
})

test_that("dropping the cluster column leaves a cluster sample's variance alone", {
  # Without the column the jackknife used to delete rows rather than clusters
  # and the degrees of freedom fell to 0, which removed the interval.
  d <- var_pop()
  s <- draw(d, design_cluster("cl", n_clusters = 6), seed = 1, weights = TRUE)
  gone <- s
  gone$cl <- NULL
  got <- ht_total(gone, "y")
  expect_equal(got$df, 5)
  expect_equal(got$ci, ht_total(s, "y")$ci)
  jk <- ht_total(gone, "y", variance = "jackknife")
  expect_match(jk$note, "Jackknife over 6 primary")
  expect_equal(jk$se, ht_total(s, "y", variance = "jackknife")$se)
})

test_that("a variance that overflows is reported as missing, with no success note", {
  d <- var_pop()
  s <- draw(d, design_simple(n = 30), seed = 1, weights = TRUE)
  huge <- rep(c(1e300, -1e300), 15)
  r <- ht_total(s, huge)
  expect_true(is.na(r$se))
  expect_equal(r$method, "none")
  expect_match(r$note, "overflowed")
  expect_false(grepl("Jackknife over", r$note))
})

test_that("an overflowing variance names its cause once, in every branch", {
  d <- var_pop()
  s <- draw(d, design_reservoir(n = 30), seed = 1, weights = TRUE)
  big <- s$y * 1e200
  for (v in c("auto", "analytic", "jackknife")) {
    r <- ht_total(s, big, variance = v)
    expect_true(is.finite(r$total), label = v)
    expect_true(is.na(r$variance), label = v)
    expect_equal(r$method, "none", label = v)
    expect_equal(lengths(regmatches(r$note, gregexpr("Rescale", r$note))), 1L,
                 label = v)
  }
})

# ---- what the dispatcher catches ---------------------------------------------

test_that("an error inside a variance estimator reaches the caller", {
  # Every error used to be read as the estimator declining, so running out of
  # memory, or a bug, came back as a jackknife figure with the error's text in
  # the note.
  pop <- data.frame(id = 1:60, y = 1:60)
  s <- draw(pop, design_simple(n = 20), seed = 1, weights = TRUE)
  # `by` computes every domain's simple random sampling variance at once, with
  # srs_variance_by().
  local_mocked_bindings(srs_variance = function(...) stop("internal failure"),
                        srs_variance_by = function(...) {
                          stop("internal failure")
                        })
  expect_error(ht_total(s, "y"), "internal failure")
  expect_error(ht_total(s, "y", variance = "analytic"), "internal failure")
  expect_error(ht_mean(s, "y"), "internal failure")
  expect_error(ht_total(s, "y", by = "id"), "internal failure")
})

test_that("an error inside the jackknife reaches the caller", {
  pop <- data.frame(id = 1:60, y = 1:60)
  s <- draw(pop, design_simple(n = 20), seed = 1, weights = TRUE)
  local_mocked_bindings(
    jackknife_groups = function(...) stop("internal failure")
  )
  expect_error(ht_total(s, "y", variance = "jackknife"), "internal failure")
})

test_that("a declining estimator still falls back to the jackknife", {
  pop <- data.frame(id = 1:60, y = 1:60)
  s <- draw(pop, design_simple(n = 20), seed = 1, weights = TRUE)
  want <- ht_total(s, "y", variance = "jackknife")
  local_mocked_bindings(srs_variance = function(...) {
    drawn:::stop_declined("No joint probabilities here.\nA second line.")
  })
  r <- ht_total(s, "y")
  expect_equal(r$method, "jackknife")
  expect_equal(r$variance, want$variance)
  expect_match(r$note,
               "jackknife was used instead. No joint probabilities here.",
               fixed = TRUE)
  expect_false(grepl("second line", r$note))

  a <- ht_total(s, "y", variance = "analytic")
  expect_true(is.na(a$variance))
  expect_equal(a$method, "none")
  expect_match(a$note, "No joint probabilities here.", fixed = TRUE)
})

test_that("a refusal from a design is a classed error outside the estimators", {
  d <- var_pop()
  err <- tryCatch(joint_prob(d, design_weighted("w", n = 10)),
                  error = function(e) e)
  expect_s3_class(err, "drawn_declined")
  expect_null(conditionCall(err))
  err <- tryCatch(inclusion_prob(d, design_systematic(8, start = 2)),
                  error = function(e) e)
  expect_s3_class(err, "drawn_declined")
})

test_that("a simple random sample of 5,000 rows gets its analytic variance", {
  # Sen-Yates-Grundy held about eleven n x n matrices of doubles at once;
  # 8,000 rows needed 6 GB. The closed form needs none.
  pop <- data.frame(id = 1:10000, y = (1:10000) %% 97)
  s <- draw(pop, design_simple(n = 5000), seed = 1, weights = TRUE)
  r <- ht_total(s, "y")
  expect_equal(r$method, "analytic")
  expect_null(r$note)
  f <- 1 - 5000 / 10000
  expect_equal(r$variance, 10000^2 * f * stats::var(s$y) / 5000,
               tolerance = 1e-8)
  expect_equal(ht_total(s, "y", variance = "analytic")$variance, r$variance)
})

test_that("Sen-Yates-Grundy declines a sample too large to hold its matrices", {
  # A design with no closed form of its own reaches the pairwise sum. It
  # declines from the projected size, before allocating anything.
  pop <- data.frame(id = 1:10000, y = (1:10000) %% 97)
  rows <- seq(1L, 10000L, by = 2L)
  fake <- structure(list(), class = c("drawn_design_nope", "drawn_design"))
  local_mocked_bindings(joint_inclusion = function(...) {
    stop("the joint matrix was built")
  })
  err <- tryCatch(drawn:::ht_variance(fake, pop, rows, pop$y[rows],
                                      rep(0.5, 5000)),
                  error = function(e) e)
  expect_s3_class(err, "drawn_declined")
  expect_match(conditionMessage(err), "1.6 GB of memory for its 5000 x 5000")

  r <- drawn:::ht_variance_dispatch(fake, pop[rows, ], pop, rows,
                                    pop$y[rows], rep(0.5, 5000), "auto")
  expect_equal(r$method, "jackknife")
  expect_match(r$note, "GB of memory")
  expect_equal(r$variance, 10000^2 * 0.5 * stats::var(pop$y[rows]) / 5000,
               tolerance = 1e-8)
  a <- drawn:::ht_variance_dispatch(fake, pop[rows, ], pop, rows,
                                    pop$y[rows], rep(0.5, 5000), "analytic")
  expect_true(is.na(a$variance))
  expect_match(a$note, "variance = \"auto\"", fixed = TRUE)
})

test_that("a jackknife that gives no figure reports method 'none'", {
  # ?ht_total says `method` is "none" when nothing produced a figure, but a
  # refused jackknife said "jackknife", in the by table too.
  pop <- data.frame(w = 1 + (1:60 * 7) %% 19, y = 100 + (1:60 * 37) %% 41,
                    g = rep(c("a", "b"), 30))
  for (des in list(design_weighted("w", n = 12, method = "poisson"),
                   design_systematic(interval = 6))) {
    s <- draw(pop, des, seed = 1, weights = TRUE)
    r <- ht_total(s, "y", variance = "jackknife")
    expect_true(is.na(r$se))
    expect_equal(r$method, "none")
    expect_true(nzchar(r$note))
    b <- ht_total(s, "y", by = "g", variance = "jackknife")
    expect_true(all(is.na(b$se)))
    expect_true(all(b$method == "none"))
  }
})

test_that("a Hajek mean of a one-row sample has no standard error", {
  # The residual of one row is 0, so the Poisson estimator gave 0 and a
  # zero-width interval.
  pop <- data.frame(w = rep(1, 40), y = 50 + (1:40 * 37) %% 23)
  des <- design_weighted("w", n = 1, method = "poisson")
  sizes <- vapply(1:40, function(i) nrow(draw(pop, des, seed = i)), integer(1))
  s <- draw(pop, des, seed = which(sizes == 1L)[1], weights = TRUE)
  m <- ht_mean(s, "y")
  expect_true(is.na(m$se))
  expect_equal(m$method, "none")
  expect_match(m$note, "one row")
  expect_true(all(is.na(m$ci)))
  expect_equal(ht_mean(s, "y", variance = "none")$note,
               "Variance not requested.")
  # A total, and the mean over the known N, are still estimated.
  expect_true(is.finite(ht_total(s, "y")$se))
  expect_true(is.finite(ht_mean(s, "y", estimator = "ht")$se))

  # A census of one row is exact.
  one <- draw(data.frame(y = 7), design_simple(n = 1), weights = TRUE)
  expect_equal(ht_mean(one, "y")$se, 0)
})

test_that("an interval missing for want of degrees of freedom says why", {
  # A one-row Poisson sample has a standard error and 0 degrees of freedom.
  # The note blamed one unit per stratum, and print left it out.
  pop <- data.frame(w = rep(1, 40), y = 50 + (1:40 * 37) %% 23)
  des <- design_weighted("w", n = 1, method = "poisson")
  sizes <- vapply(1:40, function(i) nrow(draw(pop, des, seed = i)), integer(1))
  s <- draw(pop, des, seed = which(sizes == 1L)[1], weights = TRUE)
  r <- ht_total(s, "y")
  expect_equal(r$df, 0)
  expect_true(is.finite(r$se))
  expect_true(all(is.na(r$ci)))
  expect_match(r$note, "leaves 0 degrees of freedom")
  expect_false(grepl("per stratum", r$note))
  expect_output(print(r), "leaves 0 degrees")
  expect_true(all(is.finite(ht_total(s, "y", df = Inf)$ci)))
})

test_that("a df below 1 is refused rather than dropping the interval", {
  pop <- data.frame(id = 1:100, y = 1:100)
  s <- draw(pop, design_simple(n = 20), seed = 1, weights = TRUE)
  expect_error(ht_total(s, "y", df = 0.5), "not 0.5")
  expect_error(ht_mean(s, "y", df = 0.5), "1 or more")
  r <- ht_total(s, "y", df = 1)
  expect_equal(diff(r$ci) / 2, stats::qt(0.975, 1) * r$se)
})

test_that("an exact estimate gets a zero-width interval at any df", {
  expect_equal(drawn:::interval(5, 0, 0.95, 0), c(5, 5))
  expect_equal(drawn:::interval(0.3, 0, 0.95, 0, logit = TRUE), c(0.3, 0.3))
  expect_equal(drawn:::interval(5, 1, 0.95, 0), c(NA_real_, NA_real_))
})

# ---- degrees of freedom and units taken with certainty ----------------------

dominant8 <- function() {
  # Eight units large enough that systematic PPS of 10 takes them with
  # certainty, and 40 small ones from which the other two are drawn.
  pop <- data.frame(x = c(rep(100, 8), 1 + 3 * (sin((1:40) * 2.3) + 1)))
  pop$y <- 2 * pop$x + 3 * cos(1:48)
  pop$half <- rep(c("a", "b"), 24)
  pop
}

test_that("units taken with certainty carry no degrees of freedom", {
  pop <- dominant8()
  pps <- design_weighted("x", n = 10, method = "systematic")
  cert <- design_certainty("x", 100, design_weighted("x", n = 2,
                                                     method = "systematic"))
  expect_equal(sort(inclusion_prob(pop, pps)), sort(inclusion_prob(pop, cert)))

  # The same selection built two ways had 9 and 1 degrees of freedom.
  s <- draw(pop, pps, seed = 1, weights = TRUE)
  expect_equal(sum(s$.prob == 1), 8L)
  r <- ht_total(s, "y")
  expect_equal(r$df, 1)
  expect_equal(ht_total(draw(pop, cert, seed = 1, weights = TRUE), "y")$df, 1)
  expect_equal(diff(r$ci) / 2, stats::qt(0.975, 1) * r$se)
  expect_equal(ht_mean(s, "y")$df, 1)
  expect_equal(attr(ht_total(s, "y", by = "half"), "df"), 1)
})

test_that("the t interval covers once certain units stop counting", {
  skip_on_cran()
  # On 9 degrees of freedom this interval covered 74% of the time.
  pop <- dominant8()
  truth <- sum(pop$y)
  des <- design_weighted("x", n = 10, method = "systematic")
  hit <- vapply(1:400, function(i) {
    r <- ht_total(draw(pop, des, seed = i, weights = TRUE), "y")
    r$ci[1] <= truth && truth <= r$ci[2]
  }, logical(1))
  expect_gt(mean(hit), 0.9)
})

test_that("a stratum taken whole carries no degrees of freedom", {
  pop <- data.frame(g = c(rep("all", 5), rep("big", 60)),
                    y = c(1000 * (1:5), 1:60))
  s <- draw(pop, design_stratified("g", n = 25, allocation = "equal"),
            seed = 2, weights = TRUE)
  expect_true(all(s$.prob[s$g == "all"] == 1))
  # 20 rows in the stratum that was sampled, less that one stratum
  expect_equal(ht_total(s, "y")$df, 19)
})

test_that("a sample with nothing random has infinite degrees of freedom", {
  pop <- data.frame(g = rep(c("a", "b"), 10), w = 1:20, y = (1:20)^1.5)
  for (des in list(design_simple(n = 20), design_stratified("g", n = 20),
                   design_weighted("w", n = 20, method = "systematic"))) {
    r <- ht_total(draw(pop, des, seed = 1, weights = TRUE), "y")
    expect_equal(r$df, Inf)
    expect_equal(r$variance, 0)
    expect_equal(r$ci, rep(sum(pop$y), 2))
  }
})
