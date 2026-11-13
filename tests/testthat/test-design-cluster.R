test_that("selects whole clusters", {
  df <- make_df()
  res <- draw(df, design_cluster("site", n_clusters = 2), seed = 123)
  expect_equal(length(unique(res$site)), 2)
  expect_equal(nrow(res), 50)
})

test_that("balanced = TRUE equalises cluster sizes", {
  df <- make_df()
  uneven <- df[c(1:25, 26:45, 51:65, 76:100), ]
  res <- draw(uneven, design_cluster("site", n_clusters = 4, balanced = TRUE),
              seed = 123)
  expect_equal(length(unique(as.vector(table(res$site)))), 1)
})

test_that("asking for more clusters than exist is an error", {
  df <- make_df()
  expect_error(draw(df, design_cluster("site", n_clusters = 10), seed = 1),
               "exceeds the 4 cluster")
})

test_that("NA cluster labels are reported rather than treated as a cluster", {
  df <- data.frame(x = 1:30, cl = c(rep("a", 10), rep("b", 10), rep(NA, 10)))
  expect_error(draw(df, design_cluster("cl", n_clusters = 2), seed = 1),
               "missing cluster label")

  res <- draw(df, design_cluster("cl", n_clusters = 2, na_rm = TRUE), seed = 1)
  expect_false(anyNA(res$cl))
  expect_equal(nrow(res), 20)
})

test_that("a single numeric cluster id is not read as a range", {
  # sample(5, 1) samples from 1:5, so the only real cluster was never selected.
  d <- data.frame(x = 1:10, cl = rep(5, 10))
  expect_equal(nrow(draw(d, design_cluster("cl", n_clusters = 1), seed = 1)),
               10)
})

test_that("clusters must name a single column", {
  expect_error(design_cluster(c("a", "b"), n_clusters = 2), "a single column")
})

test_that("the clusters chosen for a seed are sample(unique(labels))'s", {
  d <- data.frame(id = 1:60, cl = rep(c(5, 2, 9, 1, 7, 3), each = 10))
  for (k in c(1, 3)) {
    got <- unique(draw(d, design_cluster("cl", n_clusters = k), seed = 7)$cl)
    want <- drawn:::with_seed(7, sample(unique(d$cl), k))
    expect_setequal(got, want)
  }
})

test_that("two cluster ids that print alike stay two clusters", {
  # 0.3 and 0.1 + 0.2 both print as 0.3. Compared as text they were one
  # cluster: a two-cluster sample got no variance, and two rows of different
  # clusters got the same-cluster joint probability.
  d <- data.frame(id = 1:20, cl = rep(c(0.3, 0.1 + 0.2, 0.5, 0.7), each = 5),
                  y = c(1:5, 101:105, 21:25, 61:65))
  des <- design_cluster("cl", n_clusters = 3)
  expect_equal(inclusion_prob(d, des), rep(0.75, 20))
  m <- joint_prob(d, des, rows = c(1, 2, 6))
  expect_equal(m[1, 2], 0.75)
  expect_equal(m[1, 3], 3 * 2 / (4 * 3))

  s <- draw(d, des, seed = 2, weights = TRUE)
  expect_equal(sum(s$cl == 0.3), 5)
  expect_equal(sum(s$cl == 0.1 + 0.2), 5)
  u <- tapply(s$y / s$.prob, match(s$cl, unique(s$cl)), sum)
  r <- ht_total(s, "y")
  expect_equal(r$method, "analytic")
  expect_equal(r$variance, (1 - 3 / 4) * 3 * stats::var(u))
  expect_equal(r$df, 2)
  expect_equal(ht_total(s, "y", variance = "jackknife")$variance, r$variance)
  expect_equal(nrow(sample_summary(s)$by_group), 3)

  skip_if_not_installed("survey")
  des_sv <- as_svydesign(s)
  expect_equal(survey::degf(des_sv), 2)
  expect_equal(as.numeric(survey::SE(survey::svytotal(~y, des_sv))), r$se)
})

test_that("a factor's unused levels are not clusters", {
  d <- data.frame(id = 1:40, y = 1:40,
                  cl = factor(rep(c("a", "b", "c", "d"), each = 10),
                              levels = c("a", "b", "c", "d", "zz")))
  des <- design_cluster("cl", n_clusters = 2)
  expect_equal(inclusion_prob(d, des), rep(0.5, 40))
  s <- draw(d, des, seed = 1, weights = TRUE)
  expect_equal(nrow(s), 20)
  r <- ht_total(s, "y")
  expect_equal(r$method, "analytic")
  expect_equal(r$df, 1)
  expect_equal(nrow(sample_summary(s)$by_group), 2)
})

test_that("a missing cluster label is called that everywhere", {
  d <- data.frame(id = 1:30, cl = c(NA, rep(c("a", "b", "c"), each = 10))[1:30])
  des <- design_cluster("cl", n_clusters = 2)
  expect_error(draw(d, des, seed = 1), "1 row\\(s\\) have a missing cluster label")
  expect_error(inclusion_prob(d, des),
               "1 row\\(s\\) have a missing cluster label")
  expect_error(inclusion_prob(d, design_multistage("cl", n_clusters = 2, n = 4)),
               "1 row\\(s\\) have a missing cluster label")
})

test_that("selecting many clusters takes one pass over the frame", {
  skip_on_cran()
  # Each chosen cluster scanned the whole frame: 16 seconds for 1,000
  # clusters of a million rows.
  N <- 1e6
  d <- data.frame(cl = sprintf("c%07d", (seq_len(N) - 1) %/% 50 + 1))
  took <- system.time({
    s <- draw(d, design_cluster("cl", n_clusters = 1000), seed = 1)
    m <- draw(d, design_multistage("cl", n_clusters = 1000, n = 5000), seed = 1)
  })[["elapsed"]]
  expect_equal(nrow(s), 50000)
  expect_equal(length(unique(s$cl)), 1000)
  expect_equal(nrow(m), 5000)
  expect_lt(took, 5)
})

# ---- balanced clusters ------------------------------------------------------

test_that("balanced clusters have exact inclusion probabilities", {
  # The take is the smallest selected cluster's size. Given cluster i, the
  # other cluster is one of three, so clusters of 2/4/6/8 average takes of
  # 2, 10/3, 4 and 4.
  d <- data.frame(cl = rep(c("a", "b", "c", "d"), c(2, 4, 6, 8)))
  des <- design_cluster("cl", n_clusters = 2, balanced = TRUE)
  p <- inclusion_prob(d, des)
  expect_equal(as.numeric(tapply(p, d$cl, unique)),
               c(1 / 2, 5 / 12, 1 / 3, 1 / 4))
  expect_equal(sampling_weight(d, des), 1 / p)
  sim <- inclusion_prob(d, des, simulate = TRUE, R = 2000, seed = 1)
  expect_lt(max(abs(tapply(sim, d$cl, mean) - tapply(p, d$cl, mean))), 0.03)

  # Every selection of 3 clusters out of 7, of sizes with ties, enumerated.
  sizes <- c(3, 7, 1, 7, 4, 9, 4)
  e <- data.frame(cl = rep(seq_along(sizes), sizes))
  got <- inclusion_prob(e, design_cluster("cl", n_clusters = 3,
                                          balanced = TRUE))
  sets <- utils::combn(7, 3)
  want <- numeric(7)
  for (j in seq_len(ncol(sets))) {
    s <- sets[, j]
    want[s] <- want[s] + min(sizes[s]) / sizes[s] / ncol(sets)
  }
  expect_equal(as.numeric(tapply(got, e$cl, unique)), want)

  # The probability plot needs them too.
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit({grDevices::dev.off(); unlink(f)}, add = TRUE)
  expect_no_error(plot(des, d, type = "probability"))

  # One cluster selected is taken whole; all of them, at the smallest size.
  expect_equal(inclusion_prob(e, design_cluster("cl", n_clusters = 1,
                                                balanced = TRUE)),
               rep(1 / 7, nrow(e)))
  expect_equal(inclusion_prob(e, design_cluster("cl", n_clusters = 7,
                                                balanced = TRUE)),
               1 / rep(sizes, sizes))
})

test_that("a missing cluster label gets probability 0 under balance too", {
  d <- data.frame(cl = c(NA, rep(c("a", "b", "c"), c(2, 4, 6))))
  des <- design_cluster("cl", n_clusters = 2, balanced = TRUE, na_rm = TRUE)
  p <- inclusion_prob(d, des)
  expect_equal(p[1], 0)
  # Average takes of 2, 3 and 3, with two of three clusters selected.
  expect_equal(p[-1], rep(c(2 / 3, 1 / 2, 1 / 3), c(2, 4, 6)))
})

test_that("a balanced sample is weighted conditional on its clusters", {
  sizes <- c(3, 5, 7, 4, 6, 8, 5, 9)
  d <- data.frame(id = seq_len(sum(sizes)),
                  cl = rep(paste0("k", 1:8), sizes))
  d$y <- (d$id * 37) %% 23 + rep(c(10, 40, 25, 60, 15, 30, 50, 20), sizes)
  des <- design_cluster("cl", n_clusters = 4, balanced = TRUE)
  s <- draw(d, des, seed = 3, weights = TRUE)

  sel <- attr(s, "drawn_selection")
  expect_type(sel, "integer")
  expect_setequal(sel, match(unique(s$cl), unique(d$cl)))
  m <- min(sizes[sel])
  expect_equal(as.vector(table(s$cl)), rep(m, 4))
  expect_equal(s$.prob, 4 / 8 * m / sizes[match(s$cl, unique(d$cl))])
  pi_all <- attr(s, "drawn_pi")
  expect_true(all(is.na(pi_all[!d$cl %in% s$cl])))
  expect_equal(pi_all[s$.row], s$.prob)

  # The two-stage formula, at the take of m rows from each cluster.
  r <- ht_total(s, "y")
  N_i <- sizes[sel]
  g <- factor(s$cl, levels = unique(d$cl)[sel])
  tot <- N_i * tapply(s$y, g, mean)
  s2 <- tapply(s$y, g, stats::var)
  v <- 8^2 * (1 - 4 / 8) * stats::var(tot) / 4 +
    8 / 4 * sum(N_i^2 * (1 - m / N_i) * s2 / m)
  expect_equal(r$total, 8 / 4 * sum(tot))
  expect_equal(r$method, "analytic")
  expect_equal(r$variance, v)
  expect_equal(r$df, 3)
  expect_equal(ht_total(s, "y", variance = "jackknife")$method, "jackknife")

  sm <- sample_summary(s)
  expect_equal(sm$conditional, 4L)
  expect_equal(sm$N_reachable, nrow(d))
  expect_equal(sm$unreachable, 0L)
  expect_output(print(sm), "conditional on the 4 clusters selected")
  expect_true(is.na(sample_summary(draw(d, design_cluster("cl", 4), seed = 3,
                                        weights = TRUE))$conditional))

  # Without weights the selection is not attached.
  expect_null(attr(draw(d, des, seed = 3), "drawn_selection"))
})

test_that("conditional weights are unbiased over every selection", {
  # Given the clusters, each is estimated without bias; the clusters are a
  # simple random sample. Averaging the expected estimate over all 35
  # selections gives the total. `.prob` is a / A times the chance of a row
  # given the selection, and that chance averages to the exact probability.
  sizes <- c(3, 7, 1, 7, 4, 9, 4)
  d <- data.frame(cl = rep(seq_along(sizes), sizes))
  d$y <- seq_len(nrow(d))^2 %% 17 + d$cl * 5
  des <- design_cluster("cl", n_clusters = 3, balanced = TRUE)
  sets <- utils::combn(7, 3)
  est <- numeric(ncol(sets))
  avg <- numeric(nrow(d))
  for (j in seq_len(ncol(sets))) {
    p <- drawn:::conditional_inclusion(des, d, sets[, j])
    take <- min(sizes[sets[, j]])
    hit <- !is.na(p)
    # A row of cluster i is drawn with chance take / N_i given the clusters.
    est[j] <- sum((d$y * take / sizes[d$cl] / p)[hit])
    avg[hit] <- avg[hit] + p[hit] * 7 / 3 / ncol(sets)
  }
  expect_equal(mean(est), sum(d$y))
  expect_equal(avg, inclusion_prob(d, des))
})

test_that("a balanced design under certainty declines the two-stage form", {
  # Inside a certainty design the cluster rows are weighted with the
  # unconditional probabilities, which the two-stage estimator cannot use.
  sizes <- rep(c(3, 5, 7), 4)
  d <- data.frame(cl = rep(paste0("c", 1:12), sizes))
  d$v <- seq_len(nrow(d))
  d$y <- (d$v * 13) %% 11 + 2
  des <- design_certainty("v", 58, design_cluster("cl", n_clusters = 4,
                                                  balanced = TRUE))
  s <- draw(d, des, seed = 1, weights = TRUE)
  expect_null(attr(s, "drawn_selection"))
  r <- ht_total(s, "y")
  expect_equal(r$method, "jackknife")
  expect_match(r$note, "no closed-form joint inclusion probability")
  expect_equal(ht_total(s, "y", variance = "analytic")$method, "none")
})
