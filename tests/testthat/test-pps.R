test_that("systematic PPS gives inclusion probabilities proportional to size", {
  df <- data.frame(id = 1:8, w = c(1, 1, 2, 2, 3, 3, 4, 4))
  p <- inclusion_prob(df, design_weighted("w", n = 3, method = "systematic"))
  # Constant ratio to weight is the definition of proportional-to-size.
  expect_equal(length(unique(round(p / df$w, 10))), 1L)
  expect_equal(sum(p), 3)
})

test_that("systematic PPS draws exactly n rows", {
  df <- data.frame(id = 1:20, w = 1:20)
  for (i in 1:20) {
    expect_equal(
      nrow(draw(df, design_weighted("w", n = 6, method = "systematic"), seed = i)),
      6
    )
  }
})

test_that("units too heavy for a valid probability are taken with certainty", {
  df <- data.frame(id = 1:5, w = c(1, 1, 1, 1, 16))
  p <- inclusion_prob(df, design_weighted("w", n = 2, method = "systematic"))
  # n * p would be 1.6 for the heavy unit, which is not a probability.
  expect_equal(p[5], 1)
  expect_equal(p[1:4], rep(0.25, 4))
  expect_equal(sum(p), 2)

  # And it really is in every sample.
  ids <- vapply(1:50, function(i) {
    5 %in% draw(df, design_weighted("w", n = 2, method = "systematic"),
                seed = i)$id
  }, logical(1))
  expect_true(all(ids))
})

test_that("certainty treatment iterates when removing units creates new ones", {
  df <- data.frame(id = 1:5, w = c(1, 1, 1, 30, 60))
  p <- inclusion_prob(df, design_weighted("w", n = 3, method = "systematic"))
  expect_true(all(p <= 1))
  expect_equal(sum(p), 3)
  expect_equal(p[4:5], c(1, 1))
})

test_that("a unit whose share comes to 1 gets exactly 1", {
  # 5 * (2.3 / 11.5) is 1 in exact arithmetic and came out 1 - 1e-16 in
  # floating point, while the same sizes times 10 gave exactly 1.
  pop <- data.frame(id = 1:8, w = c(1.1, 0.2, 1.8, 1.3, 1.5, 1.6, 2.3, 1.7))
  d <- design_weighted("w", n = 5, method = "systematic")
  expect_identical(inclusion_prob(pop, d)[7], 1)
  expect_identical(inclusion_prob(transform(pop, w = w * 10), d)[7], 1)
  expect_equal(sum(inclusion_prob(pop, d)), 5)

  # A census then reads as one: every .prob is 1 and the total is exact,
  # where it had a small standard error under a Deville note.
  pop3 <- data.frame(id = 1:41, w = 0.1, y = 100 + (1:41))
  s <- draw(pop3, design_weighted("w", n = 41, method = "systematic"),
            seed = 1, weights = TRUE)
  expect_true(all(s$.prob == 1))
  r <- ht_total(s, "y")
  expect_equal(r$variance, 0)
  expect_match(r$note, "exact rather than estimated")
})

test_that("realised inclusion matches the claim for systematic PPS", {
  df <- data.frame(id = 1:8, w = c(1, 1, 2, 2, 3, 3, 4, 4))
  d <- design_weighted("w", n = 3, method = "systematic")
  claim <- inclusion_prob(df, d)

  R <- 3000
  hits <- integer(8)
  for (i in seq_len(R)) {
    hits[draw(df, d, seed = i)$id] <- hits[draw(df, d, seed = i)$id] + 1L
  }
  expect_equal(hits / R, claim, tolerance = 0.08)
})

test_that("poisson sampling is proportional but has a random size", {
  df <- data.frame(id = 1:20, w = 1:20)
  d <- design_weighted("w", n = 6, method = "poisson")
  p <- inclusion_prob(df, d)
  expect_equal(sum(p), 6)

  sizes <- vapply(1:300, function(i) nrow(draw(df, d, seed = i)), numeric(1))
  expect_gt(stats::var(sizes), 0)          # random, unlike systematic
  expect_equal(mean(sizes), 6, tolerance = 0.5)
})

test_that("successive sampling really does differ from proportional", {
  # The point of offering the other two methods. With a moderate weight spread
  # the gap is large; with one dominant unit the two nearly coincide, which is
  # why a single example is a poor way to judge it.
  df <- data.frame(id = 1:10, w = 1:10)
  target <- inclusion_prob(df, design_weighted("w", n = 5, method = "systematic"))

  R <- 2000
  hits <- integer(10)
  for (i in seq_len(R)) {
    ids <- draw(df, design_weighted("w", n = 5), seed = i)$id
    hits[unique(ids)] <- hits[unique(ids)] + 1L
  }
  obs <- hits / R
  expect_gt(max(abs(obs - target) / target), 0.2)
})

test_that("replace = TRUE cannot be combined with the PPS methods", {
  expect_error(design_weighted("w", n = 5, replace = TRUE, method = "systematic"),
               "cannot be combined")
  expect_error(design_weighted("w", n = 5, replace = TRUE, method = "poisson"),
               "cannot be combined")
  expect_s3_class(design_weighted("w", n = 5, replace = TRUE), "drawn_design")
})

test_that("PPS methods refuse to draw more rows than exist", {
  df <- data.frame(id = 1:5, w = 1:5)
  expect_error(draw(df, design_weighted("w", n = 10, method = "systematic"),
                    seed = 1),
               "cannot exceed the number of rows")
})

test_that("an unknown method is rejected at construction", {
  expect_error(design_weighted("w", n = 5, method = "pps"))
})

test_that("an empty Poisson sample estimates a total of 0 with no variance", {
  # Poisson sampling can draw nothing. Its total is then 0, which keeps the
  # estimator unbiased over repeated draws, so a simulation loop gets a number
  # instead of an error. One empty draw has no variance to report.
  pop <- data.frame(id = 1:50, w = 1 + (1:50 %% 4), y = 1:50)
  des <- design_weighted("w", n = 1, method = "poisson")
  for (i in 1:100) {
    s <- draw(pop, des, seed = i, weights = TRUE)
    if (!nrow(s)) break
  }
  expect_equal(nrow(s), 0L)

  r <- ht_total(s, "y")
  expect_equal(r$total, 0)
  expect_true(is.na(r$variance))
  expect_true(all(is.na(r$ci)))
  expect_equal(r$df, 0)
  expect_equal(r$method, "none")
  expect_match(r$note, "The sample is empty")
  expect_output(print(r), "The sample is empty")
  expect_equal(nrow(ht_total(s, "y", by = "id")), 0L)

  # The Hajek mean is 0 / 0; the Horvitz-Thompson mean divides by N
  expect_error(ht_mean(s, "y"), "Hajek mean is 0 / 0")
  expect_equal(ht_mean(s, "y", estimator = "ht")$mean, 0)

  skip_if_not_installed("survey")
  expect_error(as_svydesign(s), "`sample` has no rows", fixed = TRUE)
})

test_that("an empty sample from a design that reaches no row is refused", {
  pop <- data.frame(id = 1:20, w = 1:20, y = 1:20)
  s <- draw(pop, design_simple(n = 0), seed = 1, weights = TRUE)
  expect_error(ht_total(s, "y"), "probability of 0")
})

test_that("weights whose total overflows still give probabilities summing to n", {
  # Positive weights near 1e306 summed to Inf, which made every probability 0
  # and every sample empty.
  pop <- data.frame(id = 1:40, w = c(1:39, 100))
  pop$big <- pop$w * 1e306
  pop$small <- pop$w * 1e-300
  expect_false(is.finite(sum(pop$big)))
  for (m in c("systematic", "poisson")) {
    want <- inclusion_prob(pop, design_weighted("w", n = 10, method = m))
    for (col in c("big", "small")) {
      des <- design_weighted(col, n = 10, method = m)
      p <- inclusion_prob(pop, des)
      expect_equal(p, want, info = paste(m, col))
      expect_equal(sum(p), 10)
    }
  }
  expect_equal(nrow(draw(pop, design_weighted("big", n = 10,
                                              method = "systematic"),
                         seed = 1)), 10)
  expect_equal(nrow(draw(pop, design_weighted("big", n = 10), seed = 1)), 10)
  s <- draw(pop, design_spread("id", n = 10, size = "big"), seed = 1,
            weights = TRUE)
  expect_equal(nrow(s), 10)
  expect_equal(sum(inclusion_prob(pop, design_spread("id", n = 10,
                                                     size = "big"))), 10)
  cert <- design_certainty("id", 40, design_weighted("big", n = 5,
                                                     method = "systematic"))
  expect_equal(nrow(draw(pop, cert, seed = 1)), 6)
})

test_that("weights a double cannot scale together are refused, not zeroed", {
  # Divided by the largest, weights 1e-320 apart in scale from it underflow
  # to 0 and leave the probabilities short of n.
  w <- c(1e308, 1e308, 1e-320, 1e-320)
  expect_error(drawn:::pps_pi(w, 3), "summing to 2, where `n` is 3",
               fixed = TRUE)
  expect_equal(drawn:::pps_pi(c(1e300, 1e-300, 1e-300), 2), c(1, 0.5, 0.5))
})

test_that("certainty_set() gives the rows taken with certainty and the cut-off", {
  pop <- data.frame(id = 1:12, w = c(1:10, 60, 90))
  for (m in c("systematic", "poisson")) {
    des <- design_weighted("w", n = 4, method = m)
    cs <- certainty_set(pop, des)
    expect_identical(as.vector(cs), which(inclusion_prob(pop, des) == 1))
    expect_identical(as.vector(cs), 11:12)
    expect_identical(attr(cs, "threshold"), 60)
    # The same rows, and the same probabilities, with the cut-off explicit.
    same <- design_certainty("w", attr(cs, "threshold"),
                             rest = design_weighted("w", n = 4 - length(cs),
                                                    method = m))
    expect_equal(inclusion_prob(pop, same), inclusion_prob(pop, des))
  }

  sp <- design_spread("id", n = 4, size = "w")
  expect_identical(as.vector(certainty_set(pop, sp)), 11:12)

  none <- certainty_set(pop, design_weighted("w", n = 1, method = "poisson"))
  expect_length(none, 0)
  expect_identical(attr(none, "threshold"), NA_real_)

  na_pop <- pop
  na_pop$w[12] <- NA
  cs <- certainty_set(na_pop, design_weighted("w", n = 3, method = "poisson",
                                              na_rm = TRUE))
  expect_identical(as.vector(cs), 11L)

  expect_error(certainty_set(pop, design_simple(n = 3)),
               "not proportional to size", fixed = TRUE)
  expect_error(certainty_set(pop, design_weighted("w", n = 3)),
               "method = \"successive\"", fixed = TRUE)
  expect_error(certainty_set(pop, design_spread("id", n = 3)),
               "not proportional to size", fixed = TRUE)
  expect_error(certainty_set(pop, "w"), "must come from one of the design")
  expect_error(certainty_set(pop, design_weighted("w", n = 20,
                                                  method = "systematic")),
               "cannot exceed the number of rows (12)", fixed = TRUE)
})

# ---- within strata -----------------------------------------------------------

strata_pps_pop <- function() {
  data.frame(
    id = 1:60,
    g = rep(c("a", "b", "c"), times = c(30, 20, 10)),
    w = c(1:29, 200, (1:20) * 2, c(5, 5, 6, 6, 7, 7, 8, 8, 9, 9)),
    y = (1:60) %% 7 * 3 + c(1:29, 200, (1:20) * 2,
                            c(5, 5, 6, 6, 7, 7, 8, 8, 9, 9)) * 2
  )
}

test_that("PPS within strata draws each stratum's share of n", {
  pop <- strata_pps_pop()
  d <- design_weighted("w", n = 12, method = "systematic", strata = "g")
  tab <- allocation_table(pop, d)
  expect_identical(tab$n, c(6L, 4L, 2L))
  expect_identical(tab$status, rep("rule", 3))
  expect_true(all(is.na(tab$prob)))
  for (i in 1:20) {
    s <- draw(pop, d, seed = i)
    expect_identical(as.vector(table(factor(s$g, c("a", "b", "c")))),
                     c(6L, 4L, 2L))
  }
  p <- inclusion_prob(pop, d)
  expect_equal(as.vector(tapply(p, pop$g, sum)), c(6, 4, 2))
  # Proportional to size within each stratum, with the heavy row of `a`
  # taken with certainty and the rest of its share spread over the others.
  expect_identical(p[30], 1)
  expect_equal(p[1:29], 5 * pop$w[1:29] / sum(pop$w[1:29]))
  expect_equal(p[31:50], 4 * pop$w[31:50] / sum(pop$w[31:50]))
  expect_identical(as.vector(certainty_set(pop, d)), 30L)
  expect_identical(attr(certainty_set(pop, d), "threshold"), NA_real_)

  eq <- design_weighted("w", n = 12, method = "poisson", strata = "g",
                        allocation = "equal")
  expect_equal(as.vector(tapply(inclusion_prob(pop, eq), pop$g, sum)),
               c(4, 4, 4))
  counts <- design_weighted("w", method = "systematic", strata = "g",
                            allocation = c(a = 2, b = 3, c = 5))
  expect_identical(counts$n, 10L)
  expect_identical(allocation_table(pop, counts)$status, rep("count", 3))
  expect_equal(as.vector(tapply(inclusion_prob(pop, counts), pop$g, sum)),
               c(2, 3, 5))
})

test_that("PPS within strata includes each row as often as it claims", {
  skip_on_cran()
  pop <- strata_pps_pop()[c(1:8, 31:36, 51:55), ]
  for (m in c("systematic", "poisson")) {
    d <- design_weighted("w", n = 7, method = m, strata = "g")
    expect_equal(inclusion_prob(pop, d, simulate = TRUE, R = 2000, seed = 1),
                 inclusion_prob(pop, d), tolerance = 0.06)
  }
})

test_that("Deville within strata is sampling::varest() added over strata", {
  skip_if_not_installed("sampling")
  pop <- strata_pps_pop()
  d <- design_weighted("w", n = 18, method = "systematic", strata = "g")
  for (i in 1:3) {
    s <- draw(pop, d, seed = i, weights = TRUE)
    r <- ht_total(s, "y")
    free <- s$.prob < 1
    ref <- sum(vapply(split(which(free), s$g[free]), function(j) {
      sampling::varest(s$y[j], pik = s$.prob[j])
    }, numeric(1)))
    expect_equal(r$variance, ref)
    expect_identical(r$method, "deville")
    # Rows below certainty less the strata that hold them.
    expect_identical(r$df, sum(free) - 3)
  }
})

test_that("a stratum with one sampled row below certainty is refused by name", {
  pop <- strata_pps_pop()
  d <- design_weighted("w", method = "systematic", strata = "g",
                       allocation = c(a = 4, b = 4, c = 1))
  s <- draw(pop, d, seed = 1, weights = TRUE)
  r <- ht_total(s, "y", variance = "analytic")
  expect_true(is.na(r$se))
  expect_match(r$note, "The stratum `c` has a single sampled row",
               fixed = TRUE)
  expect_match(r$note, "count of 2 or more for each in `allocation`",
               fixed = TRUE)
  # The jackknife refuses it too, rather than leave the stratum out.
  expect_true(is.na(ht_total(s, "y")$se))
  # A domain with no unit in that stratum is measured.
  pop$dom <- ifelse(pop$g == "c", "out", "in")
  s2 <- draw(pop, d, seed = 1, weights = TRUE)
  by <- ht_total(s2, "y", by = "dom", variance = "analytic")
  expect_false(is.na(by$se[by$dom == "in"]))
  expect_true(is.na(by$se[by$dom == "out"]))
})

test_that("the jackknife of PPS within strata deletes rows within strata", {
  pop <- strata_pps_pop()
  s <- draw(pop, design_weighted("w", n = 18, method = "systematic",
                                 strata = "g"), seed = 2, weights = TRUE)
  jk <- ht_total(s, "y", variance = "jackknife")
  expect_match(jk$note, "in 3 strata", fixed = TRUE)
  expect_match(jk$note, "taken with certainty add no variance", fixed = TRUE)
  expect_gt(jk$se, 0)
})

test_that("Poisson within strata keeps the Poisson variance", {
  pop <- strata_pps_pop()
  s <- draw(pop, design_weighted("w", n = 18, method = "poisson",
                                 strata = "g"), seed = 3, weights = TRUE)
  r <- ht_total(s, "y")
  expect_equal(r$variance, sum((1 - s$.prob) / s$.prob^2 * s$y^2))
  free <- s$.prob < 1
  expect_identical(r$df,
                   as.numeric(sum(free) - length(unique(s$g[free]))))
})

test_that("strata need a closed-form PPS method; allocation needs strata", {
  expect_error(design_weighted("w", n = 10, strata = "g"),
               "cannot be combined with method = \"successive\"",
               fixed = TRUE)
  expect_error(design_weighted("w", n = 10, method = "poisson",
                               allocation = "equal"),
               "`strata` is NULL", fixed = TRUE)
  expect_error(design_weighted("w", n = 10, method = "poisson", strata = "g",
                               allocation = "neyman"),
               "`allocation` must be one of", fixed = TRUE)
  expect_error(design_weighted("w", method = "poisson", strata = "g"),
               "`n` is missing", fixed = TRUE)
  expect_error(design_weighted("w", n = 9, method = "poisson", strata = "g",
                               allocation = c(a = 2, b = 3)),
               "not the sum of the counts", fixed = TRUE)
  pop <- strata_pps_pop()
  expect_error(inclusion_prob(pop, design_weighted(
    "w", method = "systematic", strata = "g",
    allocation = c(a = 2, b = 3, c = 11))),
    "Stratum `c` has 10 row(s) but `allocation` asks for 11. Lower its count.",
    fixed = TRUE)
  expect_error(inclusion_prob(pop, design_weighted(
    "w", method = "systematic", strata = "g",
    allocation = c(a = 2, b = 3, c = 5, d = 1))), "`d`", fixed = TRUE)
  expect_error(inclusion_prob(pop, design_weighted(
    "w", method = "systematic", strata = "g",
    allocation = c(a = 2, b = 3))), "`c`", fixed = TRUE)
})

test_that("a missing stratum key follows na_rm, like a missing weight", {
  pop <- strata_pps_pop()
  pop$g[2] <- NA
  d <- design_weighted("w", n = 12, method = "systematic", strata = "g")
  expect_error(draw(pop, d, seed = 1), "1 row(s) have a missing stratum key",
               fixed = TRUE)
  kept <- design_weighted("w", n = 12, method = "systematic", strata = "g",
                          na_rm = TRUE)
  p <- inclusion_prob(pop, kept)
  expect_identical(p[2], 0)
  expect_equal(sum(p), 12)
  expect_identical(allocation_table(pop, kept)$N, c(29L, 20L, 10L))
})

test_that("samples bound stratum by stratum are refused, with strata the way", {
  pop <- strata_pps_pop()
  parts <- lapply(c("a", "b", "c"), function(h) {
    draw(pop[pop$g == h, ], design_weighted("w", n = 4,
                                            method = "systematic"),
         seed = 1, weights = TRUE)
  })
  expect_error(ht_total(do.call(rbind, parts), "y"))
  whole <- draw(pop, design_weighted("w", n = 12, method = "systematic",
                                     strata = "g"), seed = 1, weights = TRUE)
  expect_false(is.na(ht_total(whole, "y")$se))
})
