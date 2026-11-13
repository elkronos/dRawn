test_that("plan_size delivers the margin it promises, for a mean", {
  set.seed(3030)
  N <- 20000
  pop <- data.frame(id = 1:N, y = stats::rnorm(N, 100, 40))
  truth <- mean(pop$y)

  p <- plan_size(margin = 5, sd = 40, N = N)
  hit <- replicate(500, {
    s <- draw(pop, design_simple(n = p$n), weights = TRUE)
    abs(ht_mean(s, "y", variance = "none")$mean - truth) <= 5
  })
  expect_gt(mean(hit), 0.92)
  expect_lt(mean(hit), 0.99)
})

test_that("plan_size delivers the margin it promises, for a proportion", {
  set.seed(4040)
  N <- 20000
  pop <- data.frame(id = 1:N, y = stats::rbinom(N, 1, 0.45))
  truth <- mean(pop$y)

  p <- plan_size(margin = 0.03, N = N, target = "proportion")
  hit <- replicate(500, {
    s <- draw(pop, design_simple(n = p$n), weights = TRUE)
    abs(ht_mean(s, "y", variance = "none")$mean - truth) <= 0.03
  })
  expect_gt(mean(hit), 0.92)
})

test_that("the size moves the way the arithmetic says it should", {
  base <- plan_size(margin = 5, sd = 40, N = 20000)$n

  # Halving the margin quadruples the requirement
  expect_gt(plan_size(margin = 2.5, sd = 40, N = 20000)$n, 3 * base)
  # More confidence, more rows
  expect_gt(plan_size(margin = 5, sd = 40, N = 20000, level = 0.99)$n, base)
  # More spread, more rows
  expect_gt(plan_size(margin = 5, sd = 80, N = 20000)$n, base)
  # A design effect scales it directly
  expect_equal(plan_size(margin = 5, sd = 40, N = Inf, deff = 4)$n,
               4 * plan_size(margin = 5, sd = 40, N = Inf)$n, tolerance = 0.01)
  # A finite frame needs fewer than an unbounded one
  expect_lt(base, plan_size(margin = 5, sd = 40)$n)
})

test_that("non-response inflates the draw but not the analysable size", {
  full <- plan_size(margin = 5, sd = 40, N = 20000)
  part <- plan_size(margin = 5, sd = 40, N = 20000, response = 0.6)
  expect_equal(full$n_effective, part$n_effective)
  expect_equal(part$n, ceiling(part$n_effective / 0.6))
  expect_gt(part$n, full$n)
})

test_that("p = 0.5 is the most demanding proportion", {
  at_half <- plan_size(margin = 0.03, N = Inf, target = "proportion")$n
  for (p in c(0.1, 0.3, 0.7, 0.9)) {
    expect_lte(plan_size(margin = 0.03, p = p, N = Inf,
                         target = "proportion")$n, at_half)
  }
})

test_that("a margin on a total is the same problem scaled by N", {
  N <- 20000
  a <- plan_size(margin = 5 * N, sd = 40, N = N, target = "total")
  b <- plan_size(margin = 5, sd = 40, N = N, target = "mean")
  expect_equal(a$n, b$n)
  expect_error(plan_size(margin = 1e5, sd = 40, target = "total"), "finite `N`")
})

test_that("an unreachable margin is capped at a census and says so", {
  p <- plan_size(margin = 0.2, sd = 40, N = 300)
  expect_equal(p$n, 300)
  expect_true(p$capped)
  expect_output(print(p), "census")

  # Non-response cannot push the draw past the frame either.
  q <- plan_size(margin = 0.2, sd = 40, N = 300, response = 0.5)
  expect_equal(q$n, 300)
})

test_that("plan_size validates its inputs", {
  expect_error(plan_size(margin = 0, sd = 10), "positive")
  expect_error(plan_size(margin = -1, sd = 10), "positive")
  expect_error(plan_size(margin = 5), "`sd` is required")
  expect_error(plan_size(margin = 5, sd = 40, level = 0), "between 0 and 1")
  expect_error(plan_size(margin = 5, sd = 40, response = 0), "in \\(0, 1\\]")
  expect_error(plan_size(margin = 5, sd = 40, response = 1.2), "in \\(0, 1\\]")
  expect_error(plan_size(margin = 0.03, p = 2, target = "proportion"),
               "in \\[0, 1\\]")
  expect_error(plan_size(margin = 0.03, p = 0, target = "proportion"),
               "no variation")
  expect_error(plan_size(margin = 5, sd = 40, deff = 0), "positive")
  expect_error(plan_size(margin = 5, sd = 40, N = 0), "positive")
})

test_that("printing a plan states the number and the assumptions behind it", {
  out <- utils::capture.output(
    print(plan_size(margin = 5, sd = 40, N = 20000, deff = 2.5,
                    response = 0.7))
  )
  expect_true(any(grepl("draw", out)))
  expect_true(any(grepl("to analyse", out)))
  expect_true(any(grepl("70% response", out, fixed = TRUE)))
  expect_true(any(grepl("deff 2.5", out, fixed = TRUE)))
})

# ---- sample_summary --------------------------------------------------------

test_that("sample_summary counts what was drawn against what was there", {
  pop <- data.frame(
    id = 1:400,
    site = rep(c("a", "b", "c", "d"), times = c(200, 100, 60, 40))
  )
  s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
  x <- sample_summary(s)

  expect_equal(x$n, nrow(s))
  expect_equal(x$N, 400)
  expect_equal(x$design, "stratified")
  expect_equal(sum(x$by_group$drawn), nrow(s))
  expect_equal(x$by_group$in_frame, c(200, 100, 60, 40))
  expect_equal(x$by_group$rate, x$by_group$drawn / x$by_group$in_frame,
               tolerance = 1e-4)
})

test_that("the weight summary matches the weights on the sample", {
  pop <- data.frame(id = 1:300, w = stats::setNames(rep(c(1, 5), 150), NULL))
  s <- draw(pop, design_weighted("w", n = 50, method = "poisson"), seed = 2,
            weights = TRUE)
  x <- sample_summary(s)
  expect_equal(x$weight_range, range(s$.weight))
  expect_equal(x$weight_cv, stats::sd(s$.weight) / mean(s$.weight))
})

test_that("rows the design can never reach are counted", {
  # Forty rows fall outside the sampling window, so nothing can draw them.
  pop <- data.frame(
    id = 1:200,
    ts = seq(as.POSIXct("2024-01-01", tz = "UTC"), by = "hour",
             length.out = 200)
  )
  des <- design_temporal("ts", from = "2024-01-01", to = "2024-01-07",
                         interval = 1, per_interval = 2, unit = "days")
  s <- draw(pop, des, seed = 1, weights = TRUE)
  expect_equal(sample_summary(s)$unreachable, 200 - 24 * 6)
  expect_output(print(sample_summary(s)), "unreachable rows")
})

test_that("the sampling fraction is over the rows the design can reach", {
  pop <- data.frame(
    id = 1:200,
    ts = seq(as.POSIXct("2024-01-01", tz = "UTC"), by = "hour",
             length.out = 200)
  )
  des <- design_temporal("ts", from = "2024-01-01", to = "2024-01-07",
                         interval = 1, per_interval = 2, unit = "days")
  x <- sample_summary(draw(pop, des, seed = 1, weights = TRUE))
  expect_equal(x$N, 200)
  expect_equal(x$N_reachable, 144)
  expect_output(print(x), "sampling fraction  0.0833  (of the 144 reachable",
                fixed = TRUE)
})

test_that("a design with no grouping reports no group table", {
  pop <- data.frame(id = 1:100)
  s <- draw(pop, design_simple(n = 20), seed = 1, weights = TRUE)
  x <- sample_summary(s)
  expect_null(x$by_group)
  expect_equal(x$unreachable, 0)
  expect_output(print(x), "sampling fraction")
})

test_that("sample_summary reports on a cluster design by cluster", {
  pop <- data.frame(id = 1:200, cl = rep(paste0("c", 1:20), each = 10))
  s <- draw(pop, design_cluster("cl", n_clusters = 4), seed = 1, weights = TRUE)
  x <- sample_summary(s)
  expect_equal(nrow(x$by_group), 4)
  expect_true(all(x$by_group$rate == 1))
})

test_that("sample_summary counts groups from the frame, not the sample", {
  pop <- data.frame(id = 1:200, cl = rep(paste0("c", 1:20), each = 10),
                    site = rep(c("a", "b"), times = c(150, 50)))
  s <- draw(pop, design_cluster("cl", n_clusters = 4), seed = 1, weights = TRUE)
  want <- sample_summary(s)$by_group
  gone <- s
  gone$cl <- NULL
  expect_equal(sample_summary(gone)$by_group, want)
  st <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
  recoded <- st
  recoded$site <- "z"
  expect_equal(sample_summary(recoded)$by_group,
               sample_summary(st)$by_group)
})

test_that("sample_summary needs a sample carrying its design", {
  pop <- data.frame(id = 1:100)
  expect_error(sample_summary(draw(pop, design_simple(n = 10), seed = 1)),
               "weights = TRUE")
  expect_error(sample_summary(pop), "weights = TRUE")
})

test_that("printing a summary shows the group table", {
  pop <- data.frame(
    id = 1:400,
    site = rep(c("a", "b", "c", "d"), times = c(200, 100, 60, 40))
  )
  s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
  out <- utils::capture.output(print(sample_summary(s)))
  expect_true(any(grepl("Sample of 40 from 400", out)))
  expect_true(any(grepl("by site", out)))
  expect_true(any(grepl("in frame", out)))
})

test_that("the frame caps the analysable size, not just the draw", {
  # Capping n at N while leaving n_effective alone promises more rows back
  # than the cap can possibly deliver.
  p <- plan_size(margin = 1.4, sd = 40, N = 600, response = 0.5)
  expect_equal(p$n, 600)
  expect_lte(p$n_effective, floor(600 * 0.5))

  # This frame *can* reach the margin -- 504 rows would do it. Only the
  # response rate puts it out of reach, which is not a census.
  expect_false(p$capped)
  expect_true(p$short)
  expect_equal(p$n_needed, plan_size(margin = 1.4, sd = 40, N = 600)$n)
  expect_output(print(p), "would need more rows than the frame holds")
  expect_failure(expect_output(print(p), "census"))

  # A margin the frame cannot reach however you sample it is a census, and
  # says so even when the response rate is short too.
  cen <- plan_size(margin = 0.2, sd = 40, N = 300, response = 0.5)
  expect_true(cen$capped)
  expect_output(print(cen), "census")
  expect_output(print(cen), "fall short of the margin")

  # Uncapped, the relationship is the plain one
  q <- plan_size(margin = 5, sd = 40, N = 20000, response = 0.5)
  expect_lt(q$n, 20000)
  expect_false(q$capped)
  expect_equal(q$n, ceiling(q$n_effective / 0.5))
})

test_that("sample_summary refuses a subset and explains a lost design", {
  pop <- data.frame(id = 1:400,
                    site = rep(c("a", "b", "c", "d"), times = c(200, 100, 60, 40)))
  s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
  expect_equal(sample_summary(s[order(-s$id), ]), sample_summary(s))
  expect_error(sample_summary(s[s$site == "a", ]), "rows that were drawn")
  merged <- merge(s, data.frame(id = s$id, resp = 1), by = "id")
  expect_error(sample_summary(merged), "does not carry its design")
  expect_error(sample_summary(merged), "transform()", fixed = TRUE)
})

test_that("a plan prints a round margin in full", {
  out <- utils::capture.output(print(plan_size(margin = 1e5, sd = 40,
                                               N = 20000, target = "total")))
  expect_true(any(grepl("+/- 100,000 at 95% confidence", out, fixed = TRUE)))
  expect_false(any(grepl("e+", out, fixed = TRUE)))
})

test_that("sample_summary lists a stratum allocated no rows", {
  inv <- data.frame(id = 1:600, site = rep(c("north", "south", "east", "west"),
                                           times = c(300, 180, 90, 30)))
  expect_warning(
    s <- draw(inv, design_stratified("site", n = 12), seed = 1,
              weights = TRUE),
    class = "drawn_empty_stratum")
  x <- sample_summary(s)
  expect_equal(x$by_group$group, c("east", "north", "south", "west"))
  expect_equal(x$by_group$drawn, c(2, 6, 4, 0))
  expect_equal(x$by_group$in_frame, c(90, 300, 180, 30))
  out <- capture.output(print(x))
  # Every rate has three decimals, so the column lines up.
  expect_true(any(grepl("east +2 +90  0.022$", out)))
  expect_true(any(grepl("north +6 +300  0.020$", out)))
  expect_true(any(grepl("west +0 +30  0.000$", out)))
})

test_that("sample_summary gives a certainty design the strata of its rest", {
  pop <- data.frame(id = 1:200, site = rep(c("a", "b"), times = c(120, 80)),
                    w = c(rep(500, 3), rep(1, 197)))
  s <- draw(pop, design_certainty("w", 100,
                                  design_stratified("site", n = 20)),
            seed = 1, weights = TRUE)
  x <- sample_summary(s)
  expect_equal(x$group_col, "site")
  expect_equal(x$by_group$group, c("(certain)", "a", "b"))
  expect_equal(x$by_group$drawn, c(3, 12, 8))
  expect_equal(x$by_group$in_frame, c(3, 117, 80))
  expect_equal(x$n_certain, 3L)
  expect_equal(x$certain_by, c(threshold = 3L, size = 0L, whole = 0L))
  out <- capture.output(print(x))
  expect_true(any(grepl("certain rows       3  (3 above the threshold)", out,
                        fixed = TRUE)))
  expect_true(any(grepl("(certain)", out, fixed = TRUE)))
})

test_that("sample_summary says why rows were taken with certainty", {
  big <- data.frame(id = 1:100, w = c(rep(1, 95), rep(1000, 5)))
  s <- draw(big, design_weighted("w", n = 20, method = "systematic"), seed = 1,
            weights = TRUE)
  x <- sample_summary(s)
  expect_equal(x$n_certain, 5L)
  expect_equal(x$certain_by, c(threshold = 0L, size = 5L, whole = 0L))
  expect_output(print(x), "certain rows       5  (5 by size)", fixed = TRUE)

  # A stratum taken whole is certain too, and neither of the above.
  st <- data.frame(id = 1:53, site = rep(c("big", "small"), c(50, 3)))
  s <- draw(st, design_stratified("site", allocation = c(big = 10, small = 3)),
            seed = 1, weights = TRUE)
  x <- sample_summary(s)
  expect_equal(x$certain_by, c(threshold = 0L, size = 0L, whole = 3L))

  srs <- draw(st, design_simple(n = 10), seed = 1, weights = TRUE)
  x <- sample_summary(srs)
  expect_equal(x$n_certain, 0L)
  expect_false(any(grepl("certain rows", capture.output(print(x)))))
})

test_that("sample_summary gives the expected size of a sample of random size", {
  pop <- data.frame(id = 1:240)
  s <- draw(pop, design_systematic(interval = 7), seed = 1, weights = TRUE)
  x <- sample_summary(s)
  expect_equal(x$expected_n, 240 / 7)
  expect_output(print(x), "expected n         34.29", fixed = TRUE)

  fixed <- sample_summary(draw(pop, design_simple(n = 30), seed = 1,
                               weights = TRUE))
  expect_equal(fixed$expected_n, 30)
  expect_false(any(grepl("expected n", capture.output(print(fixed)))))
})

test_that("sample_summary reports Kish's effective sample size", {
  pop <- data.frame(id = 1:300, w = 1 + (1:300 * 13) %% 40)
  s <- draw(pop, design_weighted("w", n = 40, method = "systematic"), seed = 1,
            weights = TRUE)
  x <- sample_summary(s)
  expect_equal(x$n_kish, sum(s$.weight)^2 / sum(s$.weight^2))
  expect_lt(x$n_kish, nrow(s))
  expect_output(print(x), "effective n")
  srs <- sample_summary(draw(pop, design_simple(n = 30), seed = 1,
                             weights = TRUE))
  expect_equal(srs$n_kish, 30)
})

test_that("sample_summary groups a temporal sample by interval", {
  pop <- data.frame(
    id = 1:300,
    ts = seq(as.POSIXct("2024-01-01", tz = "UTC"), by = "hour",
             length.out = 300)
  )
  des <- design_temporal("ts", from = "2024-01-01", to = "2024-01-13",
                         interval = 1, per_interval = 2, unit = "days")
  x <- sample_summary(draw(pop, des, seed = 1, weights = TRUE))
  expect_equal(x$group_col, "ts")
  # Twelve intervals named by their start, in time order.
  expect_equal(x$by_group$group,
               paste(format(as.Date("2024-01-01") + 0:11), "00:00"))
  expect_equal(x$by_group$drawn, rep(2, 12))
  expect_equal(x$by_group$in_frame, rep(24, 12))
})
