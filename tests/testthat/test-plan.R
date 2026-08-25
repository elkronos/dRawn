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
