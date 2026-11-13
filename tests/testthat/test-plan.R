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

test_that("plan_size gives the smallest size that meets the margin", {
  z <- stats::qnorm(0.975)
  # The smallest n whose variance meets (margin / z)^2, found by trying every
  # size the frame allows, and never below two rows.
  smallest <- function(N, ratio) {
    n <- seq_len(N)
    ok <- ratio(n) <= 1 + 1e-12
    n <- if (any(ok)) min(n[ok]) else N
    if (N >= 2) max(n, 2) else n
  }
  for (N in c(3, 7, 10, 25, 60)) {
    for (m in c(5, 8.4, 12.5, 20, 40)) {
      mean_n <- smallest(N, function(n) (1 / n - 1 / N) * 40^2 / (m / z)^2)
      expect_equal(plan_size(margin = m, sd = 40, N = N)$n, mean_n,
                   info = paste("mean", N, m))
      total_n <- smallest(N, function(n) {
        N^2 * (1 / n - 1 / N) * 2 * 40^2 / (m * N / z)^2
      })
      expect_equal(plan_size(margin = m * N, sd = 40, N = N, deff = 2,
                             target = "total")$n, total_n,
                   info = paste("total", N, m))
    }
    for (m in c(0.05, 0.1, 0.2, 0.4)) {
      for (p in c(0.1, 0.5)) {
        prop_n <- smallest(N, function(n) {
          p * (1 - p) / n * (N - n) / (N - 1) / (m / z)^2
        })
        expect_equal(plan_size(margin = m, p = p, N = N,
                               target = "proportion")$n, prop_n,
                     info = paste("proportion", N, m, p))
      }
    }
  }
})

test_that("a margin that fewer rows than the frame meet is not a census", {
  # With sd = 40, nine of ten rows meet a margin of 8.4.
  z <- stats::qnorm(0.975)
  expect_lte((1 / 9 - 1 / 10) * 40^2, (8.4 / z)^2)
  p <- plan_size(margin = 8.4, sd = 40, N = 10)
  expect_equal(p$n, 9)
  expect_false(p$capped)

  # Nine rows do not meet a margin of 2, so only a census does.
  expect_gt((1 / 9 - 1 / 10) * 40^2, (2 / z)^2)
  cen <- plan_size(margin = 2, sd = 40, N = 10)
  expect_equal(cen$n, 10)
  expect_true(cen$capped)
  expect_output(print(cen), "No sample smaller than the frame meets")
})

test_that("a census plan is capped and not short", {
  cen <- plan_size(margin = 0.2, sd = 40, N = 300, response = 0.5)
  expect_true(cen$capped)
  expect_false(cen$short)
  expect_equal(cen$n, 300)
  expect_equal(cen$n_effective, 150)
})

test_that("non-response does not add a row through rounding", {
  # 21 / 0.7 is 30.000000000000004 in floating point.
  p <- plan_size(margin = 17.11, sd = 40, response = 0.7)
  expect_equal(p$n_effective, 21)
  expect_equal(p$n, 30)

  # Thirty rows at 70% response bring back the 21 needed, so the frame is
  # not short.
  q <- plan_size(margin = 10, sd = 40, N = 30, response = 0.7)
  expect_equal(q$n, 30)
  expect_equal(q$n_effective, 21)
  expect_false(q$short)
  expect_false(q$capped)
  expect_false(any(grepl("more rows than the frame holds",
                         utils::capture.output(print(q)), fixed = TRUE)))
})

test_that("a fractional N gives a whole n that draw() accepts", {
  p <- plan_size(margin = 0.01, sd = 40, N = 240.5)
  expect_equal(p$n, 240)
  expect_true(p$capped)
  s <- draw(data.frame(id = 1:240), design_simple(n = p$n), seed = 1)
  expect_equal(nrow(s), 240)
  expect_output(print(p), "N 240.5", fixed = TRUE)

  # An estimated N plans a total as a whole one would, below the cap.
  expect_equal(plan_size(margin = 5000, sd = 40, N = 2400.4,
                         target = "total")$n,
               ceiling(2400.4 / (1 + 2400.4 /
                                   (stats::qnorm(0.975) * 40 * 2400.4 /
                                      5000)^2)))

  expect_error(plan_size(margin = 1000, sd = 40, N = 0.5),
               "`N` is 0.5, and a frame holds at least one row")
  expect_error(plan_size(margin = 5, sd = 40, N = 1e-9), "at least one row")
  expect_equal(plan_size(margin = 5, sd = 40, N = 1)$n, 1)
})

test_that("a size too large for a number gives a census or a clear error", {
  p <- plan_size(margin = 1e-300, sd = 40, N = 240)
  expect_equal(p$n, 240)
  expect_true(p$capped)
  expect_equal(plan_size(margin = 5, sd = 40, N = 240, deff = 1e308)$n, 240)
  expect_equal(plan_size(margin = 5, sd = 1e300, N = 240)$n, 240)
  expect_equal(plan_size(margin = 5, sd = 40, N = 240, level = 1 - 1e-16)$n,
               240)
  expect_equal(plan_size(margin = 0.03, N = 240, level = 1 - 1e-16,
                         target = "proportion")$n, 240)
  expect_error(plan_size(margin = 1e-300, sd = 40),
               "`margin` = 1e-300 with a spread of 40")
  expect_error(plan_size(margin = 1e-300, sd = 40), "finite `N`")
  expect_error(plan_size(margin = 5, sd = 40, level = 1 - 1e-16),
               "lower `level`")
})

test_that("a plan has two rows at least, the fewest that give a variance", {
  p <- plan_size(margin = 80, sd = 40)
  expect_equal(p$n, 2)
  expect_equal(p$n_effective, 2)
  expect_true(p$raised)
  expect_output(print(p), "interval need\\s+two")
  expect_false(plan_size(margin = 5, sd = 40)$raised)

  # A proportion margin close to 1 would have asked for one row.
  expect_equal(plan_size(margin = 0.99, target = "proportion")$n, 2)

  two <- plan_size(margin = 80, sd = 40, N = 2)
  expect_equal(two$n, 2)
  expect_false(two$capped)
  expect_equal(plan_size(margin = 80, sd = 40, N = 1)$n, 1)
  expect_equal(plan_size(margin = 80, sd = 40, response = 0.5)$n, 4)
})

test_that("an argument the target does not use is ignored with a warning", {
  expect_warning(
    a <- plan_size(margin = 0.05, sd = 0.1, target = "proportion"),
    class = "drawn_irrelevant_argument")
  expect_equal(a$n, plan_size(margin = 0.05, target = "proportion")$n)
  expect_warning(plan_size(margin = 0.05, sd = 0.1, target = "proportion"),
                 "`sd` is not used")

  expect_warning(b <- plan_size(margin = 5, sd = 40, p = 7),
                 class = "drawn_irrelevant_argument")
  expect_equal(b$n, plan_size(margin = 5, sd = 40)$n)
  expect_warning(plan_size(margin = 5e4, sd = 40, p = 0.2, N = 2e4,
                           target = "total"),
                 "`p` is not used for target = \"total\"")

  expect_error(plan_size(margin = 0.03, p = 0.2),
               "use target = \"proportion\"", fixed = TRUE)
  expect_error(plan_size(margin = 5),
               "plausible range divided by four.$")

  expect_silent(plan_size(margin = 0.03, p = 0.2, target = "proportion"))
  expect_silent(plan_size(margin = 5, sd = 40))
})

test_that("a proportion margin of 1 or more is refused", {
  expect_error(plan_size(margin = 5, N = 20000, target = "proportion"),
               "such as 0.05 for 5 percentage points", fixed = TRUE)
  expect_error(plan_size(margin = 1, target = "proportion"),
               "0.01 for 1 percentage point.", fixed = TRUE)
  expect_equal(plan_size(margin = 0.03, target = "proportion")$n,
               plan_size(margin = 0.03, p = 0.5, target = "proportion")$n)
})

test_that("a proportion plan prints and returns p", {
  p <- plan_size(margin = 0.03, p = 0.1, N = 20000, target = "proportion")
  expect_equal(p$p, 0.1)
  expect_output(print(p), "assuming       p 0.1, N 20,000", fixed = TRUE)
  m <- plan_size(margin = 5, sd = 40, N = 20000)
  expect_true(is.na(m$p))
  expect_output(print(m), "assuming       sd 40, N 20,000", fixed = TRUE)
})

test_that("deff takes the estimate it came from, of the kind target plans", {
  pop <- data.frame(id = 1:400, cl = rep(paste0("c", 1:40), each = 10))
  pop$y <- rep(seq(10, 49), each = 10) + rep(c(-2, 2), 200)
  pilot <- draw(pop, design_cluster("cl", n_clusters = 6), seed = 1,
                weights = TRUE)
  m <- ht_mean(pilot, "y")
  tot <- ht_total(pilot, "y")
  expect_true(is.finite(deff(m)))
  expect_true(is.finite(deff(tot)))

  expect_equal(plan_size(margin = 3, sd = 10, N = 400, deff = m),
               plan_size(margin = 3, sd = 10, N = 400, deff = deff(m)))
  expect_equal(plan_size(margin = 3, sd = 10, N = 400, deff = m)$deff,
               deff(m))
  expect_equal(plan_size(margin = 1200, sd = 10, N = 400, deff = tot,
                         target = "total")$deff, deff(tot))
  expect_equal(plan_size(margin = 0.1, N = 400, deff = m,
                         target = "proportion")$deff, deff(m))

  expect_error(plan_size(margin = 1200, sd = 10, N = 400, deff = m,
                         target = "total"),
               "Pass ht_total() on the same pilot", fixed = TRUE)
  expect_error(plan_size(margin = 3, sd = 10, N = 400, deff = tot),
               "Pass ht_mean() on the same pilot", fixed = TRUE)
  expect_error(plan_size(margin = 0.1, N = 400, deff = tot,
                         target = "proportion"),
               "`target` is \"proportion\"", fixed = TRUE)

  none <- ht_mean(pilot, "y", variance = "none")
  expect_error(plan_size(margin = 3, sd = 10, N = 400, deff = none),
               "design effect of NA")
  expect_error(plan_size(margin = 3, sd = 10, deff = NA),
               "`deff` is NA, which deff() gives", fixed = TRUE)
  expect_error(plan_size(margin = 3, sd = 10, deff = NA_real_),
               "`deff` is NA", fixed = TRUE)
  expect_error(plan_size(margin = 3, sd = 10,
                         deff = ht_mean(pilot, "y", by = "cl")),
               "without `by`", fixed = TRUE)
})

test_that("plan_size meets the margin at a small N, where the correction bites", {
  skip_on_cran()
  N <- 60
  pop <- data.frame(id = seq_len(N),
                    y = 100 + 40 * stats::qnorm(stats::ppoints(N)))
  truth <- mean(pop$y)
  p <- plan_size(margin = 8, sd = stats::sd(pop$y), N = N)
  # Without the correction the plan would be more than the frame.
  expect_gt(plan_size(margin = 8, sd = stats::sd(pop$y))$n, N)
  expect_lt(p$n, N)
  hit <- vapply(seq_len(500), function(i) {
    s <- draw(pop, design_simple(n = p$n), seed = i, weights = TRUE)
    abs(ht_mean(s, "y", variance = "none")$mean - truth) <= 8
  }, logical(1))
  expect_gt(mean(hit), 0.92)
  expect_lt(mean(hit), 0.99)
})

test_that("plan_size meets a margin on a total", {
  skip_on_cran()
  N <- 400
  pop <- data.frame(id = seq_len(N),
                    y = exp(3 + 0.5 * stats::qnorm(stats::ppoints(N))))
  truth <- sum(pop$y)
  M <- 0.04 * truth
  p <- plan_size(margin = M, sd = stats::sd(pop$y), N = N, target = "total")
  expect_lt(p$n, N)
  hit <- vapply(seq_len(500), function(i) {
    s <- draw(pop, design_simple(n = p$n), seed = 1000 + i, weights = TRUE)
    abs(ht_total(s, "y", variance = "none")$total - truth) <= M
  }, logical(1))
  expect_gt(mean(hit), 0.92)
  expect_lt(mean(hit), 0.99)
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
