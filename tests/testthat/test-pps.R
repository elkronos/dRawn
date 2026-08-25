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
