test_that("weights = TRUE prepends .prob, .weight and .row", {
  df <- data.frame(id = 1:100, site = rep(letters[1:4], each = 25))
  s <- draw(df, design_stratified("site", n = 12), seed = 1, weights = TRUE)
  expect_identical(names(s), c(".prob", ".weight", ".row", names(df)))
  expect_identical(s$.row, s$id)
  expect_equal(s$.weight, 1 / s$.prob)
  expect_true(all(s$.prob > 0 & s$.prob <= 1))
})

test_that("the probabilities attached match inclusion_prob() on the population", {
  df <- data.frame(id = 1:20, g = rep(c("a", "b"), times = c(15, 5)))
  d <- design_stratified("g", n = 8)
  s <- draw(df, d, seed = 1, weights = TRUE)
  expected <- inclusion_prob(df, d)[s$id]
  expect_equal(s$.prob, expected)
})

test_that("weights are attached correctly even with duplicate rows", {
  # Two identical rows must still get their own probabilities rather than
  # being matched to the same population row.
  df <- data.frame(id = c(1, 1, 2, 2, 3, 3), g = rep(c("a", "b", "c"), each = 2))
  s <- draw(df, design_simple(n = 4), seed = 1, weights = TRUE)
  expect_equal(nrow(s), 4)
  expect_true(all(s$.prob == 4/6))
})

test_that("with-replacement designs refuse weights rather than mislead", {
  # `.prob` for a with-replacement design is the chance of appearing at least
  # once, but the sample holds duplicates. Summing y * .weight over it counts
  # every duplicate at the distinct-unit weight, high by a factor of
  # (n / N) / (1 - (1 - 1 / N)^n) -- 169% at n = 25, N = 10 -- and ht_total()
  # inherited the same bias.
  df <- data.frame(id = 1:10, y = 1:10)
  for (des in list(design_simple(n = 25, replace = TRUE),
                   design_stratified("g", n = 6, replace = TRUE))) {
    expect_error(
      draw(transform(df, g = rep(c("a", "b"), each = 5)), des, seed = 1,
           weights = TRUE),
      "with-replacement")
  }
  # Without replacement, the same designs weight fine.
  s <- draw(df, design_simple(n = 4), seed = 1, weights = TRUE)
  expect_false(anyNA(s$.prob))
  expect_equal(sum(s$.weight), nrow(df))
})

test_that("weights = TRUE fails loudly where the design has no closed form", {
  df <- data.frame(id = 1:20, w = 1:20)
  expect_error(draw(df, design_weighted("w", n = 5), seed = 1, weights = TRUE),
               "no closed-form inclusion probability")
  expect_error(draw(df, design_bootstrap(n_replicates = 3), seed = 1,
                    weights = TRUE),
               "not a probability sample")
})

test_that("weights = TRUE refuses to clobber existing columns", {
  df <- data.frame(id = 1:10, .prob = 1)
  expect_error(draw(df, design_simple(n = 3), seed = 1, weights = TRUE),
               "already has `.prob`")

  df2 <- data.frame(id = 1:10, .drawn_row_id = 1)
  expect_error(draw(df2, design_simple(n = 3), seed = 1, weights = TRUE),
               "already has a column called")

  df3 <- data.frame(id = 1:10, .row = 1, .weight = 2)
  expect_error(draw(df3, design_simple(n = 3), seed = 1, weights = TRUE),
               "already has `.weight` and `.row`; weights = TRUE would overwrite them",
               fixed = TRUE)
})

test_that("weights = TRUE needs a data frame", {
  expect_error(draw(as.list(1:10), design_reservoir(n = 3), seed = 1,
                    weights = TRUE),
               "needs a data frame")
})

test_that("weights = TRUE does not change which rows are drawn", {
  df <- data.frame(id = 1:100, site = rep(letters[1:4], each = 25))
  d <- design_stratified("site", n = 12)
  plain <- draw(df, d, seed = 42)
  with_w <- draw(df, d, seed = 42, weights = TRUE)
  expect_equal(with_w$id, plain$id)
})

test_that("Horvitz-Thompson totals are unbiased across designs", {
  set.seed(11)
  pop <- data.frame(id = 1:40,
                    g = rep(c("a", "b"), times = c(30, 10)),
                    cl = rep(letters[1:8], each = 5),
                    y = round(stats::runif(40, 1, 100)))
  pop$w <- pop$y + 20
  true <- sum(pop$y)

  designs <- list(
    simple      = design_simple(n = 10),
    stratified  = design_stratified("g", n = 12),
    cluster     = design_cluster("cl", n_clusters = 3),
    multistage  = design_multistage("cl", n_clusters = 4, n = 8),
    systematic  = design_systematic(interval = 4),
    pps_syst    = design_weighted("w", n = 10, method = "systematic"),
    pps_poisson = design_weighted("w", n = 10, method = "poisson")
  )

  for (nm in names(designs)) {
    est <- vapply(1:600, function(i) {
      s <- draw(pop, designs[[nm]], seed = i, weights = TRUE)
      if (nrow(s) == 0L) return(0)
      sum(s$y * s$.weight)
    }, numeric(1))
    se <- stats::sd(est) / sqrt(length(est))
    expect_lt(abs(mean(est) - true) / se, 4)
  }
})

test_that("a certainty design over a with-replacement rest refuses to weight", {
  # A certainty design has no `replace` of its own, so checking only the
  # outer design let duplicates through: the total came out about a third
  # high on a 240-row frame.
  df <- data.frame(id = 1:40, g = rep(c("a", "b"), each = 20),
                   w = c(rep(1, 36), rep(500, 4)), y = 1:40)
  for (rest in list(design_simple(n = 30, replace = TRUE),
                    design_stratified("g", n = 30, replace = TRUE),
                    design_multistage("g", n_clusters = 2, n = 30,
                                      replace = TRUE))) {
    d <- design_certainty("w", 100, rest)
    expect_error(draw(df, d, seed = 1, weights = TRUE),
                 "certainty design over a with-replacement `rest`",
                 fixed = TRUE)
    # Without weights the draw is legitimate and still allowed.
    expect_equal(nrow(draw(df, d, seed = 1)), 34)
  }
  s <- draw(df, design_certainty("w", 100, design_simple(n = 10)), seed = 1,
            weights = TRUE)
  expect_equal(nrow(s), 14)
})

test_that("the with-replacement refusal gives the bias factor, not one figure", {
  df <- data.frame(id = 1:10, g = rep(c("a", "b"), each = 5))
  for (des in list(design_simple(n = 25, replace = TRUE),
                   design_stratified("g", n = 6, replace = TRUE))) {
    msg <- tryCatch(draw(df, des, seed = 1, weights = TRUE),
                    error = conditionMessage)
    expect_match(msg, "(n / N) / (1 - (1 - 1 / N)^n)", fixed = TRUE)
    expect_match(msg, "N_h / n_h * sum(y_h)", fixed = TRUE)
    expect_no_match(msg, "15%\\s*high")
  }
})
