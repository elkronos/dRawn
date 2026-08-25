test_that("every stratum is represented and the total is exact", {
  df <- make_df()
  res <- draw(df, design_stratified("site", n = 40), seed = 123)
  expect_equal(nrow(res), 40)
  expect_true(all(table(res$site) >= 1))
})

test_that("equal allocation splits the total evenly", {
  df <- make_df()
  res <- draw(df, design_stratified("site", n = 20, allocation = "equal"),
              seed = 123)
  expect_equal(as.vector(table(res$site)), rep(5L, 4))
})

test_that("proportional allocation hits the requested total exactly", {
  uneven <- data.frame(
    x = 1:143,
    g = rep(letters[1:7], times = c(31, 29, 23, 19, 17, 13, 11))
  )
  expect_equal(nrow(draw(uneven, design_stratified("g", n = 50), seed = 1)), 50)
  expect_equal(nrow(draw(uneven, design_stratified("g", n = 33), seed = 1)), 33)

  # Independent per-stratum round() is subject to banker's rounding, which
  # takes four strata of 25 down to 8 rows for a request of 10.
  even <- data.frame(x = 1:100, g = rep(c("a", "b", "c", "d"), each = 25))
  expect_equal(nrow(draw(even, design_stratified("g", n = 10), seed = 1)), 10)
})

test_that("rare strata are covered only when asked for", {
  skewed <- data.frame(x = 1:1000, g = c(rep("big", 999), "rare"))

  res <- draw(skewed, design_stratified("g", n = 10), seed = 1)
  expect_equal(nrow(res), 10)

  covered <- draw(skewed, design_stratified("g", n = 10, min_per_stratum = 1),
                  seed = 1)
  expect_equal(nrow(covered), 10)
  expect_true("rare" %in% covered$g)
})

test_that("multiple strata columns are cross-classified", {
  df <- data.frame(x = 1:40,
                   g = rep(c("a", "b"), 20),
                   h = rep(c("p", "q"), each = 20))
  res <- draw(df, design_stratified(c("g", "h"), n = 12), seed = 1)
  expect_equal(nrow(res), 12)
  expect_equal(length(unique(paste(res$g, res$h))), 4)
})

test_that("missing columns and NA strata are reported clearly", {
  df <- make_df()
  expect_error(draw(df, design_stratified("nope", n = 10), seed = 1),
               "missing column")

  na_df <- data.frame(x = 1:30, g = c(rep("a", 20), rep(NA, 10)))
  expect_error(draw(na_df, design_stratified("g", n = 5), seed = 1),
               "missing stratum key")
  expect_equal(
    nrow(draw(na_df, design_stratified("g", n = 5, na_rm = TRUE), seed = 1)),
    5
  )
})

test_that("over-drawing is refused without replacement and allowed with it", {
  df <- make_df()
  expect_error(draw(df, design_stratified("site", n = 200), seed = 1),
               "cannot exceed the number of rows")
  expect_equal(
    nrow(draw(df, design_stratified("site", n = 200, replace = TRUE), seed = 1)),
    200
  )
})

test_that("min_per_stratum is checked against n", {
  df <- make_df()
  expect_error(
    draw(df, design_stratified("site", n = 2, min_per_stratum = 1), seed = 1),
    "needs at least 4 rows"
  )
})

test_that("strata must be given as column names", {
  expect_error(design_stratified(1, n = 5), "non-empty column names")
  expect_error(design_stratified(character(0), n = 5), "non-empty column names")
})

test_that("min_per_stratum still totals n exactly, at every feasible floor", {
  # A single-pass reduction can only take back one row per stratum, which
  # leaves 62 rows for a request of 60 with four strata and a floor of 8.
  d <- data.frame(id = 1:600,
                  site = rep(c("north", "south", "east", "west"),
                             times = c(300, 180, 90, 30)))
  for (m in 0:20) {
    for (nn in c(4L, 10L, 33L, 60L, 137L, 300L)) {
      if (m * 4 > nn) next
      expect_equal(nrow(draw(d, design_stratified("site", n = nn,
                                                  min_per_stratum = m), seed = 1)),
                   nn, info = paste("min_per_stratum =", m, "n =", nn))
    }
  }
})

test_that("Neyman allocation puts rows where the values vary most", {
  set.seed(2)
  d <- data.frame(id = 1:300, g = rep(c("a", "b", "c"), each = 100))
  d$v <- c(stats::rnorm(100, 50, 1), stats::rnorm(100, 50, 10),
           stats::rnorm(100, 50, 30))

  ney <- table(draw(d, design_stratified("g", n = 60, allocation = "neyman",
                                         allocation_by = "v"), seed = 1)$g)
  expect_equal(sum(ney), 60)
  expect_lt(ney[["a"]], ney[["b"]])
  expect_lt(ney[["b"]], ney[["c"]])

  # Equal-sized strata with equal spread reduce to proportional
  flat <- data.frame(id = 1:300, g = rep(c("a", "b", "c"), each = 100),
                     v = rep(stats::rnorm(100), 3))
  expect_equal(as.vector(table(draw(flat, design_stratified(
    "g", n = 60, allocation = "neyman", allocation_by = "v"), seed = 1)$g)),
    rep(20L, 3))
})

test_that("Neyman validates its auxiliary column", {
  d <- data.frame(id = 1:30, g = rep(c("a", "b"), 15), v = 1:30,
                  ch = letters[rep(1:2, 15)])
  expect_error(design_stratified("g", n = 10, allocation = "neyman"),
               "needs `allocation_by`")
  expect_error(draw(d, design_stratified("g", n = 10, allocation = "neyman",
                                         allocation_by = "ch"), seed = 1),
               "must be numeric")
  na_d <- d; na_d$v[1] <- NA
  expect_error(draw(na_d, design_stratified("g", n = 10, allocation = "neyman",
                                            allocation_by = "v"), seed = 1),
               "are missing")
})

test_that("Neyman inclusion probabilities are consistent with the draw", {
  set.seed(4)
  d <- data.frame(id = 1:200, g = rep(c("a", "b"), each = 100))
  d$v <- c(stats::rnorm(100, 0, 1), stats::rnorm(100, 0, 20))
  des <- design_stratified("g", n = 40, allocation = "neyman", allocation_by = "v")

  p <- inclusion_prob(d, des)
  hits <- integer(200)
  for (i in 1:1500) {
    ids <- draw(d, des, seed = i)$id
    hits[ids] <- hits[ids] + 1L
  }
  expect_equal(hits / 1500, p, tolerance = 0.05)
})
