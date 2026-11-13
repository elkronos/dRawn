# A data frame column can itself be a list or a matrix. Neither can serve as a
# grouping or ordering key, and before this check two designs failed silently
# rather than loudly.

odd_df <- function() {
  d <- data.frame(id = 1:20)
  d$lst <- as.list(rep(c("a", "b"), 10))
  d$mat <- matrix(1:40, nrow = 20)
  d
}

test_that("a matrix column is refused rather than silently mis-grouped", {
  d <- odd_df()
  # unique() on a matrix column works element-wise, which silently yields one
  # row from a twenty-row frame rather than an error.
  expect_error(draw(d, design_cluster("mat", n_clusters = 1), seed = 1),
               "which is a 2-dimensional matrix column")
  expect_error(draw(d, design_stratified("mat", n = 4), seed = 1),
               "which is a 2-dimensional matrix column")
  expect_error(draw(d, design_multistage("mat", n_clusters = 1, n = 2), seed = 1),
               "which is a 2-dimensional matrix column")
})

test_that("a list column is refused rather than failing inside base R", {
  d <- odd_df()
  # order() on a list column raises "unimplemented type 'list'" from base R.
  expect_error(draw(d, design_stratified("lst", n = 4), seed = 1),
               "which is a list column")
  expect_error(draw(d, design_systematic(interval = 2, order_by = "lst"), seed = 1),
               "which is a list column")
  expect_error(draw(d, design_cluster("lst", n_clusters = 1), seed = 1),
               "which is a list column")
})

test_that("the error names the offending column and the argument", {
  d <- odd_df()
  expect_error(draw(d, design_stratified(c("id", "mat"), n = 4), seed = 1),
               "`strata` names `mat`")
  expect_error(draw(d, design_systematic(interval = 2, order_by = "lst"), seed = 1),
               "`order_by` names `lst`")
})

test_that("inclusion_prob refuses them too, not just draw()", {
  d <- odd_df()
  expect_error(inclusion_prob(d, design_cluster("mat", n_clusters = 1)),
               "2-dimensional matrix column")
  expect_error(inclusion_prob(d, design_stratified("lst", n = 4)),
               "list column")
})

test_that("factors and ordinary vectors are still accepted", {
  d <- data.frame(
    id = 1:20,
    f  = factor(rep(c("a", "b"), 10)),
    ch = rep(c("x", "y"), 10),
    nu = as.numeric(rep(1:2, 10)),
    dt = rep(as.Date("2020-01-01") + 0:1, 10)
  )
  for (col in c("f", "ch", "nu", "dt")) {
    expect_equal(nrow(draw(d, design_stratified(col, n = 4), seed = 1)), 4,
                 info = col)
    expect_equal(nrow(draw(d, design_cluster(col, n_clusters = 1), seed = 1)), 10,
                 info = col)
  }
  expect_equal(nrow(draw(d, design_systematic(interval = 2, order_by = "nu"),
                         seed = 1)), 10)
})

# A blank field read from a file arrives as "". It is a key value like any
# other, so a frame whose blank is relabelled "_", which sorts to the same
# place, must give the same probabilities and the same draw.
blank_frame <- function(blank = "") {
  data.frame(
    id = 1:30,
    g = rep(c(blank, "a", "b"), each = 10),
    v = c(seq(-20, 20, length.out = 10), seq(-1, 1, length.out = 10),
          seq(-1, 1, length.out = 10)),
    w = c(1000, rep(1:29)),
    y = 1:30
  )
}

test_that("an empty-string stratum is a stratum like any other", {
  d <- blank_frame()
  ref <- blank_frame("_")
  for (des in list(
    design_stratified("g", n = 9),
    design_stratified("g", n = 9, allocation = "neyman", allocation_by = "v",
                      min_per_stratum = 1),
    design_certainty("w", threshold = 500,
                     rest = design_stratified("g", n = 6, min_per_stratum = 2))
  )) {
    expect_equal(inclusion_prob(d, des), inclusion_prob(ref, des))
    expect_equal(joint_prob(d, des, rows = c(2, 3, 12)),
                 joint_prob(ref, des, rows = c(2, 3, 12)))
    s <- draw(d, des, seed = 1, weights = TRUE)
    expect_identical(s$id, draw(ref, des, seed = 1, weights = TRUE)$id)
    expect_false(anyNA(s$.prob))
    expect_equal(ht_total(s, "y")$se, ht_total(draw(ref, des, seed = 1,
                                                    weights = TRUE), "y")$se)
  }
  # Neyman reads the blank stratum's spread, the largest of the three.
  ney <- design_stratified("g", n = 9, allocation = "neyman",
                           allocation_by = "v", min_per_stratum = 1)
  p <- inclusion_prob(d, ney)
  expect_gt(p[1], p[11])
  expect_equal(sum(p), 9)
})

test_that("an empty-string cluster is a cluster like any other", {
  d <- data.frame(id = 1:40, cl = rep(c("", "a", "b", "c"), each = 10),
                  y = 41:80)
  ref <- d
  ref$cl[ref$cl == ""] <- "_"
  des <- design_multistage("cl", n_clusters = 2, n = 6)
  expect_equal(inclusion_prob(d, des), rep(0.15, 40))
  expect_equal(joint_prob(d, des, rows = c(1, 2, 11)),
               joint_prob(ref, des, rows = c(1, 2, 11)))
  s <- draw(d, des, seed = 1, weights = TRUE)
  expect_false(anyNA(s$.prob))
  expect_true(is.finite(ht_total(s, "y")$total))

  cs <- draw(d, design_cluster("cl", n_clusters = 4), seed = 1,
             weights = TRUE)
  sm <- sample_summary(cs)
  expect_equal(sm$by_group$group, c("(blank)", "a", "b", "c"))
  expect_equal(sm$by_group$in_frame, rep(10L, 4))

  cd <- design_cluster("cl", n_clusters = 2)
  expect_equal(joint_prob(d, cd, rows = c(1, 2, 11)),
               joint_prob(ref, cd, rows = c(1, 2, 11)))
  for (seed in 1:3) {
    s <- draw(d, cd, seed = seed, weights = TRUE)
    r <- ht_total(s, "y")
    expect_equal(r$method, "analytic")
    expect_equal(r$se, ht_total(draw(ref, cd, seed = seed, weights = TRUE),
                                "y")$se)
  }
})

test_that("a carriage return inside a key value cannot merge two strata", {
  # Joined with a carriage return, ("x\ry", "z") and ("x", "y\rz") read the
  # same. Keyed on codes, they stay two strata of ten.
  d <- data.frame(id = 1:20, a = rep(c("x\ry", "x"), each = 10),
                  b = rep(c("z", "y\rz"), each = 10))
  des <- design_stratified(c("a", "b"), n = 4, allocation = "equal")
  expect_equal(inclusion_prob(d, des), rep(0.2, 20))
  s <- draw(d, des, seed = 1)
  expect_equal(as.vector(table(s$a)), c(2L, 2L))
})

test_that("blank values are shown as (blank) and columns joined with ' . '", {
  d <- data.frame(a = c("", "x", "x"), b = c("p", "", "q"))
  expect_equal(drawn:::group_label(levels(drawn:::group_key(d, c("a", "b")))),
               c("(blank) . p", "x . (blank)", "x . q"))
})

test_that("keys follow factor levels and byte order, not the locale", {
  d <- data.frame(ch = c("b", "a", "B", "_"),
                  f = factor(c("lo", "hi", "lo", "mid"),
                             levels = c("lo", "mid", "hi")),
                  n = c(10, 2, 1, 2))
  expect_equal(levels(drawn:::group_key(d, "ch")), c("B", "_", "a", "b"))
  expect_equal(levels(drawn:::group_key(d, "f")), c("lo", "mid", "hi"))
  expect_equal(levels(drawn:::group_key(d, "n")), c("1", "2", "10"))
  with_other_collation(
    expect_equal(levels(drawn:::group_key(d, "ch")), c("B", "_", "a", "b"))
  )
})

test_that("many nested strata are keyed in time linear in the rows", {
  skip_on_cran()
  # 2,000 counties by 20,000 tracts is 40 million combinations, of which
  # 20,000 occur. Crossing every level ran out of memory.
  tract <- rep(seq_len(20000), each = 2)
  d <- data.frame(county = sprintf("C%04d", tract %% 2000),
                  tract = sprintf("T%05d", tract))
  des <- design_stratified(c("county", "tract"), n = 20000,
                           min_per_stratum = 1)
  took <- system.time({
    key <- drawn:::group_key(d, c("county", "tract"))
    p <- inclusion_prob(d, des)
  })[["elapsed"]]
  expect_equal(nlevels(key), 20000L)
  expect_equal(p, rep(0.5, 40000))
  expect_lt(took, 20)
})
