weighted_df <- function() {
  set.seed(99)
  data.frame(id = 1:100, weight = stats::runif(100, 0.1, 10))
}

test_that("draws the requested number of rows", {
  df <- weighted_df()
  expect_equal(nrow(draw(df, design_weighted("weight", n = 10), seed = 123)), 10)
})

test_that("the caller's weights come back untouched", {
  d <- data.frame(id = 1:10, w = 1:10)
  out <- draw(d, design_weighted("w", n = 5), seed = 1)
  expect_equal(out$w, d$w[match(out$id, d$id)])
  expect_type(out$w, "integer")
})

test_that("heavier rows are drawn more often", {
  d <- data.frame(id = 1:4, w = c(1, 1, 1, 100))
  draws <- vapply(1:400, function(i) draw(d, design_weighted("w", n = 1),
                                          seed = i)$id, numeric(1))
  expect_gt(mean(draws == 4), 0.8)
})

test_that("even the lightest row remains selectable", {
  d <- data.frame(id = 1:5, w = c(1, 2, 3, 4, 5))
  drawn <- vapply(1:600, function(i) draw(d, design_weighted("w", n = 1),
                                          seed = i)$id, numeric(1))
  expect_true(1 %in% drawn)
})

test_that("non-positive, non-finite and non-numeric weights are rejected", {
  df <- weighted_df()

  bad <- df
  bad$weight <- -abs(bad$weight)
  expect_error(draw(bad, design_weighted("weight", n = 10), seed = 1),
               "must be positive")

  inf_df <- df
  inf_df$weight[1] <- Inf
  expect_error(draw(inf_df, design_weighted("weight", n = 10), seed = 1),
               "non-finite")

  chr <- df
  chr$weight <- as.character(chr$weight)
  expect_error(draw(chr, design_weighted("weight", n = 10), seed = 1),
               "must be numeric")
})

test_that("NA weights follow the same na_rm contract as other keys", {
  df <- weighted_df()
  df$weight[1] <- NA
  expect_error(draw(df, design_weighted("weight", n = 10), seed = 1),
               "missing weight")
  expect_equal(
    nrow(draw(df, design_weighted("weight", n = 10, na_rm = TRUE), seed = 1)),
    10
  )
})

test_that("over-drawing without replacement is refused", {
  d <- data.frame(id = 1:5, w = 1:5)
  expect_error(draw(d, design_weighted("w", n = 10), seed = 1),
               "cannot exceed the number of rows")
})

test_that("weights draw() refuses are refused on every path, naming the column", {
  # inclusion_prob() accepted logical weights and stopped inside base R on a
  # factor; joint_prob(), draw(weights = TRUE) and plot() followed it.
  pop <- data.frame(id = 1:20, size = 1:20, other = 1:20)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  bad <- list(logical = rep(TRUE, 20), factor = factor(rep(1:4, 5)),
              character = as.character(1:20))
  for (type in names(bad)) {
    pop$size <- bad[[type]]
    want <- paste0("`size` must be numeric, not ", type, ".")
    designs <- list(
      design_weighted("size", n = 5),
      design_weighted("size", n = 5, method = "systematic"),
      design_weighted("size", n = 5, method = "poisson"),
      design_certainty("other", 18, design_weighted("size", n = 5,
                                                    method = "poisson"))
    )
    for (des in designs) {
      info <- paste(type, drawn:::design_label(des))
      expect_error(draw(pop, des, seed = 1), want, fixed = TRUE, info = info)
      expect_error(inclusion_prob(pop, des), want, fixed = TRUE, info = info)
      expect_error(sampling_weight(pop, des), want, fixed = TRUE, info = info)
      expect_error(joint_prob(pop, des, rows = 1:3), want, fixed = TRUE,
                   info = info)
      expect_error(draw(pop, des, seed = 1, weights = TRUE), want,
                   fixed = TRUE, info = info)
      expect_error(plot(des, pop, type = "probability"), want, fixed = TRUE,
                   info = info)
    }
  }
})

test_that("a weight of 0 is refused with the column and what to do", {
  pop <- data.frame(w = c(0, 0, 0, 1:20))
  for (m in c("successive", "systematic", "poisson")) {
    des <- design_weighted("w", n = 5, method = m)
    expect_error(draw(pop, des, seed = 1),
                 "`w` must be positive, but 3 weight(s) are 0 or negative.",
                 fixed = TRUE)
    expect_error(inclusion_prob(pop, des), "Drop those rows from `data`",
                 fixed = TRUE)
  }
  expect_error(inclusion_prob(pop, design_spread("w", n = 5, size = "w")),
               "`w` must be positive", fixed = TRUE)
})

test_that("a design without strata draws the rows it drew before", {
  d <- data.frame(id = 1:30, w = c(1:29, 90))
  for (m in c("systematic", "poisson")) {
    des <- design_weighted("w", n = 6, method = m)
    expect_null(des$strata)
    expect_null(des$allocation)
    pi <- pps_pi(d$w, 6)
    old <- with_seed(4, if (m == "poisson") {
      which(stats::runif(30) < pi)
    } else {
      systematic_pps(pi)
    })
    expect_identical(draw(d, des, seed = 4)$id, old)
  }
  expect_false(any(grepl("strata|allocation",
                         utils::capture.output(print(
                           design_weighted("w", n = 6))))))
})

test_that("strata and allocation print, and other designs are not stratified", {
  des <- design_weighted("w", n = 6, method = "systematic", strata = "g",
                         allocation = "equal")
  out <- utils::capture.output(print(des))
  expect_true(any(grepl("strata +\"g\"", out)))
  expect_true(any(grepl("allocation +\"equal\"", out)))
  expect_error(allocation_table(data.frame(w = 1:4, g = "a"),
                                design_weighted("w", n = 2,
                                                method = "systematic")),
               "design_weighted() with `strata`", fixed = TRUE)
})

test_that("strata of a weighted design are summarised from the frame", {
  d <- data.frame(id = 1:40, g = rep(c("x", "y"), c(30, 10)),
                  w = rep(1:4, 10))
  s <- draw(d, design_weighted("w", n = 8, method = "systematic",
                               strata = "g"), seed = 1, weights = TRUE)
  sm <- sample_summary(s)
  expect_identical(sm$group_col, "g")
  expect_identical(sm$by_group$drawn, c(6L, 2L))
  expect_identical(sm$by_group$in_frame, c(30L, 10L))
  s$g <- "renamed"
  expect_identical(sample_summary(s)$by_group$drawn, c(6L, 2L))
})
