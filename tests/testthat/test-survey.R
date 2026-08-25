svy_pop <- function(N = 400) {
  set.seed(1212)
  d <- data.frame(
    id    = seq_len(N),
    site  = rep(c("a", "b", "c", "d"), times = c(N / 2, N / 4, N * 0.15,
                                                 N * 0.1)),
    cl    = rep(paste0("c", seq_len(N / 10)), each = 10),
    spend = round(stats::runif(N, 10, 500))
  )
  # A size measure distinct from the variable being totalled. Selecting
  # proportional to `spend` and then totalling `spend` makes y/pi constant, so
  # every standard error collapses to zero and the comparison proves nothing.
  d$size <- round(d$spend * stats::runif(N, 0.3, 3) + 20)
  d
}

test_that("survey and this package agree on a stratified total", {
  skip_if_not_installed("survey")
  d <- svy_pop()
  s <- draw(d, design_stratified("site", n = 80), seed = 1, weights = TRUE)
  des <- as_svydesign(s)

  ours <- ht_total(s, "spend")
  theirs <- survey::svytotal(~spend, des)
  expect_equal(as.numeric(theirs), ours$total)
  expect_equal(as.numeric(survey::SE(theirs)), ours$se, tolerance = 1e-6)
})

test_that("survey and this package agree on a simple random total", {
  skip_if_not_installed("survey")
  d <- svy_pop()
  s <- draw(d, design_simple(n = 100), seed = 2, weights = TRUE)
  ours <- ht_total(s, "spend")
  theirs <- survey::svytotal(~spend, as_svydesign(s))
  expect_equal(as.numeric(theirs), ours$total)
  expect_equal(as.numeric(survey::SE(theirs)), ours$se, tolerance = 1e-6)
})

test_that("a cluster design maps to cluster ids, and the means agree", {
  skip_if_not_installed("survey")
  d <- svy_pop()
  s <- draw(d, design_cluster("cl", n_clusters = 8), seed = 3, weights = TRUE)
  des <- as_svydesign(s)
  expect_equal(as.numeric(survey::svymean(~spend, des)),
               ht_mean(s, "spend", variance = "none")$mean)
})

test_that("the finite population correction is carried across", {
  skip_if_not_installed("survey")
  d <- svy_pop()
  s <- draw(d, design_simple(n = 200), seed = 4, weights = TRUE)
  with_fpc <- survey::SE(survey::svytotal(~spend, as_svydesign(s)))
  no_fpc <- survey::SE(survey::svytotal(
    ~spend, survey::svydesign(ids = ~1, weights = ~.weight,
                              data = as.data.frame(s))))
  # Half the frame drawn: the correction should cut the standard error clearly.
  expect_lt(as.numeric(with_fpc), as.numeric(no_fpc))
})

test_that("extra arguments reach svydesign", {
  skip_if_not_installed("survey")
  d <- svy_pop()
  s <- draw(d, design_stratified("site", n = 80), seed = 5, weights = TRUE)
  des <- as_svydesign(s, nest = TRUE)
  expect_s3_class(des, "survey.design")
})

test_that("a reservoir sample is simple random sampling, and maps as one", {
  skip_if_not_installed("survey")
  d <- svy_pop()
  r <- draw(d, design_reservoir(n = 100), seed = 6, weights = TRUE)
  ours <- ht_total(r, "spend")
  theirs <- survey::svytotal(~spend, as_svydesign(r))
  expect_equal(as.numeric(theirs), ours$total)
  expect_equal(as.numeric(survey::SE(theirs)), ours$se, tolerance = 1e-6)
})

test_that("a bootstrap design is refused", {
  skip_if_not_installed("survey")
  d <- svy_pop(100)
  # draw() will not weight a bootstrap at all, so reach the guard directly.
  expect_error(draw(d, design_bootstrap(n_replicates = 3), seed = 1,
                    weights = TRUE),
               "not a probability sample")

  fake <- d[1:10, ]
  fake$.weight <- 1
  attr(fake, "drawn_design") <- design_bootstrap(n_replicates = 3)
  attr(fake, "drawn_population") <- d
  expect_error(as_svydesign(fake), "no survey equivalent")
})

test_that("as_svydesign needs a sample carrying its design", {
  skip_if_not_installed("survey")
  d <- svy_pop(100)
  expect_error(as_svydesign(draw(d, design_simple(n = 10), seed = 1)),
               "weights = TRUE")
  expect_error(as_svydesign(d), "weights = TRUE")
})

test_that("survey reproduces this package's standard error, design by design", {
  skip_if_not_installed("survey")
  d <- svy_pop()
  designs <- list(
    simple     = design_simple(n = 100),
    stratified = design_stratified("site", n = 100),
    cluster    = design_cluster("cl", n_clusters = 8),
    reservoir  = design_reservoir(n = 100),
    pps_sys    = design_weighted("size", n = 60, method = "systematic"),
    pps_pois   = design_weighted("size", n = 60, method = "poisson"),
    certainty  = design_certainty("spend", 400,
                                  design_stratified("site", n = 60)),
    cert_srs   = design_certainty("spend", 400, design_simple(n = 60))
  )
  for (nm in names(designs)) {
    s <- draw(d, designs[[nm]], seed = 1, weights = TRUE)
    ours <- ht_total(s, "spend")
    theirs <- survey::svytotal(~spend, as_svydesign(s))
    expect_equal(as.numeric(theirs), ours$total,
                 tolerance = 1e-8, label = paste(nm, "total"))
    expect_equal(as.numeric(survey::SE(theirs)), ours$se,
                 tolerance = 1e-6, label = paste(nm, "se"))
  }
})

test_that("a temporal design maps to strata over its intervals", {
  skip_if_not_installed("survey")
  d <- data.frame(
    id = 1:480,
    ts = rep(seq(as.POSIXct("2024-01-01", tz = "UTC"), by = "day",
                 length.out = 12), each = 40),
    y = round(stats::runif(480, 5, 200))
  )
  des <- design_temporal("ts", from = "2024-01-01", to = "2024-01-13",
                         interval = 1, per_interval = 6, unit = "days")
  s <- draw(d, des, seed = 1, weights = TRUE)
  ours <- ht_total(s, "y")
  theirs <- survey::svytotal(~y, as_svydesign(s))
  expect_equal(as.numeric(survey::SE(theirs)), ours$se, tolerance = 1e-6)
})

test_that("compositions with no survey equivalent are refused, not fudged", {
  skip_if_not_installed("survey")
  d <- svy_pop()
  for (rest in list(design_cluster("cl", n_clusters = 6),
                    design_multistage("cl", n_clusters = 6, n = 30),
                    design_weighted("spend", n = 40, method = "poisson"))) {
    s <- draw(d, design_certainty("spend", 400, rest), seed = 1,
              weights = TRUE)
    expect_error(as_svydesign(s), "no single survey equivalent")
  }
})

test_that("a spatial design's population is the region, not the frame", {
  skip_if_not_installed("survey")
  skip_if_not_installed("sf")
  set.seed(3)
  pop <- data.frame(lon = stats::runif(60, -10, 10),
                    lat = stats::runif(60, -10, 10),
                    y = round(stats::runif(60, 5, 200)))
  region <- sf::st_sfc(
    sf::st_polygon(list(cbind(c(-10, 0, 0, -10, -10),
                              c(-10, -10, 10, 10, -10)))), crs = 4326)
  s <- draw(pop, design_spatial(c("lon", "lat"), region, n = 10), seed = 1,
            weights = TRUE)
  expect_equal(as.numeric(survey::SE(survey::svytotal(~y, as_svydesign(s)))),
               ht_total(s, "y")$se, tolerance = 1e-6)
})

test_that("a `rest` design that puts a row at probability 1 is not a certainty row", {
  skip_if_not_installed("survey")
  # Classifying by `.prob == 1` filed such rows into the certainty stratum,
  # whose fpc counts only the true certainty rows -- svydesign() then rejected
  # the whole design with "FPC implies >100% sampling in some strata".
  set.seed(8)
  pop <- data.frame(v = c(100, 100, 1:10), y = stats::rnorm(12, 50, 5))
  s <- draw(pop, design_certainty("v", 50, design_simple(n = 10)), seed = 1,
            weights = TRUE)
  expect_s3_class(as_svydesign(s), "survey.design")

  # And with a dominant PPS unit below the threshold
  pop2 <- data.frame(v = c(500, 500, 200, stats::runif(20, 1, 5)))
  pop2$y <- round(pop2$v * 2 + stats::rnorm(23, 0, 3))
  s2 <- draw(pop2, design_certainty("v", 400,
                                    design_weighted("v", n = 10,
                                                    method = "systematic")),
             seed = 1, weights = TRUE)
  expect_s3_class(as_svydesign(s2), "survey.design")
})
