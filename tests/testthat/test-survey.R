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
  attr(fake, "drawn_rows") <- 1:10
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

test_that("certainty over clusters or Poisson sampling agrees with survey", {
  skip_if_not_installed("survey")
  # These were refused as having no single survey design. The certainty rows
  # are one more unit beside the clusters, and Poisson sampling covers the
  # whole sample with the certainty rows at probability 1.
  d <- svy_pop()
  d$g <- ifelse(d$id %% 3 == 0, "u", "v")
  rests <- list(
    cluster = design_cluster("cl", n_clusters = 6),
    multistage = design_multistage("cl", n_clusters = 6, n = 30),
    one_row = design_multistage("cl", n_clusters = 6, n = 6),
    balanced = design_cluster("cl", n_clusters = 6, balanced = TRUE),
    poisson = design_weighted("spend", n = 40, method = "poisson"),
    poisson_strata = design_weighted("size", n = 40, method = "poisson",
                                     strata = "site")
  )
  for (nm in names(rests)) {
    s <- draw(d, design_certainty("spend", 400, rests[[nm]]), seed = 1,
              weights = TRUE)
    expect_true(any(s$.prob == 1), label = nm)
    sv <- as_svydesign(s)
    ours <- ht_total(s, "spend")
    theirs <- survey::svytotal(~spend, sv)
    expect_equal(as.numeric(theirs), ours$total, tolerance = 1e-8,
                 label = nm)
    expect_equal(as.numeric(survey::SE(theirs)), ours$se, tolerance = 1e-8,
                 label = nm)
    expect_equal(survey::degf(sv), ours$df, label = nm)
    expect_equal(as.numeric(survey::SE(survey::svymean(~spend, sv))),
                 ht_mean(s, "spend")$se, tolerance = 1e-8, label = nm)
    by_theirs <- suppressWarnings(
      survey::svyby(~spend, ~g, sv, survey::svytotal))
    expect_equal(ht_total(s, "spend", by = "g")$se,
                 as.numeric(survey::SE(by_theirs)), tolerance = 1e-8,
                 label = nm)
  }
})

test_that("a certainty sample whose Poisson `rest` drew nothing is refused", {
  skip_if_not_installed("survey")
  d <- data.frame(v = c(100, 200, rep(1, 30)), y = c(10, 20, 1:30))
  des <- design_certainty("v", 50, design_weighted("v", n = 1,
                                                   method = "poisson"))
  for (i in 1:50) {
    s <- draw(d, des, seed = i, weights = TRUE)
    if (nrow(s) == 2L) break
  }
  expect_equal(nrow(s), 2L)
  expect_error(as_svydesign(s), paste0("`sample` holds none of the 30 rows ",
                                       "with `v` below 50"), fixed = TRUE)
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

test_that("a multistage design is handed over as two stages, and agrees exactly", {
  skip_if_not_installed("survey")
  set.seed(3)
  sizes <- rep(c(8, 12, 16), length.out = 20)
  d <- data.frame(cl = rep(paste0("k", 1:20), times = sizes))
  d$y <- stats::rnorm(nrow(d), rep(stats::rnorm(20, 50, 20), times = sizes), 8)
  s <- draw(d, design_multistage("cl", n_clusters = 5, n = 20), seed = 2,
            weights = TRUE)
  ours <- ht_total(s, "y")
  theirs <- survey::svytotal(~y, as_svydesign(s))
  # Stated as one stage, `survey` saw only the first and came out ~5% low.
  expect_equal(as.numeric(survey::SE(theirs)), ours$se, tolerance = 1e-8)
})

test_that("a sample weighted conditional on its clusters agrees with survey", {
  skip_if_not_installed("survey")
  sizes <- rep(c(5, 9, 14, 7, 11), 4)
  d <- data.frame(id = seq_len(sum(sizes)),
                  cl = rep(sprintf("k%02d", 1:20), times = sizes))
  d$y <- rep((1:20 * 17) %% 23 * 3, times = sizes) + d$id %% 9
  d$g <- ifelse(d$id %% 3 == 0, "u", "v")
  designs <- list(
    design_cluster("cl", n_clusters = 5, balanced = TRUE),
    design_multistage("cl", n_clusters = 5, n = 23),
    design_multistage("cl", n_clusters = 5, n = 25,
                      allocation = "proportional", min_per_cluster = 2)
  )
  for (des in designs) {
    s <- draw(d, des, seed = 4, weights = TRUE)
    sv <- as_svydesign(s)
    ours <- ht_total(s, "y")
    theirs <- survey::svytotal(~y, sv)
    expect_equal(as.numeric(theirs), ours$total, tolerance = 1e-8)
    expect_equal(as.numeric(survey::SE(theirs)), ours$se, tolerance = 1e-8)
    expect_equal(survey::degf(sv), ours$df)
    m_ours <- ht_mean(s, "y")
    m_theirs <- survey::svymean(~y, sv)
    expect_equal(as.numeric(survey::SE(m_theirs)), m_ours$se,
                 tolerance = 1e-8)
    by_ours <- ht_total(s, "y", by = "g")
    by_theirs <- survey::svyby(~y, ~g, sv, survey::svytotal)
    expect_equal(by_ours$se, as.numeric(survey::SE(by_theirs)),
                 tolerance = 1e-8)
  }
})

test_that("a capped reservoir's population is the rows the stream reached", {
  skip_if_not_installed("survey")
  set.seed(1)
  d <- data.frame(y = stats::rnorm(300, 10, 3))
  s <- suppressWarnings(draw(d, design_reservoir(n = 30, max_items = 100),
                             seed = 1, weights = TRUE))
  expect_equal(as.numeric(survey::SE(survey::svytotal(~y, as_svydesign(s)))),
               ht_total(s, "y")$se, tolerance = 1e-8)
})

test_that("systematic PPS maps to survey's Brewer approximation, and agrees closely", {
  skip_if_not_installed("survey")
  set.seed(5)
  pop <- data.frame(x = stats::rgamma(300, 2, 0.1))
  pop$y <- 3 * pop$x + stats::rnorm(300, 0, 15)
  for (des in list(
    design_weighted("x", n = 60, method = "systematic"),
    design_certainty("x", 45, design_weighted("x", n = 40,
                                              method = "systematic"))
  )) {
    s <- draw(pop, des, seed = 1, weights = TRUE)
    theirs <- survey::svytotal(~y, as_svydesign(s))
    ours <- ht_total(s, "y")
    expect_equal(as.numeric(theirs), ours$total, tolerance = 1e-8)
    # Deville here, Brewer there: close relatives, not twins
    expect_equal(as.numeric(survey::SE(theirs)), ours$se, tolerance = 2e-3)
  }
})

test_that("deff() works whether or not survey's has masked it", {
  skip_if_not_installed("survey")
  d <- svy_pop()
  r <- ht_total(draw(d, design_stratified("site", n = 60), seed = 1,
                     weights = TRUE), "spend")
  expect_equal(survey::deff(r), drawn::deff(r))
  m <- ht_mean(draw(d, design_simple(n = 60), seed = 1, weights = TRUE),
               "spend")
  expect_equal(survey::deff(m), drawn::deff(m))
})

test_that("a reordered certainty sample hands survey the right strata", {
  skip_if_not_installed("survey")
  d <- svy_pop()
  s <- draw(d, design_certainty("size", 900, design_simple(n = 40)), seed = 1,
            weights = TRUE)
  want <- ht_total(s, "spend")$se
  r <- s[order(s$spend), ]
  se <- as.numeric(survey::SE(survey::svytotal(~spend, as_svydesign(r))))
  expect_equal(se, want)
})

test_that("survey gets its own copy of the sample", {
  skip_if_not_installed("survey")
  skip_if_not_installed("data.table")
  d <- svy_pop()
  s <- draw(d, design_stratified("site", n = 60), seed = 1, weights = TRUE)
  des <- as_svydesign(s)
  before <- as.numeric(survey::SE(survey::svytotal(~spend, des)))
  data.table::setorderv(s, "spend")
  expect_equal(as.numeric(survey::SE(survey::svytotal(~spend, des))), before)
})

test_that("as_svydesign says what dropped a sample's design", {
  skip_if_not_installed("survey")
  d <- svy_pop()
  s <- draw(d, design_stratified("site", n = 60), seed = 1, weights = TRUE)
  merged <- merge(s, data.frame(id = s$id, resp = 1), by = "id")
  expect_error(as_svydesign(merged), "does not carry its design")
  expect_error(as_svydesign(merged), "merge()", fixed = TRUE)
  expect_error(as_svydesign(s[s$site == "a", ]), "rows that were drawn")
})

test_that("survey follows `.weight`, while ht_total keeps to `.prob`", {
  skip_if_not_installed("survey")
  d <- svy_pop()
  s <- draw(d, design_stratified("site", n = 60), seed = 1, weights = TRUE)
  w <- s
  w$.weight <- w$.weight * 1.1
  theirs <- survey::svytotal(~spend, as_svydesign(w))
  expect_equal(as.numeric(theirs), 1.1 * ht_total(s, "spend")$total)
  expect_warning(ours <- ht_total(w, "spend"), class = "drawn_weight_edited")
  expect_equal(ours$total, ht_total(s, "spend")$total)

  p <- s
  p$.prob <- p$.prob / 2
  expect_error(as_svydesign(p), "`.prob` was changed after draw()",
               fixed = TRUE)
})

test_that("survey gets the strata and clusters drawn, whatever the sample says", {
  skip_if_not_installed("survey")
  d <- svy_pop()
  designs <- list(
    site = design_stratified("site", n = 60),
    cl = design_cluster("cl", n_clusters = 8),
    cl = design_multistage("cl", n_clusters = 8, n = 40)
  )
  for (i in seq_along(designs)) {
    key <- names(designs)[i]
    s <- draw(d, designs[[i]], seed = 1, weights = TRUE)
    want <- ht_total(s, "spend")$se
    renamed <- s
    names(renamed)[names(renamed) == key] <- "renamed"
    recoded <- s
    recoded[[key]] <- rev(recoded[[key]])
    for (e in list(renamed, recoded)) {
      se <- as.numeric(survey::SE(survey::svytotal(~spend, as_svydesign(e))))
      expect_equal(se, want, tolerance = 1e-8, label = key)
    }
  }
})

# ---- units taken with certainty ----------------------------------------------

dominant_svy <- function() {
  # Four units heavy enough that PPS of 8 takes them with certainty
  pop <- data.frame(id = 1:40, v = c(1000, 900, 300, 280, 1 + (1:36 %% 5)),
                    a = cos(1:40), b = sin(1.7 * (1:40)))
  pop$y <- 2 * pop$v + 2 * cos(3 * (1:40))
  pop$flag <- pop$y > stats::median(pop$y)
  pop
}

test_that("units PPS takes with certainty are a stratum of their own", {
  skip_if_not_installed("survey")
  # Left in the Brewer stratum, they put survey's standard error at 35 and 25
  # times this package's here.
  pop <- dominant_svy()
  for (des in list(
    design_weighted("v", n = 8, method = "systematic"),
    design_certainty("v", 500, design_weighted("v", n = 6,
                                               method = "systematic"))
  )) {
    s <- draw(pop, des, seed = 1, weights = TRUE)
    expect_equal(sum(s$.prob == 1), 4L)
    sv <- as_svydesign(s)
    ours <- ht_total(s, "y")
    theirs <- survey::svytotal(~y, sv)
    expect_equal(as.numeric(theirs), ours$total)
    expect_equal(as.numeric(survey::SE(theirs)), ours$se, tolerance = 5e-3)
    expect_equal(as.numeric(survey::SE(survey::svymean(~y, sv))),
                 ht_mean(s, "y")$se, tolerance = 5e-3)
    expect_equal(survey::degf(sv), ours$df)
    expect_equal(ours$df, 3)
  }
})

test_that("Poisson sampling hands its certain units over the same way", {
  skip_if_not_installed("survey")
  pop <- dominant_svy()
  s <- draw(pop, design_weighted("v", n = 8, method = "poisson"), seed = 2,
            weights = TRUE)
  expect_equal(sum(s$.prob == 1), 4L)
  sv <- as_svydesign(s)
  ours <- ht_total(s, "y")
  expect_equal(as.numeric(survey::SE(survey::svytotal(~y, sv))), ours$se,
               tolerance = 1e-8)
  expect_equal(survey::degf(sv), ours$df)
  expect_equal(ours$df, nrow(s) - 5)
})

test_that("spread with size leaves its certain units out of Brewer", {
  skip_if_not_installed("survey")
  pop <- dominant_svy()
  s <- draw(pop, design_spread(c("a", "b"), n = 8, size = "v"), seed = 1,
            weights = TRUE)
  expect_equal(sum(s$.prob == 1), 4L)
  sv <- as_svydesign(s)
  free <- as.data.frame(s)[s$.prob < 1, ]
  brewer <- survey::svydesign(ids = ~1, fpc = ~.prob, weights = ~.weight,
                              data = free, pps = "brewer")
  expect_equal(as.numeric(survey::SE(survey::svytotal(~y, sv))),
               as.numeric(survey::SE(survey::svytotal(~y, brewer))))
  expect_equal(survey::degf(sv), ht_total(s, "y")$df)
})

test_that("certainty designs agree with survey on degrees of freedom", {
  skip_if_not_installed("survey")
  # The certainty rows were each a unit of survey's, so survey::degf() counted
  # them and survey's t and logit intervals came out narrower.
  d <- svy_pop()
  d$flag <- d$spend > 250
  d$tiny <- c(rep("x", 6), rep("y", nrow(d) - 6))
  designs <- list(
    cert_srs = design_certainty("size", 900, design_simple(n = 40)),
    cert_str = design_certainty("size", 900,
                                design_stratified("site", n = 40)),
    take_all = design_stratified("tiny", n = 30, allocation = "equal")
  )
  for (nm in names(designs)) {
    s <- draw(d, designs[[nm]], seed = 1, weights = TRUE)
    expect_true(any(s$.prob == 1), label = nm)
    sv <- as_svydesign(s)
    ours <- ht_total(s, "spend")
    expect_equal(survey::degf(sv), ours$df, label = nm)
    expect_equal(as.numeric(survey::SE(survey::svytotal(~spend, sv))),
                 ours$se, tolerance = 1e-8, label = nm)
    expect_equal(unname(stats::confint(survey::svytotal(~spend, sv),
                                       df = survey::degf(sv))[1, ]),
                 ours$ci, tolerance = 1e-8, label = nm)
    ci <- attr(survey::svyciprop(~flag, sv, method = "xlogit"), "ci")
    expect_equal(unname(ci), ht_mean(s, "flag")$ci, tolerance = 1e-8,
                 label = nm)
  }

  # A systematic PPS `rest` agrees on the degrees of freedom; the standard
  # errors differ by Brewer against Deville.
  s <- draw(d, design_certainty("size", 900,
                                design_weighted("size", n = 40,
                                                method = "systematic")),
            seed = 1, weights = TRUE)
  expect_equal(survey::degf(as_svydesign(s)), ht_total(s, "spend")$df)
})

test_that("a census is handed to survey with a variance of 0", {
  skip_if_not_installed("survey")
  # The PPS mappings stated a census as probabilities of 1, and svydesign()
  # stopped with "subscript out of bounds".
  pop <- dominant_svy()
  for (des in list(design_weighted("v", n = 40, method = "systematic"),
                   design_weighted("v", n = 40, method = "poisson"),
                   design_spread(c("a", "b"), n = 40, size = "v"),
                   design_certainty("v", 500,
                                    design_weighted("v", n = 38,
                                                    method = "systematic")),
                   design_stratified("flag", n = 40))) {
    s <- draw(pop, des, seed = 1, weights = TRUE)
    expect_true(all(s$.prob == 1))
    sv <- as_svydesign(s)
    theirs <- survey::svytotal(~y, sv)
    expect_equal(as.numeric(theirs), sum(pop$y))
    expect_equal(as.numeric(survey::SE(theirs)), 0)
    # The one place the degrees of freedom differ, with nothing to scale
    expect_equal(survey::degf(sv), 39)
    expect_equal(ht_total(s, "y")$df, Inf)
  }
})

test_that("an empty sample is refused before it reaches survey", {
  skip_if_not_installed("survey")
  pop <- dominant_svy()
  s <- draw(pop, design_spread(c("a", "b"), n = 0), seed = 1, weights = TRUE)
  expect_error(as_svydesign(s), "`sample` has no rows", fixed = TRUE)
})

test_that("strata from several columns reach survey without carriage returns", {
  skip_if_not_installed("survey")
  d <- svy_pop()
  d$flag <- d$spend > 250
  s <- draw(d, design_stratified(c("site", "flag"), n = 60,
                                 min_per_stratum = 2), seed = 7, weights = TRUE)
  st <- as.character(as_svydesign(s)$strata[, 1])
  expect_false(any(grepl("\r", st, fixed = TRUE)))
  expect_setequal(unique(st), paste(s$site, s$flag, sep = " . "))
  expect_equal(as.numeric(survey::SE(survey::svytotal(~spend,
                                                      as_svydesign(s)))),
               ht_total(s, "spend")$se, tolerance = 1e-6)

  # Labels that would read alike once joined are handed over as codes, so
  # survey keeps the two strata apart.
  e <- data.frame(id = 1:40, a = rep(c("x . y", "x"), each = 20),
                  b = rep(c("z", "y . z"), each = 20),
                  spend = 1:40 + 0.5 * (1:40 %% 3))
  s2 <- draw(e, design_stratified(c("a", "b"), n = 8), seed = 1,
             weights = TRUE)
  st2 <- as.character(as_svydesign(s2)$strata[, 1])
  expect_equal(length(unique(st2)), 2L)
  expect_equal(as.numeric(survey::SE(survey::svytotal(~spend,
                                                      as_svydesign(s2)))),
               ht_total(s2, "spend")$se, tolerance = 1e-6)
})

test_that("ci_method = \"beta\" is survey's beta interval", {
  skip_if_not_installed("survey")
  # Clusters that differ in their share of flagged rows, so the effective
  # sample size, p(1 - p) / v, is well below the 60 rows and is not capped.
  d <- data.frame(id = 1:200, cl = rep(1:20, each = 10))
  d$flag <- (d$id %% 10) < rep(c(1, 9, 2, 8, 0, 10, 3, 7, 5, 4), 2)[d$cl]
  s <- draw(d, design_cluster("cl", n_clusters = 6), seed = 3, weights = TRUE)
  sv <- as_svydesign(s)
  for (level in c(0.95, 0.90)) {
    ours <- ht_mean(s, "flag", ci_method = "beta", level = level)
    expect_true(ours$mean * (1 - ours$mean) / ours$variance < ours$n)
    theirs <- survey::svyciprop(~flag, sv, method = "beta", level = level)
    expect_equal(ours$ci, unname(as.vector(attr(theirs, "ci"))),
                 tolerance = 1e-10)
  }
  expect_equal(ht_mean(s, "flag", ci_method = "beta")$ci_scale, "beta")
})

test_that("deff agrees with survey's where the design reaches part of the frame", {
  skip_if_not_installed("survey")
  check <- function(s, label) {
    sv <- as_svydesign(s)
    tot <- ht_total(s, "y")
    m <- ht_mean(s, "y")
    expect_equal(tot$deff,
                 unname(survey::deff(survey::svytotal(~y, sv, deff = TRUE))),
                 tolerance = 1e-8, label = label)
    expect_equal(m$deff,
                 unname(survey::deff(survey::svymean(~y, sv, deff = TRUE))),
                 tolerance = 1e-8, label = label)
    expect_true(abs(tot$deff - 1) < 0.35, label = label)
  }

  # The first 100 of 400 rows
  d <- data.frame(y = 50 + (1:400 * 37) %% 23)
  expect_warning(
    s <- draw(d, design_reservoir(n = 20, max_items = 100), seed = 1,
              weights = TRUE), "max_items")
  check(s, "reservoir")

  # Two days of a four-day frame
  dt <- data.frame(ts = seq(as.POSIXct("2020-01-01", tz = "UTC"), by = "hour",
                            length.out = 96),
                   y = 50 + (1:96 * 37) %% 23)
  st <- draw(dt, design_temporal("ts", "2020-01-02", "2020-01-04",
                                 interval = 12, per_interval = 4),
             seed = 2, weights = TRUE)
  check(st, "temporal")

  # Points in half the frame's extent
  skip_if_not_installed("sf")
  sp <- data.frame(lon = -10 + (1:80 * 7) %% 20, lat = -10 + (1:80 * 3) %% 20,
                   y = 50 + (1:80 * 37) %% 23)
  region <- sf::st_sfc(
    sf::st_polygon(list(cbind(c(-10.5, 0, 0, -10.5, -10.5),
                              c(-10.5, -10.5, 10, 10, -10.5)))), crs = 4326)
  ss <- draw(sp, design_spatial(c("lon", "lat"), region, n = 12), seed = 1,
             weights = TRUE)
  expect_true(sum(attr(ss, "drawn_pi") == 0) > 0)
  check(ss, "spatial")
})

test_that("a systematic design's fpc counts the rows its walk can reach", {
  skip_if_not_installed("survey")
  # Rows `na_rm` drops from the walk can never be drawn. The same draw from
  # the frame without them gets the same design.
  d <- data.frame(id = 1:60, a = (1:60 %% 7) / 7, y = 10 + (1:60 * 13) %% 17)
  d$a[5:8] <- NA
  s1 <- draw(d, design_systematic(interval = 5, order_by = "a", na_rm = TRUE),
             seed = 4, weights = TRUE)
  s2 <- draw(d[!is.na(d$a), ], design_systematic(interval = 5, order_by = "a"),
             seed = 4, weights = TRUE)
  expect_equal(s1$id, s2$id)
  sv1 <- as_svydesign(s1)
  expect_equal(unique(sv1$fpc$popsize[, 1]), 56)
  expect_equal(survey::SE(survey::svytotal(~y, sv1)),
               survey::SE(survey::svytotal(~y, as_svydesign(s2))))
})

test_that("PPS within strata maps to strata, with Brewer's in each", {
  skip_if_not_installed("survey")
  d <- svy_pop()
  d$size[c(1, 250)] <- c(9000, 7000)
  for (m in c("systematic", "poisson")) {
    s <- draw(d, design_weighted("size", n = 60, method = m, strata = "site"),
              seed = 1, weights = TRUE)
    des <- as_svydesign(s)
    ours <- ht_total(s, "spend")
    theirs <- survey::svytotal(~spend, des)
    expect_equal(as.numeric(theirs), ours$total, tolerance = 1e-8)
    expect_equal(survey::degf(des), ours$df)
    if (m == "poisson") {
      expect_equal(as.numeric(survey::SE(theirs)), ours$se, tolerance = 1e-8)
    } else {
      # Deville here, Brewer there, stratum by stratum.
      expect_equal(as.numeric(survey::SE(theirs)), ours$se, tolerance = 5e-3)
    }
  }
})

test_that("clusters selected by size go to survey in two stages", {
  skip_if_not_installed("survey")
  sizes <- c(120, 90, rep(c(8, 12, 16, 20), length.out = 28))
  d <- data.frame(cl = rep(sprintf("k%02d", 1:30), times = sizes))
  d$rows <- rep(sizes, times = sizes)
  d$y <- rep((1:30 * 37) %% 41 + 20, times = sizes) +
    seq_len(nrow(d)) %% 9 * 2
  # Ten clusters take k01 and k02 with certainty; without them, none is.
  cases <- list(list(d, 10), list(d[!d$cl %in% c("k01", "k02"), ], 6))
  for (case in cases) {
    n_clusters <- case[[2]]
    s <- draw(case[[1]], design_multistage("cl", n_clusters = n_clusters,
                                           n = 4 * n_clusters, size = "rows"),
              seed = 2, weights = TRUE)
    des <- as_svydesign(s)
    ours <- ht_total(s, "y")
    theirs <- survey::svytotal(~y, des)
    expect_equal(as.numeric(theirs), ours$total, tolerance = 1e-8)
    expect_equal(survey::degf(des), ours$df)
    # Brewer and Deville at the first stage, the same terms within clusters.
    expect_equal(as.numeric(survey::SE(theirs)), ours$se, tolerance = 5e-3)
  }
  expect_identical(length(unique(as_svydesign(s)$strata[, 1])), 1L)
})

# ---- column names --------------------------------------------------------------

test_that("cluster and strata columns of any name reach survey", {
  skip_if_not_installed("survey")
  # The cluster column's name went into the `ids` formula as it stood: a
  # space stopped svydesign(), and `cl-1` beside a column `cl` was read as
  # `cl - 1`, so survey clustered on `cl`.
  d <- data.frame(`cl-1` = rep(paste0("B", 1:8), each = 30),
                  check.names = FALSE)
  d$cl <- rep(c("north", "south"), each = 120)
  d$y <- rep(c(140, 60, 95, 130, 70, 110, 85, 120), each = 30) +
    (1:240 * 7) %% 11
  s <- draw(d, design_cluster("cl-1", n_clusters = 4), seed = 2,
            weights = TRUE)
  sv <- as_svydesign(s)
  expect_equal(length(unique(sv$cluster[, 1])), 4L)
  expect_equal(as.numeric(survey::SE(survey::svytotal(~y, sv))),
               ht_total(s, "y")$se, tolerance = 1e-8)
  expect_equal(survey::degf(sv), 3)

  names(d)[1] <- "cluster id"
  d$`site name` <- rep(c("a", "b"), 120)
  for (des in list(design_cluster("cluster id", n_clusters = 4),
                   design_multistage("cluster id", n_clusters = 4, n = 20),
                   design_stratified("site name", n = 20))) {
    s <- draw(d, des, seed = 2, weights = TRUE)
    expect_equal(as.numeric(survey::SE(survey::svytotal(~y,
                                                        as_svydesign(s)))),
                 ht_total(s, "y")$se, tolerance = 1e-8)
  }
})

test_that("a sample with a column as_svydesign() writes is refused", {
  skip_if_not_installed("survey")
  # A user column `.stratum` was overwritten with the design's strata, so
  # survey::svyby(~y, ~.stratum) estimated by site.
  d <- svy_pop()
  d$.stratum <- ifelse(d$spend > 250, "hi", "lo")
  s <- draw(d, design_stratified("site", n = 40), seed = 7, weights = TRUE)
  expect_error(as_svydesign(s),
               "`sample` already has `.stratum`; as_svydesign() would overwrite it.",
               fixed = TRUE)
  d$.fpc <- 1
  s <- draw(d, design_simple(n = 40), seed = 7, weights = TRUE)
  expect_error(as_svydesign(s), "already has `.stratum` and `.fpc`",
               fixed = TRUE)

  # Renamed in the sample, the design is unchanged.
  names(s)[names(s) %in% c(".stratum", ".fpc")] <- c("band", "f")
  expect_equal(as.numeric(survey::SE(survey::svytotal(~spend,
                                                      as_svydesign(s)))),
               ht_total(s, "spend")$se, tolerance = 1e-8)

  # A design column of a reserved name is in the sample too.
  e <- svy_pop()
  names(e)[names(e) == "cl"] <- ".cluster"
  s <- draw(e, design_cluster(".cluster", n_clusters = 8), seed = 1,
            weights = TRUE)
  expect_error(as_svydesign(s), "already has `.cluster`", fixed = TRUE)
  names(s)[names(s) == ".cluster"] <- "cl"
  expect_equal(as.numeric(survey::SE(survey::svytotal(~spend,
                                                      as_svydesign(s)))),
               ht_total(s, "spend")$se, tolerance = 1e-8)
})

test_that("unnamed arguments to as_svydesign() are refused", {
  skip_if_not_installed("survey")
  # modifyList() dropped them, so as_svydesign(s, ~cl) ignored the ~cl.
  d <- svy_pop()
  s <- draw(d, design_simple(n = 40), seed = 1, weights = TRUE)
  expect_error(as_svydesign(s, ~cl),
               "Argument 1 of `...` has no name", fixed = TRUE)
  expect_error(as_svydesign(s, nest = TRUE, ~cl, TRUE),
               "Arguments 2 and 3 of `...` have no name", fixed = TRUE)
  expect_s3_class(as_svydesign(s, ids = ~cl), "survey.design")
})

test_that("one row from a cluster of several goes to survey as ultimate clusters", {
  skip_if_not_installed("survey")
  # Two stages made survey stop with "has only one PSU at stage 2", where
  # ht_total() falls back to the delete-a-cluster jackknife.
  sizes <- rep(c(5, 9, 14, 7, 11), 4)
  d <- data.frame(id = seq_len(sum(sizes)),
                  cl = rep(sprintf("k%02d", 1:20), times = sizes))
  d$y <- rep((1:20 * 17) %% 23 * 3, times = sizes) + d$id %% 9
  d$g <- ifelse(d$id %% 3 == 0, "u", "v")
  designs <- list(
    design_multistage("cl", n_clusters = 5, n = 5),
    design_multistage("cl", n_clusters = 6, n = 9,
                      allocation = "proportional", min_per_cluster = 1)
  )
  for (des in designs) {
    s <- draw(d, des, seed = 4, weights = TRUE)
    sv <- as_svydesign(s)
    ours <- ht_total(s, "y")
    expect_equal(ours$method, "jackknife")
    expect_null(sv$fpc$popsize)
    theirs <- survey::svytotal(~y, sv)
    expect_equal(as.numeric(survey::SE(theirs)), ours$se, tolerance = 1e-8)
    expect_equal(survey::degf(sv), ours$df)
    expect_equal(as.numeric(survey::SE(survey::svymean(~y, sv))),
                 ht_mean(s, "y")$se, tolerance = 1e-8)
    by_theirs <- survey::svyby(~y, ~g, sv, survey::svytotal)
    expect_equal(ht_total(s, "y", by = "g")$se,
                 as.numeric(survey::SE(by_theirs)), tolerance = 1e-8)
  }
})

test_that("spread and systematic samples go to survey as their simple designs", {
  skip_if_not_installed("survey")
  # survey's standard error is that of the first-order design. For spread it
  # is the larger when y varies smoothly over the coordinates, and the
  # smaller here, where y has no spatial pattern and half the frame is drawn.
  i <- 1:200
  d <- data.frame(x = (i * 0.618034) %% 1, yy = (i * 0.414214) %% 1)
  d$smooth <- 100 + 80 * d$x + 60 * d$yy
  d$flat <- 50 + (i * 37) %% 23
  srs_se <- function(s, y) {
    plain <- survey::svydesign(ids = ~1, fpc = rep(nrow(d), nrow(s)),
                               weights = ~.weight, data = as.data.frame(s))
    as.numeric(survey::SE(survey::svytotal(stats::as.formula(y), plain)))
  }
  s <- draw(d, design_spread(c("x", "yy"), n = 20), seed = 1, weights = TRUE)
  se <- as.numeric(survey::SE(survey::svytotal(~smooth, as_svydesign(s))))
  expect_equal(se, srs_se(s, "~smooth"))
  expect_gt(se, ht_total(s, "smooth")$se)
  s <- draw(d, design_spread(c("x", "yy"), n = 100), seed = 1,
            weights = TRUE)
  se <- as.numeric(survey::SE(survey::svytotal(~flat, as_svydesign(s))))
  expect_equal(se, srs_se(s, "~flat"))
  expect_lt(se, ht_total(s, "flat")$se)

  # Systematic: the simple-random variance against successive differences
  d$trend <- i + 10 * sin(i)
  s <- draw(d, design_systematic(interval = 10), seed = 3, weights = TRUE)
  se <- as.numeric(survey::SE(survey::svytotal(~trend, as_svydesign(s))))
  expect_equal(se, srs_se(s, "~trend"))
  expect_gt(se, 2 * ht_total(s, "trend")$se)
})
