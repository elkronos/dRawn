dom_pop <- function() {
  set.seed(21)
  d <- data.frame(
    id = 1:400,
    site = rep(c("a", "b", "c", "d"), times = c(160, 120, 80, 40)),
    region = rep(c("east", "west"), 200),
    cl = rep(paste0("c", 1:40), each = 10)
  )
  d$y <- round(stats::runif(400, 10, 500))
  d
}

test_that("domain totals add up to the overall total", {
  d <- dom_pop()
  s <- draw(d, design_stratified("site", n = 80), seed = 1, weights = TRUE)
  by_region <- ht_total(s, "y", by = "region")
  expect_s3_class(by_region, "drawn_by")
  expect_equal(sum(by_region$total), ht_total(s, "y")$total)
  expect_equal(sort(by_region$region), c("east", "west"))
  expect_equal(sum(by_region$n), nrow(s))
})

test_that("domain estimates are unbiased across repeated samples", {
  d <- dom_pop()
  truth <- tapply(d$y, d$region, sum)
  des <- design_stratified("site", n = 80)
  est <- vapply(1:300, function(i) {
    r <- ht_total(draw(d, des, seed = i, weights = TRUE), "y", by = "region",
                  variance = "none")
    r$total[order(r$region)]
  }, numeric(2))
  expect_equal(rowMeans(est), as.numeric(truth), tolerance = 0.02)
})

test_that("domain standard errors match survey::svyby, design by design", {
  skip_if_not_installed("survey")
  d <- dom_pop()
  for (des in list(design_simple(n = 80), design_stratified("site", n = 80),
                   design_cluster("cl", n_clusters = 8))) {
    s <- draw(d, des, seed = 2, weights = TRUE)
    svd <- as_svydesign(s)

    ours <- ht_total(s, "y", by = "region")
    theirs <- survey::svyby(~y, ~region, svd, survey::svytotal)
    expect_equal(ours$total, unname(stats::coef(theirs)), tolerance = 1e-8)
    expect_equal(ours$se, unname(survey::SE(theirs)), tolerance = 1e-6)

    ours_m <- ht_mean(s, "y", by = "region")
    theirs_m <- survey::svyby(~y, ~region, svd, survey::svymean)
    expect_equal(ours_m$mean, unname(stats::coef(theirs_m)), tolerance = 1e-8)
    expect_equal(ours_m$se, unname(survey::SE(theirs_m)), tolerance = 1e-6)
  }
})

test_that("the Hajek domain variance tracks the real sampling variance", {
  d <- dom_pop()
  des <- design_simple(n = 60)
  out <- vapply(1:500, function(i) {
    r <- ht_mean(draw(d, des, seed = i, weights = TRUE), "y", by = "region")
    r <- r[r$region == "west", ]
    c(r$mean, r$se^2)
  }, numeric(2))
  ratio <- mean(out[2, ]) / stats::var(out[1, ])
  expect_gt(ratio, 0.8)
  expect_lt(ratio, 1.25)
})

test_that("the HT domain mean divides by the domain's known size", {
  d <- dom_pop()
  s <- draw(d, design_simple(n = 80), seed = 3, weights = TRUE)
  r <- ht_mean(s, "y", by = "site", estimator = "ht")
  N_a <- sum(d$site == "a")
  expect_equal(r$mean[r$site == "a"],
               sum(s$y[s$site == "a"] * s$.weight[s$site == "a"]) / N_a)
})

test_that("by accepts several columns and validates them", {
  d <- dom_pop()
  s <- draw(d, design_simple(n = 80), seed = 4, weights = TRUE)
  two <- ht_total(s, "y", by = c("site", "region"))
  expect_true(all(c("site", "region") %in% names(two)))
  expect_equal(sum(two$total), ht_total(s, "y")$total)
  expect_error(ht_total(s, "y", by = "nope"),
               "`by` names `nope`, which `sample` does not have.", fixed = TRUE)
  expect_error(ht_total(s, "y", by = 1), "`by` must name")
  expect_output(print(two), "by domain")
})

test_that("a recoded stratum column does not change domain variances", {
  d <- dom_pop()
  s <- draw(d, design_stratified("site", n = 80, min_per_stratum = 2), seed = 4,
            weights = TRUE)
  recoded <- s
  recoded$site <- "z"
  for (v in c("analytic", "jackknife")) {
    want <- ht_total(s, "y", by = "region", variance = v)
    got <- ht_total(recoded, "y", by = "region", variance = v)
    expect_equal(got$se, want$se, label = v)
    expect_equal(attr(got, "df"), attr(want, "df"), label = v)
  }
})

test_that("domains come out in key order, not sorted as text", {
  d <- dom_pop()
  d$month <- rep(1:12, length.out = 400)
  d$mfac <- factor(month.abb[d$month], levels = month.abb)
  d$case <- c("a", "B", "b")[(d$id %/% 12) %% 3 + 1]
  s <- draw(d, design_simple(n = 120), seed = 1, weights = TRUE)
  expect_equal(ht_total(s, "y", by = "month")$month, 1:12)
  expect_equal(as.character(ht_total(s, "y", by = "mfac")$mfac), month.abb)
  expect_equal(ht_total(s, "y", by = "case")$case, c("B", "a", "b"))
  with_other_collation(
    expect_equal(ht_total(s, "y", by = "case")$case, c("B", "a", "b"))
  )
  # Several columns: the first varies slowest.
  two <- ht_total(s, "y", by = c("mfac", "case"))
  expect_gt(nrow(two), 12)
  expect_identical(order(as.integer(two$mfac), match(two$case, c("B", "a", "b"))),
                   seq_len(nrow(two)))
})

test_that("domains and summaries of an sf sample carry no geometry", {
  pts <- sf_points()
  plain <- sf::st_drop_geometry(pts)
  des <- design_cluster("cl", n_clusters = 5)
  s <- draw(pts, des, seed = 1, weights = TRUE)
  by_site <- ht_total(s, "y", by = "site")
  expect_false(inherits(by_site, "sf"))
  expect_false("geometry" %in% names(by_site))
  want <- ht_total(draw(plain, des, seed = 1, weights = TRUE), "y",
                   by = "site")
  expect_equal(as.data.frame(by_site), as.data.frame(want))
  expect_equal(sample_summary(s)$by_group,
               sample_summary(draw(plain, des, seed = 1,
                                   weights = TRUE))$by_group)
  st <- draw(pts, design_stratified("site", n = 12), seed = 1, weights = TRUE)
  expect_equal(sample_summary(st)$by_group$in_frame, c(30L, 20L, 10L))
})
