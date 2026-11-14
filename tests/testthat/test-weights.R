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
  # every duplicate at the distinct-unit weight, so the total is high by a
  # factor of (n / N) / (1 - (1 - 1 / N)^n). At n = 25 and N = 10 that factor
  # is 2.69, a total 169% high. ht_total() inherited the same bias.
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
  local_seed(11)
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
  expect_error(draw(df, design_certainty("w", 100,
                                         design_simple(n = 30, replace = TRUE)),
                    seed = 1, weights = TRUE),
               "Use\nreplace = FALSE in `rest`", fixed = TRUE)
})

test_that("a sample that carries a with-replacement design is refused", {
  # drawn 0.2.0 weighted a certainty design over a with-replacement `rest`,
  # so a sample it saved can hold duplicates weighted as distinct units.
  df <- data.frame(id = 1:40, g = rep(c("a", "b"), each = 20),
                   w = c(rep(1, 36), rep(500, 4)), y = 1:40)
  s <- draw(df, design_certainty("w", 100, design_simple(n = 10)), seed = 1,
            weights = TRUE)
  old <- s
  old$.row <- NULL
  attr(old, "drawn_pi") <- NULL
  attr(old, "drawn_design") <- design_certainty(
    "w", 100, design_simple(n = 30, replace = TRUE))
  msg <- paste0("`sample` comes from a certainty design over a ",
                "with-replacement `rest`, for which `weights = TRUE` is not ",
                "available.")
  expect_error(ht_total(old, "y"), msg, fixed = TRUE)
  expect_error(ht_mean(old, "y"), "replace = FALSE in `rest`", fixed = TRUE)
  expect_error(sample_summary(old), msg, fixed = TRUE)
  # A sample with its `.row` column is refused the same way.
  attr(s, "drawn_design") <- attr(old, "drawn_design")
  expect_error(ht_total(s, "y"), msg, fixed = TRUE)
  skip_if_not_installed("survey")
  expect_error(as_svydesign(old), msg, fixed = TRUE)
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

# ---- rebuilding, joining and compacting a sample -----------------------------

field_pop <- function() {
  data.frame(id = 1:120, site = rep(c("a", "b", "c"), times = c(60, 40, 20)),
             cl = rep(sprintf("c%02d", 1:12), each = 10),
             y = (1:120 * 37) %% 101, w = 5 + (1:120 * 13) %% 50,
             note = rep(strrep("x", 30), 120), stringsAsFactors = FALSE)
}

test_that("as_drawn() rebuilds a sample read back from a CSV file", {
  pop <- field_pop()
  des <- design_stratified("site", n = 24)
  s <- draw(pop, des, seed = 1, weights = TRUE)
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  utils::write.csv(s, path, row.names = FALSE)
  back <- utils::read.csv(path)
  expect_error(ht_total(back, "y"), "as_drawn()", fixed = TRUE)

  # Through `.row`, with the rows in another order
  rebuilt <- as_drawn(back[rev(seq_len(nrow(back))), ], des, pop)
  expect_identical(names(rebuilt), names(s))
  expect_equal(ht_total(rebuilt, "y")[c("total", "se", "df")],
               ht_total(s, "y")[c("total", "se", "df")])
  expect_equal(ht_mean(rebuilt, "y", by = "site")$se,
               ht_mean(s, "y", by = "site")$se)

  # Through an id, with `.row` gone
  back$.row <- NULL
  by_id <- as_drawn(back, des, pop, by = "id")
  expect_identical(by_id$.row, s$.row)
  expect_equal(ht_total(by_id, "y")$se, ht_total(s, "y")$se)
  expect_equal(sample_summary(by_id)$n, 24L)
})

test_that("as_drawn() refuses ids that are missing, repeated or unknown", {
  pop <- field_pop()
  des <- design_stratified("site", n = 24)
  back <- as.data.frame(draw(pop, des, seed = 1, weights = TRUE))
  attributes(back) <- attributes(back)[c("names", "row.names", "class")]
  ids <- back[, c("id", "y")]

  bad <- ids
  bad$id[3] <- NA
  expect_error(as_drawn(bad, des, pop, by = "id"), "no value of `id` in 1 row")
  expect_error(as_drawn(rbind(ids, ids[2, ]), des, pop, by = "id"),
               "more than one row")
  bad <- ids
  bad$id[3] <- 999
  expect_error(as_drawn(bad, des, pop, by = "id"),
               "`id` 999, which `frame` does not have")
  twice <- rbind(pop, pop[pop$id == ids$id[1], ])
  expect_error(as_drawn(ids, des, twice, by = "id"),
               "`frame` has `id` .* in more than one row")
  expect_error(as_drawn(rbind(back, back[1, ]), des, pop), "more than once")
  expect_error(as_drawn(ids, des, pop), "no `.row` column")
  expect_error(as_drawn(ids, des, pop, by = "site_id"), "`data` does not")

  # `.row` and `by` disagree when the frame was reordered since the draw
  expect_error(as_drawn(back, des, pop[nrow(pop):1, ], by = "id"),
               "disagrees with `id`")
  # A `.prob` that the design does not give
  expect_error(as_drawn(back, design_stratified("site", n = 30), pop),
               "`data` has 24 rows, and `design` draws 30")
  other <- back
  other$.prob <- other$.prob * 1.5
  expect_error(as_drawn(other, des, pop), "differ from the probabilities")
})

test_that("as_drawn() refuses the responding rows alone", {
  pop <- field_pop()
  des <- design_stratified("site", n = 24)
  s <- draw(pop, des, seed = 1, weights = TRUE)
  expect_error(as_drawn(as.data.frame(s)[-(1:3), ], des, pop),
               "`data` has 21 rows, and `design` draws 24")

  # Weighted conditional on the clusters: a cluster short of its take
  ms <- design_multistage("cl", n_clusters = 4, n = 18)
  m <- draw(pop, ms, seed = 2, weights = TRUE)
  expect_error(as_drawn(as.data.frame(m)[-1, ], ms, pop),
               "where `design` takes")
  bal <- design_cluster("cl", n_clusters = 3, balanced = TRUE)
  b <- draw(pop, bal, seed = 2, weights = TRUE)
  expect_error(as_drawn(as.data.frame(b)[b$cl != b$cl[1], ], bal, pop),
               "holds rows of 2 cluster")
})

test_that("as_drawn() rebuilds a sample weighted on its clusters", {
  pop <- field_pop()
  pop$cl <- rep(sprintf("c%02d", 1:12), times = rep(c(6, 10, 14), 4))
  for (des in list(design_multistage("cl", n_clusters = 5, n = 23),
                   design_cluster("cl", n_clusters = 4, balanced = TRUE))) {
    for (seed in 1:5) {
      s <- draw(pop, des, seed = seed, weights = TRUE)
      plain <- as.data.frame(s)
      attributes(plain) <- attributes(plain)[c("names", "row.names", "class")]
      back <- as_drawn(plain[rev(seq_len(nrow(plain))), ], des, pop)
      expect_equal(ht_total(back, "y")$se, ht_total(s, "y")$se)
    }
  }
})

test_that("as_drawn() keeps the class of `data`", {
  skip_if_not_installed("data.table")
  pop <- field_pop()
  des <- design_simple(n = 10)
  s <- draw(pop, des, seed = 1, weights = TRUE)
  dt <- data.table::as.data.table(as.data.frame(s)[, c(".row", "id", "y")])
  # `:=` works only where data.table is attached or imported, which the
  # global environment counts as.
  env <- new.env(parent = globalenv())
  env$back <- as_drawn(dt, des, pop)
  expect_s3_class(env$back, "data.table")
  expect_equal(ht_total(env$back, "y")$se, ht_total(s, "y")$se)
  expect_no_warning(evalq(back[, z := 1], env))
})

test_that("add_fieldwork() joins results in the sample's order", {
  pop <- field_pop()
  pop_y <- pop$y
  pop$y <- NULL
  s <- draw(pop, design_stratified("site", n = 24), seed = 1, weights = TRUE)
  field <- data.frame(id = rev(s$id), y = pop_y[rev(s$id)], ok = TRUE)
  joined <- add_fieldwork(s, field, by = "id")
  expect_identical(joined$.row, s$.row)
  expect_identical(names(joined), c(names(s), "y", "ok"))
  by_hand <- s
  by_hand$y <- pop_y[s$id]
  expect_identical(ht_total(joined, "y"), ht_total(by_hand, "y"))

  # A sampled row with no result gets NA
  part <- add_fieldwork(s, field[-1, ], by = "id")
  expect_identical(sum(is.na(part$y)), 1L)
  expect_identical(is.na(part$y), part$id == field$id[1])

  expect_error(add_fieldwork(s, rbind(field, field[1, ]), "id"),
               "in more than one row")
  expect_error(add_fieldwork(s, data.frame(id = c(s$id[1], 9999), y = 1),
                             "id"), "`id` 9999, which is not in `sample`")
  expect_error(add_fieldwork(joined, field, "id"),
               "`field` has `y` and `ok`, which `sample` already has")
  expect_error(add_fieldwork(s, data.frame(id = c(s$id[1], NA), y = 1), "id"),
               "no value of `id` in 1 row")
  expect_error(add_fieldwork(s, field["id"], "id"), "nothing to add")
  expect_error(add_fieldwork(as.data.frame(as.list(s)), field, "id"),
               "does not carry its design")
})

test_that("compact_sample() shrinks the saved sample and keeps the estimates", {
  pop <- field_pop()
  pop$hit <- as.numeric(pop$y > 50)
  s <- draw(pop, design_stratified("site", n = 24), seed = 1, weights = TRUE)
  small <- compact_sample(s)
  expect_lt(length(serialize(small, NULL)), length(serialize(s, NULL)) / 2)
  expect_identical(names(attr(small, "drawn_population")), "site")
  expect_identical(ht_total(small, "y"), ht_total(s, "y"))
  expect_identical(ht_mean(small, "y", by = "site"),
                   ht_mean(s, "y", by = "site"))
  # A 0/1 column of the frame is still read as a proportion
  expect_identical(ht_mean(small, "hit"), ht_mean(s, "hit"))
  expect_identical(ht_mean(small, "hit")$ci_scale, "logit")
  expect_equal(unclass(sample_summary(small)), unclass(sample_summary(s)))
  skip_if_not_installed("survey")
  expect_equal(survey::SE(survey::svytotal(~y, as_svydesign(small))),
               survey::SE(survey::svytotal(~y, as_svydesign(s))))
})

test_that("compact_sample() keeps what each design reads", {
  pop <- field_pop()
  pop$cs <- 10 + match(pop$cl, unique(pop$cl))
  designs <- list(
    design_stratified("site", n = 24, allocation = "neyman",
                      allocation_by = "w"),
    design_stratified("site", n = 24,
                      within = design_systematic(order_by = "w")),
    design_systematic(n = 20, order_by = "w"),
    design_cluster("cl", n_clusters = 4),
    design_multistage("cl", n_clusters = 4, n = 20),
    design_multistage("cl", n_clusters = 4, n = 12, size = "cs"),
    design_weighted("w", n = 20, method = "systematic", strata = "site"),
    design_certainty("w", 50, design_stratified("site", n = 12)),
    design_spread(c("y", "w"), n = 20, strata = "site")
  )
  for (des in designs) {
    s <- draw(pop, des, seed = 4, weights = TRUE)
    small <- compact_sample(s)
    expect_false("note" %in% names(attr(small, "drawn_population")))
    expect_identical(ht_total(small, "y"), ht_total(s, "y"))
  }
})

test_that("compact_sample() names `keep` where a dropped column is needed", {
  pop <- field_pop()
  pop$g <- rep(c("x", "y"), 60)
  s <- draw(pop, design_stratified("site", n = 24), seed = 1, weights = TRUE)
  small <- compact_sample(s)
  expect_error(ht_mean(small, "y", estimator = "ht", by = "g"),
               "compact_sample() dropped `g`", fixed = TRUE)
  expect_error(frame_totals(small, ~w), "Name it in `keep`")
  kept <- compact_sample(s, keep = c("g", "w"))
  expect_identical(ht_mean(kept, "y", estimator = "ht", by = "g"),
                   ht_mean(s, "y", estimator = "ht", by = "g"))
  expect_identical(frame_totals(kept, ~w), frame_totals(s, ~w))
  expect_error(compact_sample(small, keep = "g"),
               "An earlier compact_sample() dropped `g`", fixed = TRUE)
  expect_error(compact_sample(s, keep = "nope"), "is not a column")
})

test_that("compact_sample() keeps an sf frame's geometry only when read", {
  pts <- sf_points()
  box <- flat_box(c(-50, 50), c(-30, 30))
  on_geometry <- draw(pts, design_spatial("geometry", region = box, n = 10),
                      seed = 1, weights = TRUE)
  small <- compact_sample(on_geometry)
  expect_s3_class(attr(small, "drawn_population"), "sf")
  expect_identical(ht_total(small, "y"), ht_total(on_geometry, "y"))
  on_columns <- draw(pts, design_spatial(c("lon", "lat"), region = box,
                                         n = 10), seed = 1, weights = TRUE)
  small <- compact_sample(on_columns)
  expect_identical(names(attr(small, "drawn_population")), c("lon", "lat"))
  expect_identical(ht_total(small, "y"), ht_total(on_columns, "y"))
})
