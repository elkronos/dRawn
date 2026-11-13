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

  # Left alone, the one-row stratum gets nothing -- and the user is told,
  # because its row can then never be drawn and a total would omit it.
  expect_warning(res <- draw(skewed, design_stratified("g", n = 10), seed = 1),
                 "`rare`.*allocated no rows")
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
      got <- suppressWarnings(
        draw(d, design_stratified("site", n = nn, min_per_stratum = m),
             seed = 1),
        classes = "drawn_empty_stratum")
      expect_equal(nrow(got), nn,
                   info = paste("min_per_stratum =", m, "n =", nn))
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

test_that("a tie goes to the stratum with more room, then first in byte order", {
  # Equal allocation of 8 over three strata leaves 2 rows after the whole
  # parts, with equal remainders. "b" has the most room; "B" and "a" tie on
  # room too, and "B" comes first in byte order in every locale.
  d <- data.frame(id = 1:54, g = rep(c("a", "B", "b"), times = c(12, 12, 30)))
  des <- design_stratified("g", n = 8, allocation = "equal")
  want <- ifelse(d$g == "B", 3 / 12, ifelse(d$g == "a", 2 / 12, 3 / 30))
  expect_equal(inclusion_prob(d, des), want)
  ids <- draw(d, des, seed = 1)$id
  with_other_collation({
    expect_equal(inclusion_prob(d, des), want)
    expect_identical(draw(d, des, seed = 1)$id, ids)
  })

  # A factor's level order decides instead.
  d$g <- factor(d$g, levels = c("a", "B", "b"))
  expect_equal(inclusion_prob(d, des),
               ifelse(d$g == "a", 3 / 12, ifelse(d$g == "B", 2 / 12, 3 / 30)))
})

test_that("a sample drawn in one collation estimates the same in another", {
  d <- data.frame(id = 1:54, g = rep(c("a", "B", "b"), times = c(12, 12, 30)),
                  y = 100 + 30 * rep(1:3, times = c(12, 12, 30)) + (1:54 %% 7))
  s <- draw(d, design_stratified("g", n = 8, allocation = "equal"), seed = 1,
            weights = TRUE)
  want <- ht_total(s, "y")
  with_other_collation({
    got <- ht_total(s, "y")
    expect_equal(got$se, want$se)
    expect_equal(got$df, want$df)
  })
})

test_that("an sf frame draws the same strata as the plain frame", {
  pts <- sf_points()
  plain <- sf::st_drop_geometry(pts)
  des <- design_stratified("site", n = 12)
  s <- draw(pts, des, seed = 1, weights = TRUE)
  expect_s3_class(s, "sf")
  expect_identical(s$id, draw(plain, des, seed = 1)$id)
  expect_equal(inclusion_prob(pts, des), inclusion_prob(plain, des))
  expect_equal(ht_total(s, "y")$se,
               ht_total(draw(plain, des, seed = 1, weights = TRUE), "y")$se)
})

# ---- floors and caps -------------------------------------------------------

test_that("a floor holds the allocation rule among the strata it leaves free", {
  # Raising a stratum to min_per_stratum took the rows back one at a time
  # from the others, which broke the rule between them. Under Neyman
  # allocation `a` and `b` got 19 and 101 rows where the rule gives 35.3 and
  # 84.7, which cost 14% in variance.
  bounded <- function(n, share, lower, upper) {
    f <- function(l) sum(pmin(pmax(l * share, lower), upper)) - n
    l <- stats::uniroot(f, c(0, max(upper / share)), tol = 1e-12)$root
    pmin(pmax(l * share, lower), upper)
  }
  sizes <- c(a = 5000, b = 4000,
             stats::setNames(rep(10, 20), paste0("s", 1:20)))
  sp <- c(1, 3, rep(1, 20))
  got <- drawn:::allocate(200, sizes, "neyman", 4, spread = sp)
  want <- bounded(200, sizes * sp, pmin(4, sizes), sizes)
  expect_equal(sum(got), 200L)
  expect_true(all(abs(got - want) < 1))
  expect_equal(unname(got), c(35L, 85L, rep(4L, 20)))

  got <- drawn:::allocate(200, sizes, "proportional", 4)
  expect_true(all(abs(got - bounded(200, sizes, pmin(4, sizes), sizes)) < 1))
  expect_equal(unname(got[1:2]), c(67L, 53L))

  # The same through a frame, on both the draw and the probabilities
  d <- data.frame(g = rep(c("a", "b", paste0("s", 1:20)),
                          times = c(500, 400, rep(10, 20))))
  des <- design_stratified("g", n = 200, min_per_stratum = 4)
  p <- inclusion_prob(d, des)
  expect_equal(p[match(c("a", "b", "s1"), d$g)], c(67 / 500, 53 / 400, 0.4))
  expect_equal(c(table(draw(d, des, seed = 1)$g)[c("a", "b")]),
               c(a = 67L, b = 53L))
})

test_that("a floor larger than a stratum asks only for that stratum's rows", {
  # Each of the three one-row strata is taken whole, so 7 rows cover a floor
  # of 2. The check counted 2 for every stratum and refused anything under 8.
  d <- data.frame(id = 1:103, g = c("a", "b", "c", rep("d", 100)))
  got <- draw(d, design_stratified("g", n = 7, min_per_stratum = 2), seed = 1)
  expect_equal(as.vector(table(got$g)), c(1L, 1L, 1L, 4L))
  expect_error(
    draw(d, design_stratified("g", n = 4, min_per_stratum = 2), seed = 1),
    paste0("`min_per_stratum` = 2 needs at least 5 rows across these 4 ",
           "strata, with the 3 strata smaller than 2 taken whole, but `n` ",
           "is 4."), fixed = TRUE)
  # With replacement a one-row stratum can give the floor in repeats.
  expect_error(
    inclusion_prob(d, design_stratified("g", n = 7, min_per_stratum = 2,
                                        replace = TRUE)),
    "needs at least 8 rows", fixed = TRUE)
})

test_that("Neyman refuses an infinite allocation_by", {
  # sd() of a stratum holding Inf is NaN, which was read as no spread: the
  # stratum with the most extreme values was allocated nothing.
  d <- data.frame(id = 1:60, g = rep(c("a", "b", "c"), each = 20),
                  v = c(1:20 * 50, (1:40 * 7) %% 13))
  d$v[3] <- Inf
  des <- design_stratified("g", n = 15, allocation = "neyman",
                           allocation_by = "v")
  expect_error(draw(d, des, seed = 1), "1 value(s) of `v` are missing or non-finite",
               fixed = TRUE)
  expect_error(inclusion_prob(d, des), "missing or non-finite")

  # Finite values whose standard deviation overflows split as the same
  # values on a smaller scale do.
  d$v[3] <- 150
  big <- d
  big$v <- big$v * 1e305
  expect_equal(inclusion_prob(big, des), inclusion_prob(d, des))
})

test_that("strata left with no spread after a cap are split by size", {
  # The Neyman shares are 50, 0 and 0. `a` is taken whole, and the 25 rows
  # left went to `b` and `c` equally, 12 and 13. A Neyman split over strata
  # with no spread goes by size, 6 and 19, as it does when no stratum has any.
  d <- data.frame(g = rep(c("a", "b", "c"), times = c(5, 100, 300)),
                  v = c(1:5 * 10, rep(3, 400)))
  des <- design_stratified("g", n = 30, allocation = "neyman",
                           allocation_by = "v")
  p <- inclusion_prob(d, des)
  expect_equal(p[c(1, 6, 106)], c(1, 6 / 100, 19 / 300))
})

test_that("the empty-stratum warning gives the cause it has", {
  # A Neyman stratum with constant values gets rows only once every other
  # stratum is taken whole, so "Raise `n`" pointed the wrong way.
  z <- data.frame(id = 1:60, g = rep(c("a", "b", "c"), each = 20),
                  v = c(rep(5, 20), (1:40 * 7) %% 13))
  des <- design_stratified("g", n = 12, allocation = "neyman",
                           allocation_by = "v")
  w <- tryCatch(draw(z, des, seed = 1), drawn_empty_stratum = function(w) w)
  expect_s3_class(w, "drawn_warning")
  expect_match(conditionMessage(w), paste0(
    "`v` is constant within it, so Neyman allocation gives it no share ",
    "until every other stratum is taken whole. Set min_per_stratum = 1"),
    fixed = TRUE)
  expect_no_match(conditionMessage(w), "Raise `n`", fixed = TRUE)

  # A stratum with spread but too small a share can still be reached by a
  # larger `n`, and the message keeps that advice for it.
  mixed <- data.frame(g = rep(c("a", "b", "c"), times = c(20, 2, 200)),
                      v = c(rep(5, 20), 1, 2, 1:200 %% 10))
  w <- tryCatch(
    draw(mixed, design_stratified("g", n = 10, allocation = "neyman",
                                  allocation_by = "v"), seed = 1),
    drawn_empty_stratum = function(w) w)
  expect_match(conditionMessage(w), "2 strata (`a`, `b`)", fixed = TRUE)
  expect_match(conditionMessage(w), "`v` is constant within `a`, so",
               fixed = TRUE)
  expect_match(conditionMessage(w), "Raise `n`, or set min_per_stratum = 1",
               fixed = TRUE)

  # n = 0 asks for no rows, and gets none without a warning per stratum, as
  # design_simple(n = 0) does.
  d <- data.frame(id = 1:10, g = rep(c("a", "b"), each = 5))
  expect_no_warning(got <- draw(d, design_stratified("g", n = 0), seed = 1))
  expect_equal(nrow(got), 0L)
})

# ---- counts per stratum ----------------------------------------------------

counts_frame <- function() {
  data.frame(id = 1:100,
             site = rep(c("north", "south", "east", "west"),
                        times = c(40, 30, 20, 10)),
             y = (1:100 * 37) %% 23 + 1)
}

test_that("counts given per stratum are drawn as given, on every path", {
  d <- counts_frame()
  counts <- c(north = 4, south = 6, east = 5, west = 5)
  des <- design_stratified("site", allocation = counts)
  expect_identical(des$n, 20L)
  expect_identical(design_stratified("site", n = 20, allocation = counts)$n,
                   20L)
  expect_output(print(des), "c(north = 4, south = 6, east = 5, west = 5)",
                fixed = TRUE)

  s <- draw(d, des, seed = 1, weights = TRUE)
  expect_equal(c(table(s$site))[names(counts)], counts)
  N_h <- c(table(d$site))
  want <- unname(counts[d$site] / N_h[d$site])
  expect_equal(inclusion_prob(d, des), want)
  expect_equal(s$.prob, want[s$.row])
  expect_equal(joint_prob(d, des, rows = c(1, 2, 41)),
               matrix(c(0.1, 0.1 * 3 / 39, 0.1 * 0.2,
                        0.1 * 3 / 39, 0.1, 0.1 * 0.2,
                        0.1 * 0.2, 0.1 * 0.2, 0.2), 3))

  # The variance is the stratified form over the counts as given.
  v <- sum(vapply(names(counts), function(h) {
    yh <- s$y[s$site == h]
    N_h[[h]]^2 * (1 - counts[[h]] / N_h[[h]]) * stats::var(yh) / counts[[h]]
  }, numeric(1)))
  expect_equal(ht_total(s, "y")$variance, v)
  expect_equal(ht_total(s, "y", variance = "jackknife")$variance, v)
  sm <- sample_summary(s)$by_group
  expect_equal(sm$drawn[match(names(counts), sm$group)], unname(counts))

  skip_if_not_installed("survey")
  expect_equal(
    unname(survey::SE(survey::svytotal(~y, as_svydesign(s)))[1]),
    ht_total(s, "y")$se)
})

test_that("counts for several strata columns come as a data frame", {
  d <- data.frame(id = 1:40, g = rep(c("a", "b"), 20),
                  h = rep(c(1.5, 2), each = 20))
  spec <- data.frame(h = c(1.5, 2, 1.5, 2), g = c("a", "a", "b", "b"),
                     n = c(2, 3, 4, 1))
  des <- design_stratified(c("g", "h"), allocation = spec)
  expect_identical(des$n, 10L)
  got <- draw(d, des, seed = 3)
  expect_equal(c(table(paste(got$g, got$h))),
               c("a 1.5" = 2L, "a 2" = 3L, "b 1.5" = 4L, "b 2" = 1L))
  expect_equal(inclusion_prob(d, des)[c(1, 2, 21, 22)],
               c(2, 4, 3, 1) / 10)

  # A factor column is matched on its values, not its level order.
  f <- counts_frame()
  f$site <- factor(f$site, levels = c("west", "east", "south", "north"))
  des <- design_stratified("site",
                           allocation = c(north = 4, south = 6, east = 5,
                                          west = 5))
  expect_equal(c(table(draw(f, des, seed = 1)$site)),
               c(west = 5L, east = 5L, south = 6L, north = 4L))
})

test_that("counts that do not fit the frame are refused by name", {
  d <- counts_frame()
  expect_error(
    draw(d, design_stratified("site", allocation = c(north = 4, south = 6,
                                                     east = 5)), seed = 1),
    "`allocation` has no entry for 1 stratum of `data`: `west`.",
    fixed = TRUE)
  expect_error(
    inclusion_prob(d, design_stratified(
      "site", allocation = c(north = 4, south = 6, east = 5, west = 5,
                             nrth = 1))),
    "`allocation` names 1 stratum that `data` does not have: `nrth`.",
    fixed = TRUE)
  over <- design_stratified("site", allocation = c(north = 4, south = 6,
                                                   east = 5, west = 11))
  expect_error(draw(d, over, seed = 1),
               "Stratum `west` has 10 row(s) but `allocation` asks for 11.",
               fixed = TRUE)
  expect_error(inclusion_prob(d, over), "asks for 11", fixed = TRUE)
  over$replace <- TRUE
  expect_equal(nrow(draw(d, over, seed = 1)), 26L)

  # A count of 0 leaves the stratum unreachable, and says so.
  expect_warning(
    draw(d, design_stratified("site", allocation = c(north = 4, south = 0,
                                                     east = 5, west = 5)),
         seed = 1),
    "Give each stratum a count of at least 1 in `allocation`",
    class = "drawn_empty_stratum")
})

test_that("counts are checked when the design is built", {
  four <- c(north = 4, south = 6, east = 5, west = 5)
  expect_error(design_stratified("site", n = 21, allocation = four),
               "`n` (21) is not the sum of the counts in `allocation` (20)",
               fixed = TRUE)
  expect_error(design_stratified("site", allocation = four,
                                 min_per_stratum = 2),
               "count of at least 2", fixed = TRUE)
  expect_error(design_stratified("site", allocation = c(a = 2.5, b = 4)),
               "stratum `a` has 2.5", fixed = TRUE)
  expect_error(design_stratified("site", allocation = c(a = -1, b = 4)),
               "whole numbers of 0 or more")
  expect_error(design_stratified("site", allocation = c(4, 6)),
               "must be named by stratum")
  expect_error(design_stratified("site", allocation = c(a = 1, a = 4)),
               "gives stratum `a` more than once", fixed = TRUE)
  expect_error(design_stratified(c("g", "h"), allocation = c(a = 1)),
               "a data frame with columns `g`, `h` and `n`", fixed = TRUE)
  expect_error(design_stratified("g", allocation = data.frame(g = "a",
                                                              count = 3)),
               "missing column: `n`", fixed = TRUE)
  expect_error(design_stratified("site"), "`n` is missing", fixed = TRUE)
  # Counts passed as `n` are pointed to `allocation`.
  expect_error(design_stratified("site", n = four),
               "Give the count for each stratum as `allocation`",
               fixed = TRUE)
})

# ---- standard deviations and costs per stratum -----------------------------

sd_frame <- function() {
  data.frame(id = 1:120, g = rep(c("a", "b", "c"), times = c(60, 40, 20)),
             v = c(1:60 %% 5, (1:40 * 7) %% 23, (1:20 * 11) %% 37),
             y = (1:120 * 13) %% 29 + 1)
}

test_that("Neyman from standard deviations per stratum matches a column", {
  d <- sd_frame()
  by_col <- design_stratified("g", n = 24, allocation = "neyman",
                              allocation_by = "v", min_per_stratum = 2)
  sds <- vapply(split(d$v, d$g), stats::sd, numeric(1))
  by_vec <- design_stratified("g", n = 24, allocation = "neyman",
                              allocation_by = sds, min_per_stratum = 2)
  by_df <- design_stratified(
    "g", n = 24, allocation = "neyman", min_per_stratum = 2,
    allocation_by = data.frame(g = c("c", "a", "b"),
                               sd = unname(sds[c("c", "a", "b")])))
  expect_output(print(by_vec), "allocation_by    c(a = ", fixed = TRUE)

  # Given the standard deviations, the frame needs no `v` column, on any
  # path.
  bare <- d[c("id", "g", "y")]
  p <- inclusion_prob(d, by_col)
  expect_equal(inclusion_prob(bare, by_vec), p)
  expect_equal(inclusion_prob(bare, by_df), p)
  expect_identical(draw(bare, by_vec, seed = 1)$id,
                   draw(d, by_col, seed = 1)$id)
  rows <- c(1, 2, 61, 101)
  expect_equal(joint_prob(bare, by_vec, rows = rows),
               joint_prob(d, by_col, rows = rows))
  s <- draw(bare, by_vec, seed = 1, weights = TRUE)
  ref <- draw(d, by_col, seed = 1, weights = TRUE)
  expect_equal(ht_total(s, "y")$se, ht_total(ref, "y")$se)
  expect_equal(ht_total(s, "y", variance = "jackknife")$se,
               ht_total(ref, "y", variance = "jackknife")$se)

  skip_if_not_installed("survey")
  expect_equal(
    unname(survey::SE(survey::svytotal(~y, as_svydesign(s)))[1]),
    ht_total(s, "y")$se)
})

test_that("costs per stratum give the cost-optimal split", {
  # Shares proportional to N_h * S_h / sqrt(c_h), rounded by largest
  # remainder: 27.71, 25.98, 6.93 and 29.39 rows.
  sizes <- c(a = 400, b = 300, c = 200, d = 100)
  S <- c(2, 5, 3, 12)
  cst <- c(1, 4, 9, 2)
  share <- sizes * S / sqrt(cst)
  want <- 90 * share / sum(share)
  got <- drawn:::allocate(90, sizes, "neyman", 0, spread = S, cost = cst)
  expect_equal(sum(got), 90L)
  expect_true(all(abs(got - want) < 1))
  expect_equal(unname(got), c(28L, 26L, 7L, 29L))

  # The same through a frame, with the values in one data frame
  d <- data.frame(g = rep(names(sizes), times = sizes))
  pilot <- data.frame(g = names(sizes), sd = S, cost = cst)
  des <- design_stratified("g", n = 90, allocation = "neyman",
                           allocation_by = pilot, cost = pilot)
  expect_equal(inclusion_prob(d, des),
               unname((got / sizes)[d$g]))
  expect_equal(c(table(draw(d, des, seed = 1)$g)), got)
  tab <- allocation_table(d, des)
  expect_equal(tab$n, unname(got))
  expect_equal(tab$sd, S)
  expect_equal(tab$cost, cst)

  # Costs with a column: the spread is the column's within each stratum.
  f <- sd_frame()
  sds <- vapply(split(f$v, f$g), stats::sd, numeric(1))
  costs <- c(a = 1, b = 4, c = 9)
  by_col <- design_stratified("g", n = 24, allocation = "neyman",
                              allocation_by = "v", cost = costs)
  expect_equal(allocation_table(f, by_col)$n,
               unname(drawn:::allocate(24, c(60, 40, 20), "neyman", 0,
                                       spread = sds, cost = costs)))

  # With no spread anywhere, the shares are those for equal spreads.
  expect_equal(drawn:::allocate(30, c(100, 100), "neyman", 0,
                                spread = c(0, 0), cost = c(1, 4)),
               c(20L, 10L))
})

test_that("standard deviations and costs are checked at construction", {
  ney <- function(...) {
    design_stratified("g", n = 10, allocation = "neyman", ...)
  }
  expect_error(ney(allocation_by = c(a = 1, b = -2)),
               paste0("The standard deviations in `allocation_by` must be ",
                      "finite numbers of 0 or more, and stratum `b` has -2."),
               fixed = TRUE)
  expect_error(ney(allocation_by = c(a = 1, b = Inf)), "stratum `b` has Inf",
               fixed = TRUE)
  expect_error(ney(allocation_by = c(a = 1, b = NA)), "stratum `b` has NA",
               fixed = TRUE)
  expect_error(ney(allocation_by = c(1, 2)), "must be named by stratum")
  expect_error(ney(allocation_by = TRUE),
               "`allocation_by` must be a column name, or the standard",
               fixed = TRUE)
  expect_error(ney(allocation_by = data.frame(g = c("a", "b"),
                                              sd = c("1", "2"))),
               "must be numbers, not character", fixed = TRUE)
  expect_error(ney(allocation_by = data.frame(g = "a", s = 1)),
               "missing column: `sd`", fixed = TRUE)
  expect_error(ney(allocation_by = "v", cost = c(a = 1, b = 0)),
               paste0("The costs in `cost` must be finite numbers above 0, ",
                      "and stratum `b` has 0."), fixed = TRUE)
  expect_error(ney(allocation_by = "v", cost = c(a = 1, b = -3)),
               "stratum `b` has -3", fixed = TRUE)
  expect_error(design_stratified(c("g", "h"), n = 10, allocation = "neyman",
                                 allocation_by = "v", cost = c(a = 1)),
               "a data frame with columns `g`, `h` and `cost`", fixed = TRUE)

  # Costs shape a Neyman split and nothing else.
  expect_error(design_stratified("g", n = 10, cost = c(a = 1, b = 2)),
               paste0("`cost` applies only to allocation = \"neyman\", and ",
                      "`allocation` is \"proportional\"."), fixed = TRUE)
  expect_error(design_stratified("g", allocation = c(a = 1, b = 2),
                                 cost = c(a = 1, b = 2)),
               "`allocation` gives the count for each stratum", fixed = TRUE)
})

test_that("standard deviations and costs that miss a stratum are refused", {
  d <- sd_frame()
  short <- design_stratified("g", n = 12, allocation = "neyman",
                             allocation_by = c(a = 1, b = 2))
  msg <- "`allocation_by` has no entry for 1 stratum of `data`: `c`."
  expect_error(draw(d, short, seed = 1), msg, fixed = TRUE)
  expect_error(inclusion_prob(d, short), msg, fixed = TRUE)
  expect_error(allocation_table(d, short), msg, fixed = TRUE)
  expect_error(
    inclusion_prob(d, design_stratified(
      "g", n = 12, allocation = "neyman",
      allocation_by = c(a = 1, b = 2, c = 3, cc = 4))),
    "`allocation_by` names 1 stratum that `data` does not have: `cc`.",
    fixed = TRUE)
  expect_error(
    draw(d, design_stratified("g", n = 12, allocation = "neyman",
                              allocation_by = "v", cost = c(a = 1, c = 2)),
         seed = 1),
    "`cost` has no entry for 1 stratum of `data`: `b`.", fixed = TRUE)

  # A standard deviation of 0 leaves a stratum without a share, and the
  # warning says where that came from.
  zero <- design_stratified("g", n = 12, allocation = "neyman",
                            allocation_by = c(a = 0, b = 2, c = 3))
  w <- tryCatch(draw(d, zero, seed = 1), drawn_empty_stratum = function(w) w)
  expect_match(conditionMessage(w), paste0(
    "`allocation_by` gives it a standard deviation of 0, so Neyman ",
    "allocation gives it no share"), fixed = TRUE)
})

# ---- allocation_table() ----------------------------------------------------

status_frame <- function() {
  data.frame(id = 1:262,
             site = rep(c("a", "b", "c", "d", "e"), c(150, 80, 25, 5, 2)),
             v = c(1:150 %% 7, (1:80 * 7) %% 19, rep(1:2, length.out = 25),
                   c(1, 90, 3, 70, 40), c(4, 4)))
}

test_that("allocation_table() agrees with the draw and inclusion_prob()", {
  d <- status_frame()
  des <- design_stratified("site", n = 30, allocation = "neyman",
                           allocation_by = "v", min_per_stratum = 3)
  tab <- allocation_table(d, des)
  expect_identical(names(tab),
                   c("site", "N", "n", "prob", "sd", "cost", "status"))
  expect_identical(tab$site, c("a", "b", "c", "d", "e"))
  expect_equal(tab$N, c(150L, 80L, 25L, 5L, 2L))
  expect_equal(tab$n, as.vector(table(draw(d, des, seed = 1)$site)))
  expect_equal(tab$prob,
               as.vector(tapply(inclusion_prob(d, des), d$site, unique)))
  expect_equal(tab$sd, as.vector(tapply(d$v, d$site, stats::sd)))
  expect_true(all(is.na(tab$cost)))
  # `d` is small and varies most, so the rule asks more of it than its 5
  # rows; `e` is smaller than the floor; `c` barely varies.
  expect_identical(tab$status, c("rule", "rule", "floor", "capped",
                                 "take-all"))

  # Without the floor, `c`'s share rounds to nothing and `e` has no spread.
  bare <- allocation_table(d, design_stratified(
    "site", n = 30, allocation = "neyman", allocation_by = "v"))
  expect_identical(bare$status, c("rule", "rule", "empty", "capped", "empty"))
  expect_equal(bare$prob[bare$status == "empty"], c(0, 0))

  # With replacement no stratum is capped or taken whole.
  rep_des <- design_stratified("site", n = 30, allocation = "neyman",
                               allocation_by = "v", min_per_stratum = 3,
                               replace = TRUE)
  tab <- allocation_table(d, rep_des)
  expect_identical(tab$status, c("rule", "rule", "floor", "rule", "floor"))
  expect_equal(tab$prob,
               as.vector(tapply(inclusion_prob(d, rep_des), d$site, unique)))

  # A census takes every stratum whole without capping any.
  expect_identical(
    allocation_table(d, design_stratified("site", n = 262))$status,
    rep("take-all", 5))

  # Counts as given
  cnt <- allocation_table(d, design_stratified(
    "site", allocation = c(a = 3, b = 0, c = 25, d = 2, e = 1)))
  expect_identical(cnt$status, c("count", "empty", "take-all", "count",
                                 "count"))
  expect_true(all(is.na(cnt$sd)))

  # Values whose standard deviation overflows inside sd() are split on a
  # smaller scale; the table reports the spread on the column's own scale.
  big <- d
  big$v <- big$v * 1e305
  expect_equal(allocation_table(big, des)$sd,
               allocation_table(d, des)$sd * 1e305)
})

test_that("allocation_table() keeps the strata columns and their order", {
  d <- data.frame(id = 1:40, g = factor(rep(c("b", "a"), 20),
                                        levels = c("b", "a")),
                  h = rep(c(1.5, 2), each = 20))
  tab <- allocation_table(d, design_stratified(c("g", "h"), n = 8))
  expect_identical(tab$g, factor(c("b", "b", "a", "a"), levels = c("b", "a")))
  expect_identical(tab$h, c(1.5, 2, 1.5, 2))
  expect_equal(tab$n, rep(2L, 4))
})

test_that("allocation_table() describes the rest of a certainty design", {
  d <- status_frame()
  d$w <- c(rep(1, 250), seq(100, 1200, by = 100))
  rest <- design_stratified("site", n = 20, min_per_stratum = 2)
  des <- design_certainty("w", threshold = 100, rest = rest)
  tab <- allocation_table(d, des)
  below <- d$w < 100
  expect_equal(tab$N, as.vector(table(d$site[below])))
  expect_equal(tab$prob, as.vector(tapply(inclusion_prob(d, des)[below],
                                          d$site[below], unique)))
  expect_equal(tab, allocation_table(d[below, ], rest))

  # With every row at or above the threshold, `rest` has no strata.
  none <- allocation_table(d, design_certainty("w", threshold = 0,
                                               rest = rest))
  expect_equal(nrow(none), 0L)
  expect_identical(names(none), names(tab))
})

test_that("allocation_table() refuses what it cannot describe", {
  d <- status_frame()
  expect_error(allocation_table(d, design_simple(n = 3)),
               "which has no strata to allocate over", fixed = TRUE)
  expect_error(allocation_table(d, list()),
               "`design` must come from one of the design_*() constructors",
               fixed = TRUE)
  clash <- data.frame(n = rep(c("x", "y"), 5))
  expect_error(allocation_table(clash, design_stratified("n", n = 4)),
               paste0("`strata` column `n` has the name of a column ",
                      "allocation_table() adds."), fixed = TRUE)
})
