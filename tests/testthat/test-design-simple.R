test_that("draws the requested number of rows", {
  df <- make_df()
  expect_equal(nrow(draw(df, design_simple(n = 10), seed = 123)), 10)
})

test_that("refuses to over-draw without replacement, but allows it with", {
  df <- make_df()
  expect_error(draw(df, design_simple(n = 101), seed = 1),
               "cannot exceed the number of rows")
  expect_equal(nrow(draw(df, design_simple(n = 150, replace = TRUE), seed = 1)),
               150)
})

test_that("n is validated at construction time", {
  expect_error(design_simple(n = -1), "must be non-negative")
  expect_error(design_simple(n = 2.7), "single whole number")
  expect_error(design_simple(n = NA), "single whole number")
  expect_error(design_simple(n = c(5, 10)), "single whole number")
  expect_equal(design_simple(n = 0)$n, 0L)
})

test_that("replace is validated", {
  expect_error(design_simple(n = 5, replace = "yes"), "must be TRUE or FALSE")
  expect_error(design_simple(n = 5, replace = NA), "must be TRUE or FALSE")
})

test_that("n = 0 draws nothing but keeps the schema", {
  df <- make_df()
  res <- draw(df, design_simple(n = 0), seed = 1)
  expect_equal(nrow(res), 0L)
  expect_same_schema(res, df)
})

test_that("a tibble comes back a tibble", {
  skip_if_not_installed("tibble")
  tb <- tibble::as_tibble(make_df())
  expect_s3_class(draw(tb, design_simple(n = 3), seed = 1), "tbl_df")
})

# ---- reserves -----------------------------------------------------------------

test_that("reserve is validated at construction time", {
  expect_equal(design_simple(n = 5)$reserve, 0L)
  expect_equal(design_simple(n = 5, reserve = 3)$reserve, 3L)
  expect_error(design_simple(n = 5, reserve = -1), "must be non-negative")
  expect_error(design_simple(n = 5, reserve = 1.5), "single whole number")
  expect_error(design_simple(n = 5, replace = TRUE, reserve = 2),
               "replace = FALSE with `reserve`")
  expect_error(design_certainty("value", 1, rest = design_simple(5, reserve = 2)),
               "Leave `reserve` out of `rest`")
  df <- make_df()
  expect_error(draw(df, design_simple(n = 95, reserve = 6)),
               "`n` and `reserve` together \\(101\\) cannot exceed")
  expect_error(inclusion_prob(df, design_simple(n = 95, reserve = 6)),
               "`n` and `reserve` together")
})

test_that("a list holds the main sample, then reserves disjoint from it", {
  df <- make_df()
  des <- design_simple(n = 10, reserve = 4)
  s <- draw(df, des, seed = 7)
  expect_equal(names(s), c(".reserve", names(df)))
  expect_identical(s$.reserve, c(integer(10), 1:4))
  expect_false(anyDuplicated(s$id) > 0)
  # The main sample is the one the design without reserves draws.
  expect_identical(s$id[1:10], draw(df, design_simple(n = 10), seed = 7)$id)

  w <- draw(df, des, seed = 7, weights = TRUE)
  expect_equal(names(w)[1:5], c(".prob", ".weight", ".row", ".reserve", "id"))
  expect_equal(w$.prob, rep(0.1, 14))
  expect_equal(inclusion_prob(df, des), rep(0.1, 100))
  expect_equal(attr(w, "drawn_rows"), w$.row)

  clash <- df
  clash$.reserve <- 1
  expect_error(draw(clash, des), "already has a column called `.reserve`")
  expect_equal(nrow(draw(clash, design_simple(n = 3), seed = 1)), 3L)
})

test_that("estimates refuse a list until its reserves are resolved", {
  df <- make_df()
  s <- draw(df, design_simple(n = 10, reserve = 4), seed = 1, weights = TRUE)
  expect_error(ht_total(s, "value"), "use_reserves\\(\\)")
  expect_error(ht_mean(s, "value"), "not yet a sample to estimate from")
  expect_error(sample_summary(s), "use_reserves\\(\\)")
  skip_if_not_installed("survey")
  expect_error(as_svydesign(s), "use_reserves\\(\\)")
})

test_that("use_reserves() takes the first available reserves, in order", {
  df <- make_df()
  s <- draw(df, design_simple(n = 10, reserve = 5), seed = 3, weights = TRUE)
  # Units 2 and 6 of the main sample are out of reach, and so is reserve 1;
  # reserves 2 and 3 replace them, and 4 and 5 are never visited.
  gone <- s$id[c(2, 6, 11)]
  found <- ifelse(s$.reserve <= 3, s$id %in% gone, NA)
  u <- use_reserves(s, found)
  expect_equal(nrow(u), 10L)
  expect_equal(u$.reserve, c(integer(8), 2L, 3L))
  expect_equal(u$id, s$id[c(1, 3:5, 7:10, 12, 13)])
  expect_equal(u$.prob, rep(10 / 97, 10))
  expect_equal(u$.weight, 1 / u$.prob)

  des <- attr(u, "drawn_design")
  expect_equal(des$reserve, 0L)
  expect_equal(des$unavailable, sort(s$.row[c(2, 6, 11)]))
  p <- inclusion_prob(df, des)
  expect_equal(p[gone], c(0, 0, 0))
  expect_equal(sum(p), 10)

  # A simple random sample of 10 from the 97 rows left.
  r <- ht_total(u, "value")
  expect_equal(r$se, 97 * sqrt((1 - 10 / 97) * stats::var(u$value) / 10))
  expect_equal(r$unreachable, 3L)
  expect_equal(sample_summary(u)$N_reachable, 97L)
  expect_equal(ht_mean(u, "value")$mean, mean(u$value))
  skip_if_not_installed("survey")
  svd <- as_svydesign(u)
  expect_equal(as.numeric(survey::SE(survey::svytotal(~value, svd))), r$se)
})

test_that("use_reserves() reads a column, keeps the frame copy compact", {
  df <- make_df()
  s <- draw(df, design_simple(n = 6, reserve = 3), seed = 2, weights = TRUE)
  field <- data.frame(id = s$id, gone = c(TRUE, rep(FALSE, 6), NA, NA))
  s <- add_fieldwork(s, field, by = "id")
  u <- use_reserves(s, "gone")
  expect_equal(u$id, s$id[c(2:7)])
  small <- use_reserves(compact_sample(s), "gone")
  expect_false(is.null(attr(small, "drawn_compact")))
  expect_equal(ht_total(small, "value")$se, ht_total(u, "value")$se)

  # A reordered list is matched through `.row`.
  shuffled <- s[rev(seq_len(nrow(s))), ]
  expect_equal(sort(use_reserves(shuffled, "gone")$id), sort(u$id))
})

test_that("use_reserves() refuses findings it cannot place", {
  df <- make_df()
  s <- draw(df, design_simple(n = 6, reserve = 3), seed = 2, weights = TRUE)
  expect_error(use_reserves(s, rep(NA, 9)),
               "missing for 6 unit\\(s\\) of the main sample")
  expect_error(use_reserves(s, c(TRUE, rep(FALSE, 5), NA, NA, NA)),
               "missing for reserve 1, at frame row")
  expect_error(use_reserves(s, 1:9), "must be a logical vector")
  expect_error(use_reserves(s, logical(4)), "one value for each of the 9")
  expect_error(use_reserves(s, "nope"), "no column `nope`")
  expect_error(use_reserves(s, rep(TRUE, 9)), "Every unit `sample` lists")
  expect_error(use_reserves(s[1:6, ], logical(6)), "6 of the 9 rows")
  edited <- s
  edited$.reserve <- rev(edited$.reserve)
  expect_error(use_reserves(edited, logical(9)), "`.reserve` was changed")
  plain <- draw(df, design_simple(n = 6), seed = 2, weights = TRUE)
  expect_error(use_reserves(plain, logical(6)), "lists no reserves")
})

test_that("use_reserves() warns when the reserves run out", {
  df <- make_df()
  s <- draw(df, design_simple(n = 6, reserve = 2), seed = 2, weights = TRUE)
  found <- c(TRUE, TRUE, TRUE, rep(FALSE, 3), TRUE, FALSE)
  expect_warning(u <- use_reserves(s, found), class = "drawn_reserves_short")
  expect_equal(nrow(u), 4L)
  expect_equal(u$.prob, rep(4 / 96, 4))
})

test_that("a list read back from a file is rebuilt whole", {
  df <- make_df()
  des <- design_simple(n = 6, reserve = 3)
  s <- draw(df, des, seed = 5, weights = TRUE)
  back <- as.data.frame(s)
  attributes(back) <- attributes(back)[c("names", "row.names", "class")]
  back <- back[c(9, 1:8), ]
  rebuilt <- as_drawn(back, des, df)
  expect_equal(attr(rebuilt, "drawn_rows"), s$.row)
  # The reserve order comes from `.reserve`, not from the rows' order.
  u <- use_reserves(rebuilt, rebuilt$id == s$id[1])
  expect_equal(u$id, s$id[c(2:6, 7)])
  expect_equal(ht_total(u, "value")$se,
               ht_total(use_reserves(s, s$id == s$id[1]), "value")$se)
  expect_error(as_drawn(back[1:6, ], des, df), "Rebuild the whole list")
  back$.reserve <- NULL
  expect_error(as_drawn(back, des, df), "no `.reserve` column")
})

test_that("the first n available units of a list are a simple random sample", {
  # Rows 1 to 10 cannot be reached. Each of the other 30 is used with
  # probability 8 / 30, whichever units were visited.
  pop <- data.frame(id = 1:40)
  hits <- integer(40)
  for (i in 1:1000) {
    s <- draw(pop, design_simple(n = 8, reserve = 12), seed = i,
              weights = TRUE)
    u <- use_reserves(s, s$id <= 10)
    hits[u$.row] <- hits[u$.row] + 1L
  }
  p <- hits / 1000
  expect_equal(p[1:10], rep(0, 10))
  z <- (p[11:40] - 8 / 30) / sqrt(8 / 30 * 22 / 30 / 1000)
  expect_lt(max(abs(z)), 4.5)
})
