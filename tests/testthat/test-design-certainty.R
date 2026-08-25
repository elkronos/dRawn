cert_pop <- function(N = 400) {
  set.seed(909)
  d <- data.frame(
    id    = seq_len(N),
    site  = rep(c("a", "b"), times = c(N * 0.6, N * 0.4)),
    value = round(stats::rlnorm(N, 6, 1.2))
  )
  d
}

test_that("every row at or above the threshold is taken, every time", {
  d <- cert_pop()
  thr <- stats::quantile(d$value, 0.9)
  big <- d$id[d$value >= thr]
  des <- design_certainty("value", thr, design_simple(n = 30))

  for (i in 1:20) {
    s <- draw(d, des, seed = i)
    expect_true(all(big %in% s$id))
    expect_equal(nrow(s), length(big) + 30)
  }
})

test_that("the sample keeps the frame's schema and order", {
  d <- cert_pop()
  s <- draw(d, design_certainty("value", 2000, design_simple(n = 20)), seed = 1)
  expect_same_schema(s, d)
  expect_false(is.unsorted(s$id))
})

test_that("inclusion probabilities are exactly 1 above and correct below", {
  d <- cert_pop()
  thr <- 3000
  des <- design_certainty("value", thr, design_simple(n = 40))
  p <- inclusion_prob(d, des)

  above <- d$value >= thr
  expect_true(all(p[above] == 1))
  expect_equal(unique(p[!above]), 40 / sum(!above))
  expect_equal(sum(p), sum(above) + 40)
})

test_that("simulated inclusion matches the closed form", {
  d <- cert_pop(80)
  des <- design_certainty("value", 1500, design_stratified("site", n = 20))
  exact <- inclusion_prob(d, des)

  set.seed(21)
  hits <- integer(nrow(d))
  R <- 3000
  for (i in seq_len(R)) {
    s <- draw(d, des)
    hits[s$id] <- hits[s$id] + 1L
  }
  expect_lt(max(abs(hits / R - exact)), 0.03)
})

test_that("the Horvitz-Thompson total is unbiased and its variance honest", {
  d <- cert_pop(300)
  truth <- sum(d$value)
  des <- design_certainty("value", 2500, design_stratified("site", n = 50))

  out <- vapply(1:400, function(i) {
    s <- draw(d, des, seed = i, weights = TRUE)
    e <- ht_total(s, "value")
    c(e$total, e$variance)
  }, numeric(2))

  se <- stats::sd(out[1, ]) / sqrt(ncol(out))
  expect_lt(abs(mean(out[1, ]) - truth), 4 * se)

  ratio <- mean(out[2, ]) / stats::var(out[1, ])
  expect_gt(ratio, 0.75)
  expect_lt(ratio, 1.35)
})

test_that("certainty rows add nothing to the variance", {
  d <- cert_pop(200)
  thr <- 3000
  des <- design_certainty("value", thr, design_simple(n = 40))
  s <- draw(d, des, seed = 4, weights = TRUE)

  # Doubling the certainty rows' y changes the estimate but not its variance.
  y1 <- s$value
  y2 <- s$value * ifelse(s$.prob == 1, 2, 1)
  expect_gt(sum(y2 / s$.prob), sum(y1 / s$.prob))
  expect_equal(ht_total(s, y2)$variance, ht_total(s, y1)$variance)
})

test_that("joint probabilities involving a certainty row collapse correctly", {
  d <- data.frame(id = 1:10, value = c(rep(1, 8), 100, 100))
  des <- design_certainty("value", 50, design_simple(n = 4))
  m <- joint_prob(d, des)
  p <- inclusion_prob(d, des)

  expect_equal(m[9, 10], 1)             # two certainty rows always co-occur
  expect_equal(m[9, 1], p[1])           # certainty with a sampled row
  expect_equal(m[1, 2], 4 * 3 / (8 * 7))  # both from the sampled part
  expect_equal(diag(m), p)
})

test_that("joint probabilities hold at both edges of the threshold", {
  d <- data.frame(id = 1:10, value = c(rep(1, 8), 100, 100))
  inner <- design_simple(n = 4)

  # Nothing certain: identical to the inner design
  expect_equal(joint_prob(d, design_certainty("value", Inf, inner)),
               joint_prob(d, inner))

  # Everything certain: every pair co-occurs with probability 1
  expect_true(all(joint_prob(d, design_certainty("value", 0, inner)) == 1))

  # Exactly one certain row, which is the case a zero-length branch could break
  d1 <- data.frame(id = 1:10, value = c(rep(1, 9), 100))
  m <- joint_prob(d1, design_certainty("value", 50, inner))
  expect_equal(m[10, 10], 1)
  expect_equal(m[1, 10], 4 / 9)          # certain row, so just the other's pi
  expect_equal(m[1, 2], 4 * 3 / (9 * 8)) # both from the nine below
  expect_equal(diag(m)[1], 4 / 9)
})

test_that("a threshold that takes everything still works", {
  d <- cert_pop(50)
  s <- draw(d, design_certainty("value", 0, design_simple(n = 10)), seed = 1)
  expect_equal(nrow(s), nrow(d))
  expect_equal(inclusion_prob(d, design_certainty("value", 0,
                                                  design_simple(n = 10))),
               rep(1, nrow(d)))
})

test_that("a threshold that takes nothing degrades to the inner design", {
  d <- cert_pop(50)
  inner <- design_simple(n = 10)
  s <- draw(d, design_certainty("value", Inf, inner), seed = 7)
  expect_equal(nrow(s), 10)
  expect_equal(inclusion_prob(d, design_certainty("value", Inf, inner)),
               inclusion_prob(d, inner))
})

test_that("missing size values error, or drop when asked", {
  d <- cert_pop(60)
  d$value[c(3, 40)] <- NA
  des <- design_certainty("value", 1000, design_simple(n = 10))
  expect_error(draw(d, des, seed = 1), "missing")

  ok <- draw(d, design_certainty("value", 1000, design_simple(n = 10),
                                 na_rm = TRUE), seed = 1)
  expect_false(any(is.na(ok$value)))
})

test_that("the constructor rejects what it cannot compose", {
  expect_error(design_certainty("value", 10, "simple"), "must be a design")
  nested <- design_certainty("value", 10, design_simple(n = 5))
  expect_error(design_certainty("value", 20, nested),
               "cannot itself be a certainty design")
  expect_error(design_certainty("value", "big", design_simple(n = 5)),
               "single number")
  expect_error(design_certainty(c("a", "b"), 10, design_simple(n = 5)))
})

test_that("a non-numeric size column is caught with the column named", {
  d <- cert_pop(40)
  des <- design_certainty("site", 1, design_simple(n = 5))
  expect_error(draw(d, des, seed = 1), "site")
})

test_that("the design prints its parts", {
  des <- design_certainty("value", 20000, design_stratified("site", n = 40))
  out <- utils::capture.output(print(des))
  expect_true(any(grepl("certainty", out)))
  expect_true(any(grepl("20000|20,000", out)))
})

# ---- composition, which is where the bugs live -----------------------------

test_that("the variance defers to `rest` instead of assuming a fixed size", {
  # A Poisson `rest` has a random sample size and independent units. Handling
  # the certainty design as one flat fixed-size design instead silently applies
  # Sen-Yates-Grundy to pi_ij = pi_i * pi_j, where every term cancels and the
  # standard error comes out as exactly zero.
  d <- cert_pop(400)
  des <- design_certainty("value", 2500,
                          design_weighted("value", n = 50, method = "poisson"))

  out <- vapply(1:400, function(i) {
    e <- ht_total(draw(d, des, seed = i, weights = TRUE), "value")
    c(e$total, e$variance)
  }, numeric(2))

  expect_true(all(out[2, ] > 0))
  ratio <- mean(out[2, ]) / stats::var(out[1, ])
  expect_gt(ratio, 0.75)
  expect_lt(ratio, 1.35)
})

test_that("a systematic `rest` refuses a variance rather than returning one", {
  d <- cert_pop(300)
  des <- design_certainty("value", 2500, design_systematic(interval = 8))
  s <- draw(d, des, seed = 1, weights = TRUE)

  # Analytic is impossible here; taken flat it would return a negative number.
  strict <- ht_total(s, "value", variance = "analytic")
  expect_true(is.na(strict$variance))
  expect_match(strict$note, "no design-unbiased variance")

  # The jackknife cannot rescue it either: a systematic `rest` has a single
  # primary sampling unit. Saying so beats reporting a number from a method
  # that declined.
  auto <- ht_total(s, "value")
  expect_true(is.na(auto$variance))
  expect_equal(auto$method, "none")
  expect_match(auto$note, "one primary sampling unit")
  expect_match(auto$note, "no design-unbiased variance")
  expect_true(is.na(deff(auto)))
})

test_that("every sampled row certain means an exact total, variance 0", {
  d <- cert_pop(50)
  s <- draw(d, design_certainty("value", 0, design_simple(n = 10)), seed = 1,
            weights = TRUE)
  r <- ht_total(s, "value")
  expect_equal(r$total, sum(d$value))
  expect_equal(r$variance, 0)
  expect_match(r$note, "exact rather than estimated")
})

test_that("dropping missing sizes does not shift the probabilities", {
  # certainty_split() removes rows, so positions in the reduced frame are not
  # positions in the original one. Getting that wrong shifts every probability
  # by the number of rows dropped above it.
  d <- data.frame(id = 1:10, value = c(NA, rep(1, 7), 100, 100))
  des <- design_certainty("value", 50, design_simple(n = 3), na_rm = TRUE)

  p <- inclusion_prob(d, des)
  expect_equal(p, c(0, rep(3 / 7, 7), 1, 1))

  m <- joint_prob(d, des, rows = c(2, 3, 9))
  expect_equal(diag(m), c(3 / 7, 3 / 7, 1))
  expect_equal(m[1, 2], 3 * 2 / (7 * 6))
  expect_equal(m[1, 3], 3 / 7)

  # A dropped row is in no sample, so every joint probability involving it is
  # 0 -- which is the answer. `rows` defaults to the whole frame, so erroring
  # here would make the documented default unusable.
  full <- joint_prob(d, des)
  expect_equal(dim(full), c(10L, 10L))
  expect_true(all(full[1, ] == 0))
  expect_true(all(full[, 1] == 0))
  expect_equal(full[2, 3], 3 * 2 / (7 * 6))

  # And a draw never selects the dropped row
  for (i in 1:20) expect_false(1 %in% draw(d, des, seed = i)$id)
})
