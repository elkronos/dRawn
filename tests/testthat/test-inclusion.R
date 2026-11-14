test_that("closed forms match the textbook values", {
  d1 <- data.frame(id = 1:10)
  expect_equal(inclusion_prob(d1, design_simple(n = 3)), rep(0.3, 10))

  d2 <- data.frame(id = 1:20, g = rep(c("a", "b"), times = c(15, 5)))
  # n = 8 proportional over 15/5 allocates 6 and 2
  expect_equal(inclusion_prob(d2, design_stratified("g", n = 8)),
               c(rep(6/15, 15), rep(2/5, 5)))

  d4 <- data.frame(id = 1:20, cl = rep(letters[1:5], each = 4))
  expect_equal(inclusion_prob(d4, design_cluster("cl", n_clusters = 2)),
               rep(0.4, 20))

  expect_equal(inclusion_prob(data.frame(id = 1:12), design_reservoir(n = 4)),
               rep(1/3, 12))
})

test_that("systematic inclusion is 1/interval for every row, at any N", {
  for (N in c(20, 23, 29)) {
    p <- inclusion_prob(data.frame(id = seq_len(N)),
                        design_systematic(interval = 7))
    expect_equal(p, rep(1/7, N), info = paste("N =", N))
  }
})

test_that("a fixed systematic start is refused as a probability sample", {
  expect_error(
    inclusion_prob(data.frame(id = 1:20), design_systematic(interval = 5, start = 2)),
    "not a probability sample"
  )
})

test_that("systematic inclusion is n / N with `n`, at any N", {
  for (N in c(20, 23, 29)) {
    p <- inclusion_prob(data.frame(id = seq_len(N)), design_systematic(n = 7))
    expect_equal(p, rep(7 / N, N), info = paste("N =", N))
  }
})

test_that("a fixed systematic start is refused by every weighting path alike", {
  d <- data.frame(id = 1:20)
  des <- design_systematic(interval = 5, start = 2)
  want <- paste0("A systematic design with a fixed `start` is not a ",
                 "probability sample: each row is selected with probability ",
                 "0 or 1. Leave `start` as NULL for a random start.")
  expect_error(inclusion_prob(d, des), want, fixed = TRUE)
  expect_error(sampling_weight(d, des), want, fixed = TRUE)
  expect_error(joint_prob(d, des, rows = 1:3), want, fixed = TRUE)
  expect_error(draw(d, des, weights = TRUE), want, fixed = TRUE)
  # The design is deterministic, and the simulation reports its 0s and 1s.
  expect_identical(inclusion_prob(d, des, simulate = TRUE, R = 20, seed = 1),
                   as.numeric(d$id %in% c(2, 7, 12, 17)))
})

test_that("sampling_weight is the reciprocal, and 0 becomes NA not Inf", {
  d <- data.frame(id = 1:10)
  expect_equal(sampling_weight(d, design_simple(n = 5)), rep(2, 10))

  # A row outside the window has probability 0 and no finite weight.
  ts <- data.frame(id = 1:6,
                   t = seq(as.POSIXct("2020-01-01", tz = "UTC"), by = "hour",
                           length.out = 6))
  w <- sampling_weight(ts, design_temporal("t", "2020-01-01", "2020-01-01 03:00:00",
                                         interval = 1, per_interval = 1,
                                         unit = "hours"))
  expect_true(all(is.na(w[4:6])))
  expect_true(all(is.finite(w[1:3])))
})

test_that("rows the design can never reach get probability 0", {
  ts <- data.frame(id = 1:6,
                   t = seq(as.POSIXct("2020-01-01", tz = "UTC"), by = "hour",
                           length.out = 6))
  p <- inclusion_prob(ts, design_temporal("t", "2020-01-01",
                                          "2020-01-01 03:00:00", interval = 1,
                                          per_interval = 1, unit = "hours"))
  expect_equal(p, c(1, 1, 1, 0, 0, 0))
})

test_that("designs with no closed form refuse rather than invent one", {
  df <- data.frame(id = 1:20, cl = rep(letters[1:4], each = 5), w = 1:20)

  expect_error(inclusion_prob(df, design_weighted("w", n = 3)),
               "no closed-form inclusion probability")
  # Balanced clusters have one: with clusters of equal size every selected
  # cluster is taken whole.
  expect_equal(
    inclusion_prob(df, design_cluster("cl", n_clusters = 2, balanced = TRUE)),
    rep(0.5, 20)
  )
  expect_error(
    inclusion_prob(df, design_multistage("cl", n_clusters = 2, n = 4,
                                         allocation = "proportional")),
    "no closed-form inclusion probability"
  )
  expect_error(inclusion_prob(df, design_bootstrap(n_replicates = 5)),
               "not a probability sample of a finite population")
})

test_that("simulate = TRUE works where the closed form refuses", {
  df <- data.frame(id = 1:12, cl = rep(letters[1:4], each = 3))
  p <- inclusion_prob(df, design_cluster("cl", n_clusters = 2, balanced = TRUE),
                      simulate = TRUE, R = 400, seed = 1)
  expect_length(p, 12)
  expect_true(all(p >= 0 & p <= 1))
  # Whole clusters of equal size: every row shares its cluster's rate.
  expect_equal(mean(p), 0.5, tolerance = 0.1)
})

test_that("simulation agrees with the closed form where both exist", {
  d <- data.frame(id = 1:20, g = rep(c("a", "b"), times = c(15, 5)))
  design <- design_stratified("g", n = 8)
  exact <- inclusion_prob(d, design)
  sim <- inclusion_prob(d, design, simulate = TRUE, R = 3000, seed = 1)
  expect_equal(sim, exact, tolerance = 0.06)
})

test_that("inclusion_prob validates its arguments", {
  df <- data.frame(id = 1:10)
  expect_error(inclusion_prob(df, "simple"), "must come from one of the design")
  expect_error(inclusion_prob(df, design_simple(n = 2), simulate = "yes"),
               "TRUE or FALSE")
  expect_error(inclusion_prob(df, design_simple(n = 2), simulate = TRUE, R = 0),
               "at least 1")
})

# ---- probabilities must describe the draw that actually happens ------------

test_that("multistage refuses where the per-cluster take is not constant", {
  # allocate() runs at draw time over the SELECTED clusters, so a remainder
  # goes to the largest of them and an undersized cluster is capped. Averaging
  # over the possibilities was out by 13% and could exceed 1.
  uneven <- data.frame(cl = rep(c("A", "B", "C", "D"), c(20, 10, 8, 6)))
  expect_error(inclusion_prob(uneven, design_multistage("cl", n_clusters = 2,
                                                        n = 7)),
               "no closed-form")
  tiny <- data.frame(cl = rep(c("x", "b", "c", "e"), c(1, 10, 10, 10)))
  expect_error(inclusion_prob(tiny, design_multistage("cl", n_clusters = 2,
                                                      n = 8)),
               "fewer than")

  # Where the take really is constant, it is exact.
  even <- data.frame(cl = rep(paste0("c", 1:8), each = 10))
  des <- design_multistage("cl", n_clusters = 4, n = 20)
  p <- inclusion_prob(even, des)
  expect_equal(unique(p), 0.25)
  sim <- inclusion_prob(even, des, simulate = TRUE, R = 4000, seed = 1)
  expect_lt(max(abs(sim - p)), 0.03)
})

test_that("max_items bounds what a reservoir design can reach", {
  d <- data.frame(y = 1:20)
  des <- design_reservoir(n = 4, max_items = 10)
  p <- inclusion_prob(d, des)
  expect_equal(p, c(rep(0.4, 10), rep(0, 10)))
  expect_equal(sum(p), 4)

  m <- joint_prob(d, des, rows = c(1, 2, 15))
  expect_equal(diag(m), c(0.4, 0.4, 0))
  expect_equal(m[1, 3], 0)               # can never co-occur
  expect_equal(m[1, 2], 4 * 3 / (10 * 9))
})

test_that("systematic probabilities follow `order_by`, not frame order", {
  # draw_design() sorts first. Reading residue classes off the unsorted frame
  # reported 0 for pairs that always co-occur and 1/k for pairs that never can.
  pop <- data.frame(id = 1:6, v = c(2, 5, 1, 6, 3, 4))
  des <- design_systematic(interval = 3, order_by = "v")
  m <- joint_prob(pop, des)
  sim <- joint_prob(pop, des, simulate = TRUE, R = 6000, seed = 1)
  expect_lt(max(abs(m - sim)), 0.04)
  expect_equal(m[1, 4], 0)

  # And a dropped row cannot be selected, so its probability is 0.
  na_pop <- data.frame(v = c(1:12, rep(NA, 8)))
  na_des <- design_systematic(interval = 3, order_by = "v", na_rm = TRUE)
  p <- inclusion_prob(na_pop, na_des)
  expect_equal(p, c(rep(1 / 3, 12), rep(0, 8)))
  expect_equal(sum(p), 4)
})

test_that("weights the draw would reject are rejected here too", {
  # pps_pi()'s rescaling loop runs off the rails on a non-positive weight and
  # returns a vector of 0s and 1s summing to the wrong total.
  bad <- data.frame(w = c(1, 2, 3, -4, 5))
  for (m in c("systematic", "poisson")) {
    expect_error(inclusion_prob(bad, design_weighted("w", n = 2, method = m)),
                 "must be positive")
  }
  expect_error(inclusion_prob(data.frame(w = c(1, 2, Inf)),
                              design_weighted("w", n = 2,
                                              method = "systematic")),
               "non-finite")
})

test_that("na_rm means the same thing on the probability path as on the draw", {
  cases <- list(
    list(d = data.frame(g = c(rep("a", 4), rep("b", 4), NA, NA)),
         des = function(na) design_stratified("g", n = 4, na_rm = na)),
    list(d = data.frame(g = c(rep("a", 4), rep("b", 4), NA, NA)),
         des = function(na) design_cluster("g", n_clusters = 1, na_rm = na)),
    list(d = data.frame(w = c(1:8, NA, NA)),
         des = function(na) design_weighted("w", n = 3, method = "poisson",
                                            na_rm = na)),
    list(d = data.frame(v = c(1:8, NA, NA)),
         des = function(na) design_systematic(interval = 3, order_by = "v",
                                              na_rm = na))
  )
  for (case in cases) {
    expect_error(draw(case$d, case$des(FALSE), seed = 1), "na_rm = TRUE")
    expect_error(inclusion_prob(case$d, case$des(FALSE)), "na_rm = TRUE")
    expect_equal(length(inclusion_prob(case$d, case$des(TRUE))), nrow(case$d))
  }
})

test_that("Neyman allocation reports probabilities for a design it can draw", {
  # stratum_spread() checked the whole frame while the draw checked only the
  # rows it kept, so a drawable design had no weights.
  d <- data.frame(g = c(rep("a", 5), rep("b", 5), NA, NA),
                  v = c(seq(1, 5), seq(2, 10, length.out = 5), NA, NA))
  des <- design_stratified("g", n = 4, allocation = "neyman",
                           allocation_by = "v", na_rm = TRUE)
  expect_equal(nrow(draw(d, des, seed = 1)), 4)
  p <- inclusion_prob(d, des)
  expect_equal(length(p), nrow(d))
  expect_equal(p[11:12], c(0, 0))
  expect_equal(nrow(draw(d, des, seed = 1, weights = TRUE)), 4)
})

test_that("joint_prob refuses a matrix past about 1 GB and says how to avoid it", {
  # 8 bytes a pair: every row of a 12,000-row frame is 1.2 GB, and a frame
  # some times larger could end the session.
  pop <- data.frame(id = seq_len(12000))
  d <- design_simple(n = 10)
  expect_error(joint_prob(pop, d),
               "12,000 x 12,000 matrix for every row of `data`, which needs about 1.2 GB",
               fixed = TRUE)
  expect_error(joint_prob(pop, d), "rows = s$.row", fixed = TRUE)
  expect_error(joint_prob(pop, d, rows = seq_len(12000)),
               "for the 12,000 `rows` given, which needs about 1.2 GB of memory. Pass fewer `rows`.",
               fixed = TRUE)
  expect_error(joint_prob(pop, d, simulate = TRUE, R = 2), "1.2 GB",
               fixed = TRUE)
  expect_equal(dim(joint_prob(pop, d, rows = c(5, 9000))), c(2L, 2L))
})

test_that("joints filled block by block equal the pairwise forms", {
  # Same group: a simple random sample of n from N. Different groups:
  # independent. A row in no group co-occurs with nothing. Built here pair by
  # pair.
  pairwise <- function(group, n_g, size_g, pi_i) {
    k <- length(group)
    out <- outer(pi_i, pi_i)
    for (i in seq_len(k)) for (j in seq_len(k)) {
      if (is.na(group[i]) || is.na(group[j])) {
        out[i, j] <- 0
      } else if (group[i] == group[j]) {
        n <- n_g[group[i]]
        N <- size_g[group[i]]
        out[i, j] <- if (N > 1) n * (n - 1) / (N * (N - 1)) else 0
      }
    }
    diag(out) <- pi_i
    out
  }
  group <- c(2L, 1L, NA, 3L, 2L, 1L, 2L, 4L, 1L)
  n_g <- c(2, 3, 1, 1)
  size_g <- c(5, 3, 4, 1)
  pi_i <- ifelse(is.na(group), 0, (n_g / size_g)[group])
  expect_equal(drawn:::joint_from_groups(group, n_g, size_g, pi_i),
               pairwise(group, n_g, size_g, pi_i))

  # Cluster: the same cluster a / A, different clusters a(a-1) / (A(A-1)).
  pop <- data.frame(cl = c("a", "a", "b", NA, "c", "b", "c", "c", "d"))
  d <- design_cluster("cl", n_clusters = 2, na_rm = TRUE)
  same <- outer(pop$cl, pop$cl, "==")
  want <- ifelse(same, 2 / 4, 2 * 1 / (4 * 3))
  want[is.na(pop$cl), ] <- 0
  want[, is.na(pop$cl)] <- 0
  diag(want) <- inclusion_prob(pop, d)
  expect_equal(joint_prob(pop, d), want)

  # Multistage, a = 2 of A = 4 clusters and m = 2 rows from each: the stages
  # multiplied.
  pop <- data.frame(cl = c(rep(c("a", "b", "c", "d"), times = c(3, 4, 5, 4)),
                           NA))
  d <- design_multistage("cl", n_clusters = 2, n = 4, na_rm = TRUE)
  N_i <- as.numeric(table(pop$cl)[pop$cl])
  same <- outer(pop$cl, pop$cl, "==")
  want <- ifelse(same, (2 / 4) * 2 / (N_i * (N_i - 1)),
                 2 / (4 * 3) * outer(2 / N_i, 2 / N_i))
  want[is.na(pop$cl), ] <- 0
  want[, is.na(pop$cl)] <- 0
  diag(want) <- inclusion_prob(pop, d)
  expect_equal(joint_prob(pop, d), want)
})

test_that("the probability paths refuse what draw() refuses, with its message", {
  # They answered instead: 1.25 for every row of a simple random sample of 30
  # from 24, a census of 1s for a stratified one, 0 for every row when no
  # interval fits the window.
  d <- data.frame(id = 1:24, site = rep(c("a", "b", "c"), c(12, 8, 4)),
                  cl = rep(sprintf("c%02d", 1:6), each = 4), w = c(1:23, 100),
                  ts = seq(as.POSIXct("2024-01-01", tz = "UTC"), by = "hour",
                           length.out = 24),
                  stringsAsFactors = FALSE)
  none <- d
  none$site <- NA_character_
  none$cl <- NA_character_
  none$w <- NA_real_
  none$ts[] <- NA
  cases <- list(
    list(d, design_simple(n = 30), "cannot exceed the number of rows (24)"),
    list(d, design_stratified("site", n = 30),
         "cannot exceed the number of rows (24)"),
    list(none, design_stratified("site", n = 3, na_rm = TRUE),
         "Every row has a missing stratum key."),
    list(d, design_cluster("cl", n_clusters = 7),
         "`n_clusters` (7) exceeds the 6 cluster(s) available."),
    list(none, design_cluster("cl", n_clusters = 1, na_rm = TRUE),
         "Every row has a missing cluster label."),
    list(d, design_multistage("cl", n_clusters = 7, n = 14),
         "`n_clusters` (7) exceeds the 6 cluster(s) available."),
    list(d, design_multistage("cl", n_clusters = 2, n = 10),
         "`n` (10) exceeds the 8 row(s) available across the 2 selected"),
    list(d, design_weighted("w", n = 30, method = "systematic"),
         "`n` (30) cannot exceed the number of rows (24) for method"),
    list(d, design_weighted("w", n = 30, method = "poisson"),
         "`n` (30) cannot exceed the number of rows (24) for method"),
    list(none, design_weighted("w", n = 3, method = "poisson", na_rm = TRUE),
         "Every row has a missing weight in `w`."),
    list(d, design_certainty("w", 20, design_simple(n = 30)),
         "applied to the 19 rows with `w` below 20: `n` (30) cannot exceed"),
    list(d, design_certainty("w", 20, design_weighted("w", n = 30,
                                                      method = "poisson")),
         "applied to the 19 rows with `w` below 20: `n` (30) cannot exceed"),
    list(d, design_temporal("ts", "2024-01-01", "2024-01-01 02:00:00",
                            interval = 3, per_interval = 1),
         "shorter than one interval"),
    list(none, design_temporal("ts", "2024-01-01", "2024-01-02", interval = 3,
                               per_interval = 1, na_rm = TRUE),
         "Every row has a missing or unparseable timestamp in `ts`."),
    list(d, design_systematic(interval = 3, order_by = "W"),
         "`data` is missing column: `W`."),
    list(none, design_systematic(interval = 3, order_by = "w", na_rm = TRUE),
         "Every row has a missing `w`.")
  )
  if (requireNamespace("sf", quietly = TRUE)) {
    d$lon <- seq(-10, 10, length.out = 24)
    d$lat <- seq(-5, 5, length.out = 24)
    none$lon <- NA_real_
    none$lat <- d$lat
    gap <- d
    gap$lat[3] <- NA
    box <- function(lo, hi) {
      sf::st_sfc(sf::st_polygon(list(cbind(c(lo, hi, hi, lo, lo),
                                           c(lo, lo, hi, hi, lo)))),
                 crs = 4326)
    }
    cases <- c(cases, list(
      list(d, design_spatial(c("lon", "lat"), box(-2, 2), n = 10),
           "`n` (10) exceeds the 4 row(s) inside `region`."),
      list(d, design_spatial(c("lon", "lat"), box(50, 60), n = 1),
           "No rows fall inside `region`."),
      list(none, design_spatial(c("lon", "lat"), box(-2, 2), n = 1,
                                na_rm = TRUE),
           "Every row has missing coordinates in `lon` or `lat`."),
      list(gap, design_spatial(c("lon", "lat"), box(-2, 2), n = 1),
           "1 row(s) have missing coordinates in `lon` or `lat`")
    ))
  }
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  # The 10-degree box at latitude 50 to 60 bows past the edge warning's
  # limit; only the errors are compared here.
  message_of <- function(code) {
    tryCatch({
      withCallingHandlers(force(code), drawn_geodesic_edges = function(w) {
        invokeRestart("muffleWarning")
      })
      NA_character_
    }, error = conditionMessage)
  }
  for (case in cases) {
    pop <- case[[1]]
    des <- case[[2]]
    info <- drawn:::design_label(des)
    want <- message_of(draw(pop, des, seed = 1))
    expect_true(grepl(case[[3]], want, fixed = TRUE), info = info)
    expect_identical(message_of(inclusion_prob(pop, des)), want, info = info)
    expect_identical(message_of(sampling_weight(pop, des)), want, info = info)
    expect_identical(message_of(joint_prob(pop, des, rows = 1:2)), want,
                     info = info)
    expect_identical(message_of(draw(pop, des, seed = 1, weights = TRUE)),
                     want, info = info)
    expect_true(grepl(want, message_of(plot(des, pop, type = "probability")),
                      fixed = TRUE), info = info)
  }
})

test_that("joint_prob pairs a repeated row with itself at its own probability", {
  # A row given twice was two units: the pair it made with itself got
  # n(n-1)/(N(N-1)), and simulate = TRUE gave it 0.
  d <- data.frame(id = 1:20, g = rep(c("a", "b"), each = 10))
  for (des in list(design_simple(n = 5), design_stratified("g", n = 6))) {
    p <- inclusion_prob(d, des)
    m <- joint_prob(d, des, rows = c(3, 3, 12))
    expect_equal(dim(m), c(3L, 3L))
    expect_equal(m[1, 2], p[3])
    expect_equal(m[1:2, 1:2], matrix(p[3], 2, 2))
    expect_equal(m[c(1, 3), c(1, 3)], joint_prob(d, des, rows = c(3, 12)))
    expect_equal(m[2, 3], m[1, 3])
    sim <- joint_prob(d, des, rows = c(3, 3, 12), simulate = TRUE, R = 200,
                      seed = 1)
    expect_equal(sim[1, 2], sim[1, 1])
    expect_equal(sim[2, 2], sim[1, 1])
    expect_gt(sim[1, 1], 0)
  }
})

test_that("joint_prob refuses rows that are not whole numbers", {
  d <- data.frame(id = 1:20)
  expect_error(joint_prob(d, design_simple(n = 5), rows = c(1, 1.5)),
               "`rows` must be whole row numbers, not 1.5.", fixed = TRUE)
  expect_error(joint_prob(d, design_simple(n = 5), rows = 1.7,
                          simulate = TRUE, R = 10),
               "not 1.7", fixed = TRUE)
  # Past five values the list counts the rest, so it ends in one full stop.
  expect_error(joint_prob(d, design_simple(n = 5), rows = 1:8 + 0.5),
               "not 1.5, 2.5, 3.5, 4.5, 5.5 and 3 more. Pass", fixed = TRUE)
})

# ---- where a refusal points ---------------------------------------------------

test_that("a refusal names the function that simulates", {
  df <- data.frame(id = 1:20, cl = rep(letters[1:4], each = 5), w = 1:20)
  succ <- design_weighted("w", n = 3)
  # draw() and plot() take no `simulate`, so the pointer names the function
  # that does.
  for (call in list(quote(inclusion_prob(df, succ)),
                    quote(draw(df, succ, seed = 1, weights = TRUE)))) {
    err <- tryCatch(eval(call), error = conditionMessage)
    expect_match(err, "inclusion_prob(data, design, simulate = TRUE)",
                 fixed = TRUE)
    expect_no_match(err, "Or pass simulate", fixed = TRUE)
  }
  # Joint refusals point to a simulation where one works ...
  for (des in list(succ,
                   design_cluster("cl", n_clusters = 2, balanced = TRUE),
                   design_multistage("cl", n_clusters = 2, n = 4,
                                     allocation = "proportional"),
                   design_simple(n = 3, replace = TRUE))) {
    expect_error(joint_prob(df, des, rows = 1:3),
                 "joint_prob(data, design, simulate = TRUE)", fixed = TRUE)
  }
  # ... and not for a bootstrap, where it is refused.
  err <- tryCatch(joint_prob(df, design_bootstrap(n_replicates = 2),
                             rows = 1:3), error = conditionMessage)
  expect_no_match(err, "simulate = TRUE", fixed = TRUE)
})

# ---- warnings and time in a simulation --------------------------------------

test_that("a warning each draw raises is given once, with the count", {
  short <- data.frame(y = 1:5)
  caught <- list()
  p <- withCallingHandlers(
    inclusion_prob(short, design_systematic(interval = 8), simulate = TRUE,
                   R = 400, seed = 1),
    warning = function(w) {
      caught[[length(caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  expect_length(caught, 1L)
  expect_s3_class(caught[[1L]], "drawn_empty_start")
  expect_match(conditionMessage(caught[[1L]]),
               "^The random start \\([678]\\) is past the last row \\(5\\)")
  # The start shown is the first draw's, and the count covers every draw
  # with a start past the last row, whichever it was. Each other draw takes
  # one row, so the empty draws are the 400 less the rows drawn.
  empty <- 400L - as.integer(round(sum(p) * 400))
  expect_match(conditionMessage(caught[[1L]]),
               paste0("This is the first of these warnings. ", empty,
                      " of the 400 simulated draws gave one like it, with ",
                      "details that differ from draw to draw."),
               fixed = TRUE)
  expect_no_match(conditionMessage(caught[[1L]]), "This happened in",
                  fixed = TRUE)
  expect_lt(max(abs(p - 1 / 8)), 0.05)

  # A message that is the same in every draw is given with its count.
  expect_warning(
    inclusion_prob(short, design_systematic(interval = 8, start = 7),
                   simulate = TRUE, R = 30, seed = 1),
    paste("`start` (7) is past the last row (5); returning no rows.",
          "This happened in 30 of the 30 simulated draws."),
    fixed = TRUE, class = "drawn_empty_start")

  # simulate_joint() counts them the same way.
  caught <- list()
  withCallingHandlers(
    joint_prob(short, design_systematic(interval = 8), simulate = TRUE,
               R = 50, seed = 1),
    warning = function(w) {
      caught[[length(caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  expect_length(caught, 1L)
  expect_s3_class(caught[[1L]], "drawn_empty_start")

  # A fixed start keeps its wording.
  expect_warning(draw(short, design_systematic(interval = 8, start = 7)),
                 "`start` (7) is past the last row (5); returning no rows.",
                 fixed = TRUE, class = "drawn_empty_start")
  expect_warning(draw(short[1, , drop = FALSE],
                      design_systematic(interval = 50, n_starts = 2),
                      seed = 1),
                 "Every random start", class = "drawn_empty_start")
})

test_that("a long simulation says how long it will take", {
  df <- data.frame(id = 1:20)
  des <- design_simple(n = 3)
  # Each reading of the clock is 10 seconds after the last, so the first draw
  # takes 10 seconds and 40 draws project 400.
  clock <- 0
  local_mocked_bindings(elapsed_time = function() {
    clock <<- clock + 10
    clock
  })
  expect_message(p <- inclusion_prob(df, des, simulate = TRUE, R = 40,
                                     seed = 1),
                 "Simulating 40 draws will take about 6.7 minutes, at 10 seconds a draw.",
                 fixed = TRUE)
  # Timing and messages leave the draws as they were.
  expect_identical(p, suppressMessages(
    inclusion_prob(df, des, simulate = TRUE, R = 40, seed = 1)))
  expect_message(joint_prob(df, des, rows = 1:3, simulate = TRUE, R = 40,
                            seed = 1), "Simulating 40 draws")

  # An interactive session also gets a progress bar.
  local_mocked_bindings(is_interactive = function() TRUE)
  shown <- capture.output(suppressMessages(
    inclusion_prob(df, des, simulate = TRUE, R = 40, seed = 1)))
  expect_match(paste(shown, collapse = ""), "100%", fixed = TRUE)
})

test_that("a short simulation is quiet", {
  df <- data.frame(id = 1:20)
  local_mocked_bindings(is_interactive = function() TRUE)
  expect_silent(inclusion_prob(df, design_simple(n = 3), simulate = TRUE,
                               R = 200, seed = 1))
})

test_that("durations read in the unit that suits them", {
  expect_identical(format_duration(0.0123), "0.012 seconds")
  expect_identical(format_duration(1), "1 second")
  expect_identical(format_duration(400), "6.7 minutes")
  expect_identical(format_duration(3600 * 3.3), "3.3 hours")
})
