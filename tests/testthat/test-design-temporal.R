test_that("returns only rows inside the window", {
  df <- make_df()
  d <- design_temporal("time", from = "2020-01-01 00:00:00",
                       to = "2020-01-01 18:00:00", interval = 6,
                       per_interval = 2, unit = "hours")
  res <- draw(df, d, seed = 123)
  tv <- as.POSIXct(res$time, tz = "UTC")
  expect_true(all(tv >= as.POSIXct("2020-01-01 00:00:00", tz = "UTC") &
                    tv < as.POSIXct("2020-01-01 18:00:00", tz = "UTC")))
  expect_same_schema(res, df)
})

test_that("each interval contributes at most per_interval rows", {
  df <- make_ts()
  d <- design_temporal("ts", "2020-01-01", "2020-01-03", interval = 6,
                       per_interval = 2, unit = "hours")
  res <- draw(df, d, seed = 1)
  bucket <- as.integer(difftime(res$ts, as.POSIXct("2020-01-01", tz = "UTC"),
                                units = "hours")) %/% 6
  expect_true(all(table(bucket) <= 2))
  expect_equal(nrow(res), 16)
})

test_that("POSIXct and Date columns are used as they are", {
  expect_silent(
    res <- draw(make_ts(),
                design_temporal("ts", "2020-01-01", "2020-01-02", interval = 6,
                                per_interval = 2, unit = "hours"), seed = 1)
  )
  expect_equal(nrow(res), 8)

  dates <- data.frame(id = 1:10,
                      d = seq(as.Date("2020-01-01"), by = "day", length.out = 10))
  res_d <- draw(dates, design_temporal("d", "2020-01-01", "2020-01-06",
                                       interval = 1, per_interval = 1,
                                       unit = "days"), seed = 1)
  expect_equal(nrow(res_d), 5)
})

test_that("calendar units step by calendar boundaries", {
  df <- data.frame(
    id = 1:365,
    ts = seq(as.POSIXct("2020-01-01", tz = "UTC"), by = "day", length.out = 365)
  )
  res <- draw(df, design_temporal("ts", "2020-01-01", "2020-04-01", interval = 1,
                                  per_interval = 1, unit = "months"), seed = 1)
  expect_equal(nrow(res), 3)
  expect_equal(sort(as.integer(format(res$ts, "%m"))), 1:3)
})

test_that("tz shifts the window rather than relabelling the output", {
  df <- make_ts()
  utc <- draw(df, design_temporal("ts", "2020-01-01 00:00:00",
                                  "2020-01-01 06:00:00", interval = 6,
                                  per_interval = 100, unit = "hours"), seed = 1)
  # The same wall-clock bounds in a zone 5 hours behind cover different instants.
  shifted <- draw(df, design_temporal("ts", "2020-01-01 00:00:00",
                                      "2020-01-01 06:00:00", interval = 6,
                                      per_interval = 100, unit = "hours",
                                      tz = "America/New_York"), seed = 1)
  expect_false(identical(utc$id, shifted$id))
})

test_that("a window shorter than one interval is reported clearly", {
  df <- make_df()
  expect_error(
    draw(df, design_temporal("time", "2020-01-01 00:00:00",
                             "2020-01-01 03:00:00", interval = 6,
                             per_interval = 2, unit = "hours")),
    "shorter than one interval"
  )
})

test_that("an empty window keeps the schema", {
  df <- make_df()
  res <- draw(df, design_temporal("time", "2021-01-01", "2021-01-02",
                                  interval = 6, per_interval = 2,
                                  unit = "hours"), seed = 1)
  expect_equal(nrow(res), 0L)
  expect_same_schema(res, df)
})

test_that("bad bounds and units are rejected", {
  df <- make_df()
  expect_error(design_temporal("time", "2020-01-01", "2020-01-02", interval = 6,
                               per_interval = 2, unit = "fortnights"))
  expect_error(
    draw(df, design_temporal("time", "2020-01-02", "2020-01-01", interval = 6,
                             per_interval = 2, unit = "hours")),
    "must be after"
  )
})

test_that("unparseable timestamps follow the na_rm contract", {
  df <- make_df()
  df$time[1:3] <- "not a date"
  expect_error(
    draw(df, design_temporal("time", "2020-01-01", "2020-01-02", interval = 6,
                             per_interval = 2, unit = "hours"), seed = 1),
    "unparseable timestamp"
  )
  expect_silent(
    draw(df, design_temporal("time", "2020-01-01", "2020-01-02", interval = 6,
                             per_interval = 2, unit = "hours", na_rm = TRUE),
         seed = 1)
  )
})

test_that("a backwards window is refused on every path, not just the draw", {
  d <- data.frame(ts = seq(as.POSIXct("2020-01-01", tz = "UTC"), by = "hour",
                           length.out = 10))
  expect_error(design_temporal("ts", from = "2020-01-02", to = "2020-01-01",
                               interval = 1, per_interval = 1, unit = "hours"),
               "`to` must be after `from`")
  # A design built before the constructor checked its bounds
  bad <- design_temporal("ts", from = "2020-01-01", to = "2020-01-02",
                         interval = 1, per_interval = 1, unit = "hours")
  bad$from <- "2020-01-02"
  bad$to <- "2020-01-01"
  expect_error(draw(d, bad), "`to` must be after `from`")
  expect_error(inclusion_prob(d, bad), "`to` must be after `from`")
  expect_error(joint_prob(d, bad), "`to` must be after `from`")
})

test_that("a column mixing dates with date-times keeps both", {
  # The date-only fallback fired only when nothing at all parsed, so a mixed
  # column silently lost every date-only row.
  parse_time <- drawn:::parse_time
  got <- parse_time(c("2020-01-01", "2020-01-01 05:00:00", "2020-01-02"),
                    "UTC", "x")
  expect_false(anyNA(got))
  expect_equal(format(got[1], "%Y-%m-%d %H:%M:%S"), "2020-01-01 00:00:00")

  d <- data.frame(ts = c("2024-03-01", "2024-03-01 06:00:00",
                         "2024-03-02 09:00:00", "2024-03-02"))
  des <- design_temporal("ts", from = "2024-03-01", to = "2024-03-03",
                         interval = 1, per_interval = 1, unit = "days")
  expect_equal(sum(inclusion_prob(d, des) > 0), 4)
})

sparse_days <- function() {
  # Nine days with three rows each, which per_interval = 3 takes whole, and
  # one busy day of 30 rows from which three are drawn.
  ts <- c(as.POSIXct("2024-01-01", tz = "UTC") + rep(0:8, each = 3) * 86400 +
            rep(c(1, 2, 3) * 3600, 9),
          as.POSIXct("2024-01-10", tz = "UTC") + (1:30) * 600)
  data.frame(ts = ts, v = c(rep(1, 27), 50, 60, rep(1, 28)),
             y = 100 + 20 * sin(1:57))
}

test_that("intervals taken whole carry no degrees of freedom", {
  # Counting the 27 rows of the sparse days gave 20 degrees of freedom where
  # the three random rows support 2.
  d <- sparse_days()
  des <- design_temporal("ts", from = "2024-01-01", to = "2024-01-11",
                         interval = 24, per_interval = 3, unit = "hours")
  s <- draw(d, des, seed = 1, weights = TRUE)
  expect_equal(sum(s$.prob == 1), 27L)
  expect_equal(ht_total(s, "y")$df, 2)

  # The same holds under a certainty design over the temporal one
  cert <- design_certainty("v", 50, des)
  sc <- draw(d, cert, seed = 1, weights = TRUE)
  expect_equal(sum(sc$v >= 50), 2L)
  expect_equal(ht_total(sc, "y")$df, 2)

  skip_if_not_installed("survey")
  expect_equal(survey::degf(as_svydesign(s)), 2)
  expect_equal(survey::degf(as_svydesign(sc)), 2)
  expect_equal(as.numeric(survey::SE(survey::svytotal(~y, as_svydesign(sc)))),
               ht_total(sc, "y")$se, tolerance = 1e-8)
})

day_frame <- function(from, days) {
  data.frame(day = as.Date(from) + seq_len(days) - 1L, y = seq_len(days))
}

test_that("days and weeks start at local midnight across a clock change", {
  # A fixed 86,400-second step drifted off midnight after the change, so two
  # calendar days shared one interval and the last day was lost.
  for (zone in c("America/New_York", "Europe/London")) {
    for (start in c("2024-03-07", "2024-10-24", "2024-10-31")) {
      d <- day_frame(start, 10)
      des <- design_temporal("day", from = start,
                             to = as.Date(start) + 10, interval = 1,
                             per_interval = 1, unit = "days", tz = zone)
      b <- drawn:::temporal_bucket(des, d)
      expect_identical(as.integer(b), 1:10)
      expect_identical(levels(b),
                       paste(format(as.Date(start) + 0:9), "00:00"))
      expect_equal(inclusion_prob(d, des), rep(1, 10))
      expect_equal(nrow(draw(d, des, seed = 1)), 10L)
    }
  }

  # Hourly rows: every interval opens on a local midnight.
  zone <- "America/New_York"
  h <- data.frame(ts = seq(as.POSIXct("2024-03-08", tz = zone),
                           as.POSIXct("2024-03-13", tz = zone), by = "hour"))
  des <- design_temporal("ts", from = "2024-03-08", to = "2024-03-13",
                         interval = 1, per_interval = 2, unit = "days")
  b <- drawn:::temporal_bucket(des, h)
  first <- h$ts[!duplicated(b) & !is.na(b)]
  expect_equal(format(first, "%H:%M"), rep("00:00", 5))
  expect_equal(as.vector(table(b)), c(24, 24, 23, 24, 24))

  # Weeks are seven calendar days each.
  d <- day_frame("2024-03-04", 28)
  des <- design_temporal("day", from = "2024-03-04", to = "2024-04-01",
                         interval = 1, per_interval = 7, unit = "weeks",
                         tz = "America/New_York")
  expect_equal(as.vector(table(drawn:::temporal_bucket(des, d))), rep(7, 4))
})

test_that("a Date on a day with no midnight is its first instant", {
  # Santiago moved its clocks from 00:00 to 01:00 on 2022-09-11. The Date
  # became 23:00 the day before, and the same day as a string did not parse.
  zone <- "America/Santiago"
  start <- drawn:::parse_time(as.Date("2022-09-11"), zone, "x")
  expect_equal(format(start, "%Y-%m-%d %H:%M %z"), "2022-09-11 01:00 -0300")
  expect_equal(drawn:::parse_time("2022-09-11", zone, "x"), start)

  d <- day_frame("2022-08-11", 61)
  months <- design_temporal("day", from = "2022-08-11", to = "2022-10-11",
                            interval = 1, per_interval = 1, unit = "months",
                            tz = zone)
  expect_equal(as.vector(table(drawn:::temporal_bucket(months, d))),
               c(31, 30))
  expect_equal(inclusion_prob(d, months)[d$day == as.Date("2022-09-11")],
               1 / 30)

  # A Date `to` on that day ends the window at the end of the day before.
  days <- design_temporal("day", from = as.Date("2022-09-08"),
                          to = as.Date("2022-09-11"), interval = 1,
                          per_interval = 1, unit = "days", tz = zone)
  d <- day_frame("2022-09-08", 5)
  expect_equal(inclusion_prob(d, days), c(1, 1, 1, 0, 0))

  # Rows of strings on that day parse under na_rm = FALSE.
  s <- data.frame(day = format(as.Date("2022-09-08") + 0:4), y = 1)
  from_strings <- design_temporal("day", from = "2022-09-08",
                                  to = "2022-09-13", interval = 1,
                                  per_interval = 1, unit = "days", tz = zone)
  expect_equal(inclusion_prob(s, from_strings), rep(1, 5))
})

test_that("frame order does not move a Date on a repeated midnight", {
  # The Azores repeated 00:00-01:00 on 2019-10-27, and which midnight a Date
  # became depended on the rows beside it.
  zone <- "Atlantic/Azores"
  d <- day_frame("2019-10-24", 7)
  des <- design_temporal("day", from = "2019-10-27", to = "2019-10-31",
                         interval = 1, per_interval = 1, unit = "days",
                         tz = zone)
  p <- inclusion_prob(d, des)
  expect_equal(p, c(0, 0, 0, 1, 1, 1, 1))
  expect_equal(rev(inclusion_prob(d[7:1, ], des)), p)
  expect_equal(
    format(drawn:::parse_time(as.Date(c("2019-10-28", "2019-10-27")), zone,
                              "x"), "%H:%M %z"),
    c("00:00 -0100", "00:00 +0000"))
})

test_that("a day the zone skipped does not stop the draw", {
  # Samoa went from 2011-12-29 to 2011-12-31. A Date column holding
  # 2011-12-30 stopped every path, even with na_rm = TRUE.
  d <- data.frame(day = seq(as.Date("2011-12-01"), as.Date("2012-01-31"),
                            by = "day"), y = 1)
  des <- design_temporal("day", from = "2011-12-01", to = "2012-02-01",
                         interval = 1, per_interval = 3, unit = "months",
                         tz = "Pacific/Apia", na_rm = TRUE)
  expect_equal(inclusion_prob(d, des), rep(3 / 31, 62))
  expect_equal(nrow(draw(d, des, seed = 1)), 6L)
  expect_equal(dim(joint_prob(d, des, rows = 1:3)), c(3L, 3L))
  expect_equal(drawn:::parse_time(as.Date("2011-12-30"), "Pacific/Apia", "x"),
               drawn:::parse_time(as.Date("2011-12-31"), "Pacific/Apia", "x"))
})

test_that("tz = NULL takes the zone of a POSIXct column", {
  # The window was read in UTC, five hours off for New York data.
  zone <- "America/New_York"
  ny <- data.frame(id = 1:72, ts = seq(as.POSIXct("2023-12-31", tz = zone),
                                       by = "hour", length.out = 72))
  des <- design_temporal("ts", from = "2024-01-01", to = "2024-01-02",
                         interval = 24, per_interval = 100, unit = "hours")
  s <- draw(ny, des, seed = 1)
  expect_equal(format(range(s$ts), "%Y-%m-%d %H:%M"),
               c("2024-01-01 00:00", "2024-01-01 23:00"))
  expect_equal(inclusion_prob(as.data.frame(list(ts = as.POSIXlt(ny$ts))), des),
               inclusion_prob(ny, des))

  # Naming UTC still reads the window in UTC.
  utc <- draw(ny, design_temporal("ts", from = "2024-01-01", to = "2024-01-02",
                                  interval = 24, per_interval = 100,
                                  unit = "hours", tz = "UTC"), seed = 1)
  expect_equal(format(min(utc$ts), "%H:%M"), "19:00")

  # A column with an empty zone, shown in the session's zone, gets UTC.
  local <- ny
  attr(local$ts, "tzone") <- ""
  p <- inclusion_prob(local, des)
  expect_equal(format(range(local$ts[p > 0]), "%H:%M", tz = "UTC"),
               c("00:00", "23:00"))

  expect_output(print(des), "tz +from the time column, else UTC")
  expect_output(print(design_temporal("ts", "2024-01-01", "2024-01-02", 1, 1,
                                      tz = zone)),
                "tz +\"America/New_York\"")
})

test_that("bad zones and bounds are refused when the design is built", {
  expect_error(design_temporal("ts", "2024-01-01", "2024-01-02", 6, 2,
                               tz = "Mars/Olympus"),
               "`tz` is \"Mars/Olympus\", which is not a time zone R knows")
  expect_error(design_temporal("ts", "2024-01-01", "2024-01-02", 6, 2,
                               tz = c("UTC", "GMT")),
               "`tz` must be a single time zone name")
  expect_error(design_temporal("ts", c("2024-01-01", "2024-01-02"),
                               "2024-01-03", 6, 2),
               "`from` must be a single date-time, not 2 values")
  expect_error(design_temporal("ts", "yesterday", "2024-01-03", 6, 2),
               "`from` is \"yesterday\", which does not read as a date-time")
  expect_error(design_temporal("ts", "2024-01-01", NA, 6, 2),
               "`to` is NA")
  expect_error(design_temporal("ts", 5, "2024-01-03", 6, 2),
               "`from` must hold date-times")

  d <- data.frame(ts = .POSIXct(1704067200 + 3600 * 0:3, tz = "Mars/Olympus"))
  expect_error(
    inclusion_prob(d, design_temporal("ts", "2024-01-01", "2024-01-02", 6, 2)),
    "`ts` is in time zone \"Mars/Olympus\", which R does not know")
})

test_that("strings with hours and minutes but no seconds parse", {
  d <- data.frame(ts = c("2024-01-01 05:00", "2024-01-01T07:30",
                         "2024-01-01 09", "2024-01-01 11:00:00",
                         "2024-01-01"))
  des <- design_temporal("ts", from = "2024-01-01 00:00", to = "2024-01-02 00",
                         interval = 6, per_interval = 1)
  expect_equal(inclusion_prob(d, des), c(1 / 2, 1 / 3, 1 / 3, 1 / 3, 1 / 2))

  # A date written day or month first stays unread rather than guessed.
  expect_true(is.na(drawn:::parse_time("01/02/2024", "UTC", "x")))
  expect_error(
    inclusion_prob(data.frame(ts = c("01/02/2024", "2024-01-01")), des),
    "strings must be year-month-day with an optional time")

  # A string with its own offset keeps its instant in another zone.
  got <- drawn:::parse_time(c("2024-01-01T05:00:00Z", "2024-01-01 05:00"),
                            "America/New_York", "x")
  expect_equal(format(got, "%H:%M"), c("00:00", "05:00"))
})

test_that("months and years count from the anchor", {
  # seq() rolled 31 January on to 2 March, and 29 February on to 1 March.
  utc <- function(x) as.POSIXct(x, tz = "UTC")
  expect_equal(
    drawn:::make_breaks(utc("2024-01-31"), utc("2024-07-01"), 1, "months"),
    utc(c("2024-01-31", "2024-02-29", "2024-03-31", "2024-04-30",
          "2024-05-31", "2024-06-30")))
  expect_equal(
    drawn:::make_breaks(utc("2024-02-29"), utc("2028-03-01"), 1, "years"),
    utc(c("2024-02-29", "2025-02-28", "2026-02-28", "2027-02-28",
          "2028-02-29")))

  d <- day_frame("2024-01-31", 182)
  des <- design_temporal("day", from = "2024-01-31", to = "2024-07-31",
                         interval = 1, per_interval = 1, unit = "months")
  expect_equal(as.vector(table(drawn:::temporal_bucket(des, d))),
               c(29, 31, 30, 31, 30, 31))
})

test_that("inclusion probabilities of a long hourly frame take one pass", {
  skip_on_cran()
  # One pass per interval took about 6 seconds here.
  d <- data.frame(ts = seq(as.POSIXct("2024-01-01", tz = "UTC"), by = 150,
                           length.out = 2e5))
  des <- design_temporal("ts", from = "2024-01-01", to = "2024-12-01",
                         interval = 1, per_interval = 2, unit = "hours")
  time <- system.time(p <- inclusion_prob(d, des))[["elapsed"]]
  expect_lt(time, 2)
  inw <- d$ts < as.POSIXct("2024-12-01", tz = "UTC")
  expect_equal(p, ifelse(inw, 1 / 12, 0))
})

test_that("intervals are named by their start, in time order", {
  h <- data.frame(ts = seq(as.POSIXct("2024-01-01", tz = "UTC"), by = "hour",
                           length.out = 12))
  des <- design_temporal("ts", from = "2024-01-01", to = "2024-01-01 12:00",
                         interval = 1, per_interval = 1, unit = "hours")
  b <- drawn:::temporal_bucket(des, h)
  expect_identical(as.integer(b), 1:12)
  expect_identical(levels(b), sprintf("2024-01-01 %02d:00", 0:11))
  expect_identical(sort(levels(b), method = "radix"), levels(b))

  # Seconds appear when an interval starts between minutes.
  secs <- design_temporal("ts", from = "2024-01-01", to = "2024-01-01 00:02",
                          interval = 30, per_interval = 1, unit = "seconds")
  expect_identical(levels(drawn:::temporal_bucket(secs, h)),
                   paste("2024-01-01", c("00:00:00", "00:00:30", "00:01:00",
                                         "00:01:30")))

  # New York repeats 01:00 when its clocks go back, and the two intervals
  # starting then are told apart by the zone.
  zone <- "America/New_York"
  fall <- data.frame(ts = seq(as.POSIXct("2024-11-03 00:00", tz = zone),
                              by = "hour", length.out = 5))
  des <- design_temporal("ts", from = "2024-11-03 00:00",
                         to = "2024-11-03 04:00", interval = 1,
                         per_interval = 1, unit = "hours")
  b <- drawn:::temporal_bucket(des, fall)
  expect_identical(levels(b), c("2024-11-03 00:00", "2024-11-03 01:00 EDT",
                                "2024-11-03 01:00 EST", "2024-11-03 02:00",
                                "2024-11-03 03:00"))
  expect_identical(as.integer(b), 1:5)
})

test_that("notes and survey's strata name intervals by their start", {
  d <- data.frame(ts = rep(seq(as.POSIXct("2024-01-01", tz = "UTC"),
                               by = "day", length.out = 12), each = 4),
                  y = (1:48 * 7) %% 11)
  des <- design_temporal("ts", from = "2024-01-01", to = "2024-01-13",
                         interval = 1, per_interval = 1, unit = "days")
  s <- draw(d, des, seed = 1, weights = TRUE)
  # In time order: the intervals were numbered and sorted as text, 1, 10, 11.
  expect_match(ht_total(s, "y")$note, paste0(
    "12 intervals `2024-01-01 00:00`, `2024-01-02 00:00`, ",
    "`2024-01-03 00:00`, `2024-01-04 00:00`, `2024-01-05 00:00`, ... each"),
    fixed = TRUE)

  skip_if_not_installed("survey")
  two <- draw(d, design_temporal("ts", from = "2024-01-01", to = "2024-01-13",
                                 interval = 1, per_interval = 2,
                                 unit = "days"),
              seed = 1, weights = TRUE)
  sv <- as_svydesign(two)
  expect_setequal(unique(sv$variables$.stratum),
                  paste(format(as.Date("2024-01-01") + 0:11), "00:00"))
})
