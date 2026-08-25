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
  bad <- design_temporal("ts", from = "2020-01-02", to = "2020-01-01",
                         interval = 1, per_interval = 1, unit = "hours")
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
