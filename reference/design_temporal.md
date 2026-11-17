# Temporal sampling

Divides a time window into intervals of `interval` units and samples
within each one. `per_interval` is named for what it is: unlike `n`
elsewhere in the package, it is a per-group figure, not a total.

## Usage

``` r
design_temporal(
  time,
  from,
  to,
  interval,
  per_interval,
  unit = c("hours", "seconds", "minutes", "days", "weeks", "months", "years"),
  tz = NULL,
  na_rm = FALSE
)
```

## Arguments

- time:

  Column of timestamps: `POSIXct`, `POSIXlt`, `Date`, character or
  factor. Strings are read as year-month-day with an optional time, in
  the forms
  [`lubridate::ymd_hms()`](https://lubridate.tidyverse.org/reference/ymd_hms.html)
  accepts, such as `"2024-03-01 14:30:00"` or `"2024-03-01T14:30:00Z"`,
  and also `"2024-03-01 14:30"`, `"2024-03-01 14"` and `"2024-03-01"`.

- from, to:

  Window bounds, each a single `POSIXct`, `Date` or string read as
  `time` is. The window is half-open: `[from, to)`. Trailing partial
  intervals are dropped. Both are checked when the design is built, and
  `to` must come after `from`. Without `tz`, a bound that names an
  instant (a `POSIXct`, or a string ending in `Z` or an offset) paired
  with a clock time (a `Date`, or a string without an offset) can fall
  in either order depending on the zone. Such a pair is refused when the
  design is built if, with the clock time read in UTC, `to` comes a day
  or more before `from`, since no zone's offset from UTC can undo that.
  Otherwise its order is checked when the design is used, in the zone
  the window is laid out in (see "Time zones and calendar steps").

- interval:

  Interval width, in units of `unit`.

- per_interval:

  Rows to draw from each interval. Intervals holding fewer rows
  contribute all of them.

- unit:

  One of `"seconds"`, `"minutes"`, `"hours"`, `"days"`, `"weeks"`,
  `"months"` or `"years"`. See "Time zones and calendar steps".

- tz:

  Time zone, a name from
  [`OlsonNames()`](https://rdrr.io/r/base/timezones.html), in which the
  bounds, any `Date` or string timestamps and the calendar steps are
  read. `NULL`, the default, takes the zone of a `POSIXct` `time` column
  and UTC otherwise. It shifts the window rather than relabelling the
  output.

- na_rm:

  Drop rows whose timestamp is missing or unparseable instead of raising
  an error.

## Value

A design object, for use with
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

## Details

Each row's interval is found in a single pass over the data, so the draw
and the inclusion probabilities take time that grows with the number of
rows rather than with rows times intervals.

## Time zones and calendar steps

The window and the intervals are laid out in one time zone: `tz` when it
is given, otherwise the zone of a `POSIXct` `time` column, otherwise
UTC. A `POSIXct` column with an empty `tzone`, which R shows in the
session's zone, also gets UTC, so a design gives the same intervals on
every machine.

A `Date`, in the column or as a bound, is its midnight in that zone, and
a string without an offset is a clock time there. A clock time that a
daylight-saving change skips moves to the moment the clocks jump, and
one that happens twice takes its first occurrence, so a `Date` stands
for the first instant of its day: 01:00 where the change skips midnight,
as in Santiago or Cairo, and the first of two midnights where it repeats
one, as in the Azores. A string ending in `Z` or an offset such as
`+01:00` keeps the instant it names.

`"seconds"`, `"minutes"` and `"hours"` are fixed durations counted from
`from`, so hourly intervals ignore the clock change. `"days"` and
`"weeks"` step by calendar days, so each interval starts at the same
clock time as `from` and lasts 23 or 25 hours across a change.
`"months"` and `"years"` step by calendar units counted from `from`, and
a day the month lacks moves to its last day: a window starting on 31
January 2024 has intervals starting on 29 February, 31 March and 30
April.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md),
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md),
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md),
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
[`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md),
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)

## Examples

``` r
df <- data.frame(
  id = 1:48,
  ts = seq(as.POSIXct("2020-01-01", tz = "UTC"), by = "hour", length.out = 48)
)
d <- design_temporal("ts", from = "2020-01-01", to = "2020-01-02",
                     interval = 6, per_interval = 2, unit = "hours")
draw(df, d, seed = 1)$ts
#> [1] "2020-01-01 00:00:00 UTC" "2020-01-01 03:00:00 UTC"
#> [3] "2020-01-01 06:00:00 UTC" "2020-01-01 07:00:00 UTC"
#> [5] "2020-01-01 14:00:00 UTC" "2020-01-01 16:00:00 UTC"
#> [7] "2020-01-01 19:00:00 UTC" "2020-01-01 23:00:00 UTC"
```
