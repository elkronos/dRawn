# Temporal sampling

Divides a time window into equal intervals and samples within each one.
`per_interval` is named for what it is: unlike `n` elsewhere in the
package, it is a per-group figure, not a total.

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

  Column of timestamps. `POSIXct` and `Date` columns are used as they
  are; character columns are parsed with
  [`lubridate::ymd_hms()`](https://lubridate.tidyverse.org/reference/ymd_hms.html).

- from, to:

  Window bounds, as timestamps or strings. The window is half-open:
  `[from, to)`. Trailing partial intervals are dropped.

- interval:

  Interval width, in units of `unit`.

- per_interval:

  Rows to draw from each interval. Intervals holding fewer rows
  contribute all of them.

- unit:

  One of `"seconds"`, `"minutes"`, `"hours"`, `"days"`, `"weeks"`,
  `"months"` or `"years"`. `"months"` and `"years"` step by calendar
  units, not by fixed 30.44-day durations.

- tz:

  Time zone the timestamps are expressed in. Applied when the column
  *and* both bounds are parsed, so it genuinely shifts the window rather
  than relabelling the output.

- na_rm:

  Drop rows whose timestamp is missing or unparseable instead of raising
  an error.

## Value

A design object, for use with
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

## Details

Buckets are assigned arithmetically in a single pass rather than by
re-filtering the data per interval, so cost scales with the data rather
than with the number of intervals.

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
