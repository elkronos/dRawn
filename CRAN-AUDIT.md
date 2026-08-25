# sampleR — pre-CRAN audit

**52 findings across 10 sampling routines.** Every claim below was
reproduced by sourcing `R/utils.R` and `R/sampling_functions.R` into a
live R session — none are inferred from reading.

|  |  |
|----|----|
| Verified against | R 4.3.3 · dplyr 1.1.4 · sf 1.0-15 · lubridate 1.9.3 · testthat 3.2.1 |
| Severity split | 22 blockers · 22 should-fix · 8 polish |
| Status | 50 of 52 fixed in `drawn` 0.1.0; F9 and F11 partly done |
| `R CMD check --as-cran` | 2 substantive WARNINGs, 1 NOTE (on a hand-built skeleton) |
| Tests discovered by testthat | 0 — `UAT.R` does not match `^test` |
| Current suite | 8 of 10 blocks pass; temporal fails and spatial errors, both every run |

Two things I expected to be problems were **not**, and are deliberately
absent below:

- The reservoir sampler’s Algorithm R is a correct uniform sample (χ² p
  = 0.64 over 20,000 replications).
- `do.call(rbind, …)` in `multi_stage_sampling()` cost only 0.15 s at
  500 clusters.

The findings were then re-derived independently by a second pass told to
break them. That pass added seven — A15, A16, E10 and the extra cases
folded into A8, A11, B2 and C2 — and corrected four claims drawn too
narrowly.

------------------------------------------------------------------------

## A. Silent wrong answers (16)

These do not error. They return a data frame that looks right and is
not.

### A1 · Blocker — the `sample()` length-one trap silently invents clusters

`cluster_sampling`, `multi_stage_sampling`

`sample(clusters, num_clusters)` changes meaning when `clusters` is a
length-one numeric: R samples from `1:clusters` instead of from the
vector.

``` r

d <- data.frame(x = 1:10, cl = rep(5, 10))
cluster_sampling(d, "cl", num_clusters = 1, seed = 1)
#> [1] x  cl
#> <0 rows> (or 0-length row.names)

multi_stage_sampling(d, "cl", num_clusters = 1, stage_two_sample_size = 2, seed = 1)
#> Error: Sample size requested for cluster 1 exceeds available rows.   # cluster 1 does not exist
```

**Fix** —
`sampled_clusters <- clusters[sample.int(length(clusters), num_clusters)]`

### A2 · Blocker — `max_rows` takes `head()` of sorted output, discarding whole strata

`apply_max_rows` → stratified · cluster · multi_stage · temporal

`stratified_sampling()` returns rows grouped and sorted by stratum, so
capping keeps only the first one.

``` r

table(stratified_sampling(df, "group", sample_size = 40, max_rows = 10, seed = 1)$group)
#>  a
#> 10        # groups b, c, d vanished entirely

table(cluster_sampling(df, "group", num_clusters = 4, max_rows = 10, seed = 1)$group)
#>  a
#> 10
```

Same in `temporal_sampling()`: a 24-hour window capped at 4 rows returns
only the first two intervals.

**Fix** — take a random subset rather than the head, or better, push the
cap into the allocation step so each stratum scales down proportionally.

``` r

apply_max_rows <- function(data, max_rows = NULL) {
  if (is.null(max_rows)) return(data)
  if (nrow(data) <= max_rows) return(data)
  data[sort(sample.int(nrow(data), max_rows)), , drop = FALSE]
}
```

### A3 · Blocker — proportional stage-two divides by the wrong denominator

`multi_stage_sampling`

`n_sample` is scaled by `nrow(cluster_data) / nrow(data)` — the whole
population — but only the *selected* clusters are available. The error
runs both directions: undershooting when few clusters are selected, and
— because `max(1, …)` floors every cluster at one row — overshooting
badly when many small clusters are selected.

``` r

# undershoot
nrow(multi_stage_sampling(df, "group", num_clusters = 2,
     stage_two_sample_size = 20, proportional_stage_two = TRUE, seed = 1))
#> [1] 10      # asked for 20

# overshoot — 100 clusters of 10
d <- data.frame(x = 1:1000, cl = as.character(rep(1:100, each = 10)))
nrow(multi_stage_sampling(d, "cl", num_clusters = 20,
     stage_two_sample_size = 5, proportional_stage_two = TRUE, seed = 1))
#> [1] 20      # asked for 5 — every cluster floored up to 1
# num_clusters = 50, stage_two_sample_size = 10  ->  50 rows (asked for 10)
```

**Fix** — divide by the row count across the selected clusters, not
`nrow(data)`.

### A4 · Blocker — weighted sampling hands back rewritten weights

`weighted_sampling`

Normalisation is applied to `data_copy`, and `data_copy` is what gets
returned.

``` r

d   <- data.frame(id = 1:10, w = 1:10)
out <- weighted_sampling(d, "w", sample_size = 5, normalization = "min-max", seed = 1)
out$w
#> [1] 0.8888889 0.7777778 0.5555556 0.2222222 1.0000000
d$w[match(out$id, d$id)]
#> [1]  9  8  6  3 10       # what should have come back
```

**Fix** — normalise into a local vector, then index the *original*
`data`.

### A5 · Blocker — both normalisations make the lightest unit unselectable

`weighted_sampling`

Min–max maps the smallest weight to exactly 0; `+ .Machine$double.eps`
leaves it at ~2.2e-16. **The z-score branch has the identical defect** —
line 129 also subtracts the minimum before adding the same epsilon.
Subtracting the minimum is an affine shift, so it also changes every
*relative* probability — this is not a neutral rescaling.

``` r

d <- data.frame(id = 1:5, w = c(1, 2, 3, 4, 5))

table(replicate(2000, weighted_sampling(d, "w", 1, normalization = "min-max")$id))
#>   2   3   4   5
#> 188 398 595 819       # id 1 drawn 0 times in 2,000

table(replicate(2000, weighted_sampling(d, "w", 1, normalization = "z-score")$id))
#>   2   3   4   5
#> 213 390 586 811       # id 1 drawn 0 times; its prob is 2.22e-16
```

**Fix** — drop `normalization` (`prob=` is already scale-invariant so it
buys nothing legitimate), or replace with a ratio-preserving positive
map and document that min–max is a deliberate re-weighting.

### A6 · Blocker — balanced cluster sampling returns nothing for factor columns

`cluster_sampling(balanced = TRUE)`

`min(as.vector(table(...)))` counts unused factor levels as zero, so
`slice_sample(n = 0)` empties the result.

``` r

d <- data.frame(x = 1:30,
      cl = factor(rep(c("a","b","c"), each = 10), levels = c("a","b","c","zz")))
nrow(cluster_sampling(d, "cl", num_clusters = 2, balanced = TRUE, seed = 1))
#> [1] 0
```

**Fix** —
`min(table(droplevels(as.factor(sampled_data[[cluster_column]]))))`,
then guard `< 1`.

### A7 · Blocker — the block bootstrap draws every row from a single block

`bootstrap_sampling(method = "block")`

It picks *one* block and resamples `sample_size` rows with replacement
from inside it. A real block bootstrap concatenates ⌈n/ℓ⌉ independently
chosen blocks — that is what preserves the serial dependence you are
bootstrapping over.

``` r

r <- bootstrap_sampling(data.frame(id = 1:100), num_samples = 3,
                        sample_size = 10, method = "block", seed = 1)
lapply(r, function(z) range(z$id))
#> [[1]] 81  87
#> [[2]] 95 100
#> [[3]] 11  20      # each replicate confined to one 10-row block
```

**Fix** — separate `block_length` from `sample_size`, use moving blocks:

``` r

l      <- block_length %||% max(1L, floor(nrow(data)^(1/3)))
nb     <- ceiling(sample_size / l)
starts <- sample.int(nrow(data) - l + 1L, nb, replace = TRUE)
idx    <- unlist(lapply(starts, function(s) s:(s + l - 1L)))[seq_len(sample_size)]
data[idx, , drop = FALSE]
```

### A8 · Blocker — `max_rows` silently changes the number of bootstrap replicates

`bootstrap_sampling`

The branch rbinds all replicates, de-duplicates, truncates and
re-splits. Deduplicating a bootstrap sample is wrong on its face —
repeated draws *are* the bootstrap.

``` r

r <- bootstrap_sampling(df, num_samples = 5, sample_size = 10, max_rows = 20, seed = 1)
length(r)
#> [1] 2      # asked for 5 replicates

sapply(bootstrap_sampling(d, num_samples = 5, sample_size = 10, max_rows = 25, seed = 1), nrow)
#>  1  2  3
#> 10 10  5      # and the replicates come back ragged
```

**Fix** — remove `max_rows` here; cap `sample_size` per replicate
instead.

### A9 · Blocker — temporal sampling re-parses columns that are already date-times

`temporal_sampling`

[`lubridate::ymd_hms()`](https://lubridate.tidyverse.org/reference/ymd_hms.html)
is applied unconditionally. A `POSIXct` column round-trips through
[`as.character()`](https://rdrr.io/r/base/character.html), midnight
values lose their time component and become `NA`. **This is why your own
UAT currently fails.**

``` r

d$time <- as.POSIXct(d$time, tz = "UTC")
temporal_sampling(d, "time", "2020-01-01 00:00:00", "2020-01-02 00:00:00",
                  interval = 6, sample_size = 2, unit = "hours", seed = 1)
#> Warning: 2 failed to parse.

d$time <- seq(as.Date("2020-01-01"), by = "day", length.out = 10)   # Date
#> Warning: All formats failed to parse. No formats found.
```

**Fix**

``` r
tv <- data[[time_column]]
tv <- if (inherits(tv, "POSIXt")) tv
      else if (inherits(tv, "Date")) as.POSIXct(tv, tz = tz %||% "UTC")
      else lubridate::ymd_hms(tv, tz = tz %||% "UTC")
if (anyNA(tv)) stop("`", time_column, "` contains unparseable date-times.", call. = FALSE)
```

### A10 · Blocker — `complex_region` breaks with more than two polygons

`spatial_sampling`

`do.call(sf::st_union, region)` maps the third and later elements onto
`st_union()`’s `by_feature` / `...` arguments.

``` r

spatial_sampling(sdf, "lat", "lon", region = list(p1, p2, p3),
                 sample_size = 30, complex_region = TRUE, seed = 1)
#> Error: `options` must be created using s2_options()
```

**Fix** —
`region <- sf::st_union(do.call(c, lapply(region, sf::st_geometry)))`

### A11 · Blocker — `max_rows` accepts three kinds of nonsense without complaint

`apply_max_rows`

`head(x, -3)` drops the last three rows; a non-numeric value is ignored
entirely; and because `head.data.frame()` accepts a vector `n`, a
length-two `max_rows` silently *drops columns*.

``` r

nrow(apply_max_rows(data.frame(x = 1:10), -3))    #> [1] 7    # trimmed from the end
nrow(apply_max_rows(data.frame(x = 1:10), "5"))   #> [1] 10   # silently ignored
dim(apply_max_rows(data.frame(a=1:10, b=1:10, c=1:10,
                   d=1:10, e=1:10, f=1:10), c(3, 2)))
#> [1] 3 2      # four columns gone
```

**Fix** — validate: single, non-negative, whole, non-`NA`, numeric.

### A12 · Blocker — missing `drop = FALSE` mangles single-column inputs

`cluster_sampling`, `multi_stage_sampling`, `systematic_sampling`

``` r

d <- data.frame(cl = rep(c("a","b"), each = 5))
cluster_sampling(d, "cl", num_clusters = 1, seed = 1)
#>   sampled_data
#> 1            a
#> 2            a
#> 3            a
#> 4            a
#> 5            a           # the column is renamed to the internal variable

systematic_sampling(data.frame(v = c(3,1,2)), interval = 1, sort_column = "v", start = 1)
#> Error: 'to' must be of length 1
```

**Fix** — add `, drop = FALSE` to all three subset calls.

### A13 · Should fix — the reservoir is padded with `NULL`s when the stream is short

`reservoir_sampling`

``` r

r <- reservoir_sampling(data.frame(x = 1:5), sample_size = 10, seed = 1)
c(length = length(r), n_null = sum(vapply(r, is.null, logical(1))))
#> length n_null
#>     10      5
```

**Fix** — `reservoir <- reservoir[seq_len(min(count, sample_size))]`

### A14 · Should fix — an empty time window returns a 0×0 frame with the columns gone

`temporal_sampling`

``` r

r <- temporal_sampling(tdf, "time", "2021-01-01 00:00:00", "2021-01-02 00:00:00",
                       interval = 6, sample_size = 2, unit = "hours", seed = 1)
dim(r)
#> [1] 0 0      # expected 0 x ncol(tdf)
```

**Fix** — `if (length(parts) == 0L) return(data[0L, , drop = FALSE])`.
Same guard in `multi_stage_sampling()`.

### A15 · Blocker — `==` instead of `%in%` fabricates rows that were never in the data

`multi_stage_sampling`

Stage two selects with `data[data[[cluster_column]] == cluster, ]`. In
R, `NA == "a"` is `NA`, and indexing a data frame with `NA` returns a
row of all-`NA` values. Every missing cluster label becomes a phantom
row inside *every* cluster’s stage-two frame — and those rows can be
sampled and returned. `cluster_sampling()` uses `%in%` and is safe; this
one is not.

``` r

d <- data.frame(x = 1:30, cl = c(rep("a",10), rep("b",10), rep(NA,10)))
c(truth = nrow(d[d$cl %in% "a", ]), buggy = nrow(d[d$cl == "a", ]))
#> truth buggy
#>    10    20

multi_stage_sampling(d, "cl", num_clusters = 1, stage_two_sample_size = 5)
#>        x   cl
#> NA.18 NA <NA>
#> NA.2  NA <NA>
#> NA.6  NA <NA>
#> NA.11 NA <NA>
#> NA.29 NA <NA>      # all five rows invented
```

The inflated count also defeats the oversampling guard, so
`replace = FALSE` stops protecting you:

``` r

r <- multi_stage_sampling(d, "cl", num_clusters = 1, stage_two_sample_size = 15, replace = FALSE)
c(nrow = nrow(r), n_fabricated = sum(is.na(r$x)))
#>         nrow n_fabricated
#>           15            9      # 15 rows drawn from a 10-row cluster
```

**Fix** —
`cluster_data <- data[data[[cluster_column]] %in% cluster, , drop = FALSE]`,
then decide explicitly what `NA` means as a cluster label (C7).

### A16 · Blocker — regions crossing the antimeridian collapse to a sliver

`spatial_sampling`

With `sf_use_s2() == TRUE`, consecutive polygon vertices are joined by
the *shortest* great-circle path. An edge from longitude −179 to +179
therefore spans the 2° across the antimeridian, not the 358° the
coordinates suggest. A “whole world” rectangle collapses into a 2°-wide
pole-to-pole strip and contains almost nothing. No warning is emitted.

This is not hypothetical: it is the polygon in your own UAT, and it is
why the spatial test errors on every run.

``` r

gp <- st_sfc(st_polygon(list(cbind(c(-179,179,179,-179,-179),
                                   c(-89,-89,89,89,-89)))), crs = 4326)
as.numeric(st_area(gp)) / 1e12
#> [1] 2.833269   # million km2; a 2-degree pole-to-pole strip is 2/360 * 510 = 2.83

sum(lengths(st_intersects(world_points_sf, gp)) > 0)
#> [1] 0          # of 100 points spread over the whole world
```

**Correction to an earlier draft of this finding.** I first attributed
this to ring winding and suggested reversing the orientation. That is
wrong, and testing it is how I found out: `sf` normalises winding, so
clockwise and counter-clockwise give the identical 27.9 M km2 strip for
a -170..170 box. Adding intermediate vertices does not reliably fix it
either (7 of 200 world points matched). The two things that do work are
splitting the region at the antimeridian, or `sf::sf_use_s2(FALSE)`,
which puts all 200 back inside.

**Fix — detect on the coordinates, not on the area**

``` r

# An area-vs-bounding-box test does NOT work: a collapsed ring and its bbox collapse
# together, so the antimeridian polygon scores a ratio of 1.00, while three legitimate
# disjoint boxes score 0.12 -- exactly backwards.
span <- max(abs(diff(sf::st_coordinates(region)[, "X"])))   # per ring
if (isTRUE(sf::sf_use_s2()) && sf::st_is_longlat(region) && span > 180)
  warning("An edge of `region` spans ", round(span), " degrees of longitude. Under ",
          "spherical geometry that edge is drawn the short way, across the antimeridian, ",
          "so `region` is much smaller than its coordinates suggest.", call. = FALSE)
if (nrow(region_data) == 0L) stop("No rows fall inside `region`.", call. = FALSE)
```

The zero-points case currently reports “Sample size cannot be larger
than the number of points within the region”, which reads as a sizing
problem rather than an empty selection.

------------------------------------------------------------------------

## B. Statistical accuracy (6)

The code runs and returns plausible output, but the sample does not have
the properties its name implies.

### B1 · Blocker — weighted sampling without replacement is not πPS

`weighted_sampling(replace = FALSE)`

`sample(prob = w, replace = FALSE)` performs *successive* sampling: the
weights govern each sequential draw, not the probability that a row ends
up in the sample.

**Correction to an earlier draft of this finding.** My first version
compared the realised rates against `n × p_i` — 0.10 for the light
units, **1.60** for the heavy one — and reported the light units as
included “2.6× their nominal rate”. That target is not reachable: a
probability cannot exceed 1, so a fixed-size-2 design must take the
heavy unit with certainty and spread the remaining slot over the other
four at 0.25 each. Against the *feasible* πPS values, successive
sampling in that example is about 4% off. I was measuring the gap to an
impossible target.

    30,000 replications, n = 2 from 5, weights (1,1,1,1,16)

    naive  n*p     0.100  0.100  0.100  0.100  1.600     # not a probability
    feasible piPS  0.250  0.250  0.250  0.250  1.000
    successive     0.261  0.255  0.259  0.256  0.968     # 4.2% off

The finding stands; it needed a fair example. The gap depends on the
weight spread, and one dominant unit is the case where the two nearly
coincide:

    weights 1..10, n = 5
    feasible piPS  0.091 0.182 0.273 0.364 0.455 0.545 0.636 0.727 0.818 0.909
    successive     0.122 0.229 0.339 0.431 0.502 0.579 0.634 0.684 0.720 0.759
                   # up to 34.5% away -- light units over-included, heavy under

    weights 1,2,4,8,16,32, n = 3   max 30.8% away
    weights 1,1,1,1,1,50,  n = 3   max  1.3% away -- nearly identical

And it does bias estimation:

    30-row population, cor(w, y) = 0.84, 6,000 replications, n = 15, true = 815

    successive treated as piPS   823.1   (+1.0%, 11.0 SE)
    systematic PPS               814.7   ( -0.0%,  0.4 SE)

Modest in size, unmistakable in sign. An earlier version of this check
set `w = y`, which makes `y/pi` constant and the estimator exact for
*any* pi — it could not have detected bias at all.

**Fix** — add a `method` argument: keep `"successive"` as the default,
and add systematic PPS and Poisson sampling, whose inclusion
probabilities are exactly `n * p_i`. Handle units with `n * p_i > 1` by
taking them with certainty and rescaling the rest, iterating, since
removing units raises everyone else’s share.

### B2 · Should fix — proportional allocation misses the target and floors rare strata upward

`stratified_sampling`

``` r

d <- data.frame(x = 1:143, g = rep(letters[1:7], times = c(31,29,23,19,17,13,11)))
nrow(stratified_sampling(d, "g", sample_size = 50, seed = 1))
#> [1] 51

d <- data.frame(x = 1:1000, g = c(rep("big", 999), "rare"))
table(stratified_sampling(d, "g", sample_size = 10, seed = 1)$g)
#>  big rare
#>   10    1      # 11 rows; "rare" at 9% vs 0.1% of the population
```

It misses low as readily as high —
[`round()`](https://rdrr.io/r/base/Round.html) uses banker’s rounding,
so exact halves all go down together:

``` r

d <- data.frame(x = 1:100, g = rep(c("a","b","c","d"), each = 25))
nrow(stratified_sampling(d, "g", sample_size = 10, seed = 1))
#> [1] 8      # round(2.5) = 2 in all four strata
```

**Fix** — largest-remainder allocation hits the target exactly:

``` r

alloc <- function(n, sizes, min_n = 0L) {
  raw   <- n * sizes / sum(sizes)
  base  <- pmax(min_n, floor(raw))
  short <- n - sum(base)
  if (short > 0) {
    ord <- order(raw - floor(raw), decreasing = TRUE)
    base[ord[seq_len(short)]] <- base[ord[seq_len(short)]] + 1L
  }
  pmin(base, sizes)
}
```

Expose `min_n` (default `0`) so the coverage-vs-bias tradeoff is the
caller’s, and add `allocation = c("proportional", "equal", "neyman")`.

### B3 · Should fix — `unit = "months"` uses 30.4375-day durations, not calendar months

`temporal_sampling`

``` r

as.numeric(lubridate::dmonths(1), "days")
#> [1] 30.4375

seq(ymd_hms("2020-01-01 00:00:00"), ymd_hms("2020-04-01 00:00:00") - dmonths(1), by = dmonths(1))
#> "2020-01-01"  "2020-01-31 10:30:00"     # second bucket starts mid-January
```

**Fix** — use `seq(start, end, by = "month")` for calendar units. Add
`"seconds"`, `"minutes"`, `"years"` too.

### B4 · Should fix — `time_zone` cannot change which rows fall in the window

`temporal_sampling`

`with_tz()` relabels an instant; it does not move it.
`start_time`/`end_time` are parsed as UTC regardless, so the filter is
identical for every `time_zone`.

``` r

identical(utc$id, chicago$id)   #> TRUE
identical(utc$id, tokyo$id)     #> TRUE
# only the display changes:
# UTC "2020-01-01"  Chicago "2019-12-31 18:00:00"  Tokyo "2020-01-01 09:00:00"
```

**Fix** — parse the column *and* both bounds with `tz = time_zone`.

### B5 · Should fix — no design weights or inclusion probabilities come back

all ten functions

Without inclusion probabilities the user cannot form an unbiased
estimate from a stratified, cluster, multi-stage or weighted sample —
which is the reason to draw one. This is the biggest gap between the
package as it stands and something a survey statistician would adopt.

``` r

structure(out,
  sampling_design = list(method = "stratified", strata = strata_column,
                         n_pop = N_h, n_samp = n_h),
  prob   = pi_i,          # inclusion probability per returned row
  weight = 1 / pi_i,      # design weight
  class  = c("sampleR_sample", class(out)))
```

With [`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) methods on that class
you get a coherent package identity and a natural bridge to
`survey::svydesign()`.

### B6 · Should fix — block length conflated with sample size; last block is short

`bootstrap_sampling(method = "block")`

n = 95, `sample_size` = 10:

       start end len
    9     81  90  10
    10    91  95   5      # 10 rows drawn with replacement from 5

Fix follows A7.

------------------------------------------------------------------------

## C. Input validation and error messages (7)

Bad input currently surfaces from deep inside
[`base::sample()`](https://rdrr.io/r/base/sample.html),
[`seq()`](https://rdrr.io/r/base/seq.html) or `sf`.

### C1 · Should fix — `sample_size` is never validated

``` r

dim(simple_random_sampling(data.frame(x = 1:10), 0))       #> [1] 0 1   silent
nrow(simple_random_sampling(data.frame(x = 1:10), 2.7))    #> [1] 2     silently truncated
simple_random_sampling(data.frame(x = 1:10), -1)           #> Error: invalid 'size' argument
simple_random_sampling(data.frame(x = 1:10), NA)           #> Error: missing value where TRUE/FALSE needed
simple_random_sampling(data.frame(x = 1:20), c(5, 10))     #> Error: 'length = 2' in coercion to 'logical(1)'
```

**Fix** — one helper called at the top of all ten functions:

``` r

check_count <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 || x != trunc(x))
    stop("`", arg, "` must be a single non-negative whole number.", call. = FALSE)
  as.integer(x)
}
```

### C2 · Should fix — `interval` and `start` are unchecked

`systematic_sampling`

``` r

systematic_sampling(df, interval = 0, start = 1)      #> Error: invalid '(to - from)/by'
systematic_sampling(df, interval = 10, start = 150)   #> Error: wrong sign in 'by' argument
```

The silent case matters more than either error — a non-integer
`interval` produces an invalid design with uneven gaps, and says
nothing:

``` r

systematic_sampling(data.frame(id = 1:20), interval = 2.5, start = 1)$id
#> [1]  1  3  6  8 11 13 16 18      # gaps of 2, 3, 2, 3, … — not a fixed interval
```

Require `interval` to be a whole number `>= 1`; require `start` in
`1:interval` (that is what makes it a valid systematic design); return
zero rows with a warning when `start > nrow(data)`.

### C3 · Should fix — a window shorter than one interval throws from `seq()`

`temporal_sampling`

``` r

temporal_sampling(tdf, "time", "2020-01-01 00:00:00", "2020-01-01 03:00:00",
                  interval = 6, sample_size = 2, unit = "hours", seed = 1)
#> Error: wrong sign in 'by' argument
```

Detect and explain it. Also document that partial trailing intervals are
dropped by design — a reasonable choice that is currently invisible.

### C4 · Should fix — CRS mismatch and NA coordinates surface as raw `sf` errors

`spatial_sampling`

``` r

region <- st_transform(p1, 3857)   #> Error: st_crs(x) == st_crs(y) is not TRUE
region <- st_set_crs(p1, NA)       #> Error: st_crs(x) == st_crs(y) is not TRUE
d$lat[1] <- NA                     #> Error: missing values in coordinates not allowed
```

**Fix** — validate the CRS, transform when it differs, check for NA
coordinates up front. Add a `crs` argument; hard-coding 4326 rules out
projected data.

### C5 · Should fix — zero-variance weights produce `NaN` or `Inf`

`weighted_sampling` — all weights = 3:

``` r

normalization = "min-max"   #> Error: NA in probability vector
normalization = "z-score"   #> Warning: no non-missing arguments to min; returning Inf
```

Detect constant weights; reject `NA` weights up front (`na.rm = TRUE` in
the min/max hides them until
[`sample()`](https://rdrr.io/r/base/sample.html) fails).

### C6 · Polish — the same condition produces different messages in different functions

Both: `sample_size` 10 from 5 rows, `replace = FALSE`.

    simple_random  "Sample size cannot be larger than the population when sampling without replacement."
    weighted       "cannot take a sample larger than the population when 'replace = FALSE'"   # base R leaked

Related: `group_modify()` strips the grouping column, so the stratified
error reads `"…rows in stratum "` with nothing after it. Route every
check through shared helpers and add `call. = FALSE`.

### C7 · Polish — `NA` handling in key columns is undefined

``` r

d <- data.frame(x = 1:30, cl = c(rep("a",10), rep("b",10), rep(NA,10)))
table(cluster_sampling(d, "cl", num_clusters = 3, seed = 1)$cl, useNA = "ifany")
#>    a    b <NA>
#>   10   10   10

systematic_sampling(data.frame(v = c(5,NA,4,2,3)), interval = 2, sort_column = "v", start = 1)$v
#> [1]  2  4 NA
```

Add `na_rm = FALSE`; drop with a message or stop. Either way, document
it.

------------------------------------------------------------------------

## D. Performance (2)

### D1 · Blocker — reservoir sampling materialises the entire stream first

`reservoir_sampling`

`split(data_stream, seq_len(nrow(data_stream)))` builds one data frame
per row before a single item is processed. This inverts the whole
purpose of reservoir sampling.

200,000 rows, `sample_size` = 100:

    reservoir_sampling()      6.54 – 7.92 s
    vectorised equivalent     < 0.001 s        ~7,000x

The algorithm itself is correct (χ² p = 0.64 over 20,000 replications) —
only the input handling is wrong.

**Fix**

``` r

# 1. a data.frame is not a stream — take the fast path
if (is.data.frame(data_stream))
  return(data_stream[sort(sample.int(nrow(data_stream), sample_size)), , drop = FALSE])
# 2. accept what a real stream looks like: a connection or a generator function
# 3. use Algorithm L — O(k log(n/k)) skips instead of one RNG draw per item
```

### D2 · Should fix — temporal sampling re-scans the full table once per interval

`temporal_sampling`

The `lapply` over `time_seq` runs a fresh `dplyr::filter()` across all
rows for every bucket.

50,000 rows, 833 hourly buckets:

    current                 3.11 s
    bucket index + split    0.043 s     ~72x

**Fix**

``` r

keep   <- tv >= start_dt & tv < end_dt
d      <- data[keep, , drop = FALSE]
bucket <- findInterval(tv[keep], breaks)
idx    <- unlist(lapply(split(seq_len(nrow(d)), bucket), function(i)
            if (length(i) <= sample_size) i else sample(i, sample_size)))
d[sort(idx), , drop = FALSE]
```

------------------------------------------------------------------------

## E. API and design consistency (10)

Unifying ten scripts into one package is mostly this section.

### E1 · Should fix — three different return types across ten functions

`simple_random_sampling()` and `weighted_sampling()` preserve the input
class; `stratified`, `cluster`, `multi_stage`, `temporal` force
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html);
`reservoir` and `bootstrap` return lists.

``` r

class(simple_random_sampling(tibble(x = 1:10), 3, seed = 1))
#> [1] "tbl_df" "tbl" "data.frame"
class(stratified_sampling(tibble(...), ...))
#> [1] "data.frame"      # tibble silently downgraded
```

Pick one rule — “same class in, same class out” is least surprising —
and apply it everywhere.

### E2 · Should fix — `max_rows` is redundant where it’s safe and harmful where it isn’t

``` r

identical(simple_random_sampling(d, 20, max_rows = 5, seed = 3)$x,
          simple_random_sampling(d,  5, seed = 3)$x)
#> [1] TRUE
```

Recommendation: remove `max_rows` from the public API entirely. If you
keep it, fix A2 and A11 first and document precisely what it does per
function.

### E3 · Should fix — `equal_samples` redefines what `sample_size` means

`stratified_sampling`

With `FALSE`, `sample_size` is the total across strata; with `TRUE` it
becomes per-stratum. Same argument, two units. Replace with
`allocation = c("proportional", "equal", "neyman")` and keep
`sample_size` as the total. `multi_stage_sampling()` has the same
problem with `proportional_stage_two`.

### E4 · Polish — only one strata or cluster column is supported

``` r

stratified_sampling(d, c("g", "h"), sample_size = 10, seed = 1)
#> Error: Must subset the data pronoun with a string, not a character vector.
```

**Fix** —
`dplyr::group_by(data, dplyr::across(dplyr::all_of(strata_column)))`

### E5 · Blocker — passing `seed` permanently alters the caller’s RNG stream

Every function calls [`set.seed()`](https://rdrr.io/r/base/Random.html)
on the global RNG and never restores it.

``` r

set.seed(42); a <- runif(1)
set.seed(42); invisible(simple_random_sampling(d, 2, seed = 99)); b <- runif(1)
c(a, b)
#> 0.9148060   0.5349936      # the caller's stream moved
```

**Fix**

``` r

with_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  if (exists(".Random.seed", .GlobalEnv)) {
    old <- get(".Random.seed", .GlobalEnv)
    on.exit(assign(".Random.seed", old, .GlobalEnv), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)
  force(code)
}
```

([`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html)
does exactly this if you’d rather take the dependency.)

### E6 · Blocker — three placeholder utilities would ship as public API

`utils.R`

`setup_logging()` prints a message and does nothing. `set_seed()` is a
one-line wrapper. `get_random_state()` reads `.Random.seed` from the
global environment and errors if the RNG has never been used.
`validate_dataframe()` and `apply_max_rows()` are internal helpers.

Delete the first three; keep the last two internal. That also removes
five entries from the “undocumented objects” warning in F3.

### E7 · Polish — row names carry over from the population

A sample of three rows is numbered `"2" "5" "1"`. Reset them, or add an
option to keep the original index as a real column — which is what
people actually want it for.

### E8 · Polish — `handle_infinite_stream` hides a hard-coded 10,000

`reservoir_sampling` — 25,000-row input:

``` r
range(ids drawn with handle_infinite_stream = TRUE)
#> [1] 2653 8364      # nothing past row 10,000 was ever seen
```

Replace with `max_items = Inf`, and warn when the cap truncates.

### E9 · Polish — `sample_size` should default to `nrow(data)` for the bootstrap

The standard nonparametric bootstrap resamples n from n. Default it, and
default `num_samples` to something conventional like 1,000.

### E10 · Should fix — stratified sampling silently reorders your columns

`dplyr::group_modify()` returns the grouping key first, so the strata
column jumps to position 1. Positional indexing breaks, and
[`rbind()`](https://rdrr.io/r/base/cbind.html)-ing the sample against
the population frame fails or misaligns. No other function does this.

``` r

names(d)                                        #> "id" "value" "group"
names(stratified_sampling(d, "group", 8, seed = 1))
#> [1] "group" "id" "value"
```

**Fix** — `out[, names(data), drop = FALSE]`. Worth a regression test:
“column order and class are preserved” is exactly the invariant a
sampling function should guarantee, for all ten.

------------------------------------------------------------------------

## F. Getting to a CRAN package (11)

### F1 · Blocker — the name collides with `sampler`, already on CRAN

[sampler
0.2.0](https://cran.r-project.org/web/packages/sampler/index.html) —
“Sample Design, Drawing & Data Analysis Using Data Frames” — has been on
CRAN since 2018, in the same problem domain. *Writing R Extensions* is
explicit:

> Because some file systems \[…\] are not case-sensitive, to maintain
> portability it is strongly recommended that case distinctions not be
> used to distinguish different packages. For example, if you have a
> package named `foo`, do not also create a package named `Foo`.

`sampleR` and `sampler` differ only in the case of the final letter.
Expect this to be rejected. Decide now — the name appears in
`DESCRIPTION`, `NAMESPACE`, the repo, every `@examples` block and the
pkgdown site.

Also unavailable: `sampling`, `samplesize`, `sampleSelection`,
`SamplingStrata`. Check candidates with
`available::available("yourname")`.

### F2 · Blocker — there is no package skeleton yet

No `DESCRIPTION`, no `NAMESPACE`, no `man/`, no `.Rbuildignore`, no
roxygen comments.

    Package: <new name>
    Title: Design-Based Sampling from Data Frames
    Version: 0.1.0
    Authors@R: person("Justin", "Chase", email = "jchase.msu@gmail.com",
                      role = c("aut", "cre"))
    Description: Draws simple random, stratified, systematic, cluster, multi-stage,
        weighted, reservoir, bootstrap, temporal and spatial samples from data
        frames, returning design weights alongside the sampled rows.
    License: GPL-3
    Encoding: UTF-8
    Depends: R (>= 4.1)
    Imports: dplyr (>= 1.1.0), rlang, lubridate, stats, utils
    Suggests: sf, testthat (>= 3.0.0), knitr, rmarkdown, withr
    Config/testthat/edition: 3
    RoxygenNote: 7.3.1

`Title` is title case with no trailing period; `Description` is a full
sentence. CRAN checks both.

### F3 · Blocker — what `R CMD check --as-cran` says today

I wrapped your two files in a minimal DESCRIPTION/NAMESPACE and ran the
real check. It reported *3 WARNINGs, 2 NOTEs*, but two are artefacts of
my sandbox — a locale warning under “syntax errors”, and a “future file
timestamps” NOTE from having no network to check the clock against. The
substantive result is **2 WARNINGs and 1 NOTE**, all three
CRAN-blocking. Unabridged:

    * checking R code for possible problems ... NOTE
    apply_max_rows: no visible global function definition for 'head'
    bootstrap_sampling: no visible global function definition for 'head'
    cluster_sampling: no visible global function definition for '%>%'
    cluster_sampling: no visible binding for global variable '.data'
    stratified_sampling: no visible global function definition for '%>%'
    stratified_sampling: no visible binding for global variable '.data'
    temporal_sampling: no visible binding for global variable '.data'
    temporal_sampling : <anonymous>: no visible binding for global variable
      '.data'
    weighted_sampling: no visible global function definition for 'sd'
    Undefined global functions or variables:
      %>% .data head sd
    Consider adding
      importFrom("stats", "sd")
      importFrom("utils", "head")
    to your NAMESPACE file.

    * checking for missing documentation entries ... WARNING
    Undocumented code objects:
      'apply_max_rows' 'bootstrap_sampling' 'cluster_sampling'
      'get_random_state' 'multi_stage_sampling' 'reservoir_sampling'
      'set_seed' 'setup_logging' 'simple_random_sampling'
      'spatial_sampling' 'stratified_sampling' 'systematic_sampling'
      'temporal_sampling' 'validate_dataframe' 'weighted_sampling'

    * checking for code which exercises the package ... WARNING
    No examples, no tests, no vignettes

**Fix**

``` r

#' @importFrom stats sd
#' @importFrom utils head
#' @importFrom rlang .data
#' @importFrom magrittr %>%     # or drop the pipe for R 4.1's native |>
```

Roxygen every exported function with `@param`, `@return` and a runnable
`@examples` block; that clears both WARNINGs at once.

### F4 · Blocker — zero tests are discoverable

testthat only collects files matching `^test`. `UAT.R` matches nothing,
and there is no `tests/testthat.R` runner — so `R CMD check` runs none
of the tests you wrote.

``` r

list.files("tests/testthat", pattern = "^test.*\\.[rR]$")
#> character(0)
```

**Fix**

``` r
# tests/testthat.R
library(testthat); library(<pkg>); test_check("<pkg>")

# split UAT.R into tests/testthat/test-<topic>.R, one per function
```

### F5 · Should fix — two of the ten test blocks fail on every run, and `context()` is retired

Eight blocks pass. **Temporal fails and spatial errors** — both
deterministically, 3 runs out of 3. Running the file top-to-bottom hides
the second one, because execution halts at the first failure; you only
see both under a reporter.

       test                                    nb failed error passed warning
     9 temporal_sampling: time window, …        3      1 FALSE      1       1
    10 spatial_sampling: region filtering, …    0      0  TRUE      0       0

    Failure ('UAT.R:188'): all(...) is not TRUE                    # cause: A9
    Error ('UAT.R:199'): Sample size cannot be larger than
      the number of points within the region.                     # cause: A16

Neither is a flaky test — both are the library telling you the truth
about A9 and A16. Fix those two and the suite goes green.

`context()` is deprecated in testthat 3e — the file name serves that
role. Remove it and set `Config/testthat/edition: 3`.

### F6 · Should fix — move `sf` to Suggests

`sf` pulls in GDAL, GEOS and PROJ. As a hard dependency, everyone
installing for stratified sampling pays for a geospatial toolchain, and
any CRAN platform where `sf` is unavailable takes your package down with
it. One function out of ten needs it.

``` r

if (!requireNamespace("sf", quietly = TRUE))
  stop("Package 'sf' is required for spatial sampling.", call. = FALSE)
```

Guard the spatial tests with `skip_if_not_installed("sf")` and examples
with `@examplesIf`.

### F7 · Should fix — the bare GPL-3 text in `LICENSE` will draw a NOTE

The file is the full 674-line GPL-3. CRAN expects `License: GPL-3` in
DESCRIPTION and reads a `LICENSE` file only when the field says
`+ file LICENSE`. Rename to `LICENSE.md`, add to `.Rbuildignore`, set
`License: GPL-3` — `usethis::use_gpl3_license()` does all three.

### F8 · Should fix — the test fixture is built from unseeded random data

`rnorm(100)`, `runif(100, -90, 90)` and `runif(100, -180, 180)` run with
no [`set.seed()`](https://rdrr.io/r/base/Random.html). Nothing currently
depends on the draw in a way that flips a result — I checked 20
independent fixtures and the outcomes were stable — but an unseeded
fixture means any future failure arrives with no way to reproduce it,
which is the worst position to be in when it is CRAN’s machine failing
and not yours.

Seed the fixture, or better, build a deterministic one under
`tests/testthat/fixtures/`. Note that `runif(100, -90, 90)` and
`runif(100, -180, 180)` can put points exactly at the poles and on the
antimeridian, where s2 is least forgiving — inset those ranges.

### F9 · Should fix — no examples, vignette, NEWS or cran-comments

Examples are checked by CRAN and are the main thing reviewers read. A
vignette — “choosing a sampling design”, walking all ten functions over
one dataset — turns this from ten utilities into a package with a
thesis. Add `R/<pkg>-package.R`, `NEWS.md`, `cran-comments.md`, and a
GitHub Actions `R-CMD-check` workflow across Linux/macOS/Windows.

### F10 · Polish — `here` is a test-only dependency that shouldn’t exist

In a package, `test_check()` loads the namespace for you. Remove
[`library(here)`](https://here.r-lib.org/), `here::i_am()` and both
[`source()`](https://rdrr.io/r/base/source.html) lines.

### F11 · Polish — housekeeping

README is three lines and still says “Testing to be finalized” — it
needs an installation block, a two-function example and a short design
table. Add `.Rbuildignore` (`^\.github$`, `^README\.Rmd$`,
`^LICENSE\.md$`, `^cran-comments\.md$`), a `_pkgdown.yml` grouping
functions by design family, and credit
[sample_py](https://github.com/elkronos/sample_py) in DESCRIPTION rather
than only the README.

------------------------------------------------------------------------

## Suggested order of work

Each step makes the next cheaper; out of order means redoing work.

1.  **Settle the name** (F1). It propagates into every other file. Run
    `available::available()` on your shortlist first.
2.  **Build the skeleton and get tests running** (F2, F4, F5, F10). You
    need a red suite before you can trust a green one.
3.  **Fix the silent wrong answers** (all of A, plus E5). Write a
    failing test for each first — the reproductions above are ready to
    paste in. A9 and A16 also turn your own suite green.
4.  **Settle the API, then document it** (E1–E3, E6, E9, then roxygen).
    Decide the `max_rows` question and the return-type rule *before*
    writing docs.
5.  **Fix the statistics and add design weights** (B1–B6). B5 changes
    what the package *is*; do it before users depend on the current
    return shape.
6.  **Validation, messages and performance** (C, D). Both mechanical
    once the API is settled.
7.  **Green `R CMD check --as-cran` on three platforms** (F3, F6–F9,
    F11), then vignette, `cran-comments.md`, submit.

------------------------------------------------------------------------

### One structural suggestion

Ten functions with ten different signatures is the shape of a script
collection. The shape of a package is a small number of verbs over a
design object —
`draw(data, design_stratified(strata = "region", n = 500))` — where
every design shares validation, seeding, weight calculation and return
type. You don’t have to go that far, but the ten functions will need to
agree on those four things anyway, and section E is mostly the cost of
them currently disagreeing.

------------------------------------------------------------------------

*Reviewed against `sampleR-main` as of 18 Feb 2025 · R 4.3.3, dplyr
1.1.4, sf 1.0-15, lubridate 1.9.3, testthat 3.2.1 · 52 findings, each
reproduced in a live session and re-derived by an independent second
pass.*
