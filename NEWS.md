# drawn 0.1.0

First packaged release. The code was previously two loose script files under the
name `sampleR`, exporting ten `*_sampling()` functions.

## Continuous integration and a vignette

* `R-CMD-check` runs on macOS, Windows and Linux, against the current R release,
  the previous one, and R-devel, failing on warnings as well as errors. A
  `pkgdown` workflow publishes the reference site.
* A `Choosing a sampling design` vignette works one audit population through all
  ten designs, Horvitz-Thompson estimation, and the cases with no closed-form
  inclusion probability.
* `tools/check-duplicates.R` guards against a function being defined in two
  files under `R/`. R sources that directory alphabetically, so a copy left
  behind by a refactor silently wins over the intended one -- which is exactly
  what happened to `allocate()` and is how the bug below survived a green suite.

## Allocation with a floor

`min_per_stratum` could return more rows than requested. The reduction branch of
the allocator made a single pass, so it could take back at most one row per
stratum; four strata with a floor of 8 returned 62 rows for a request of 60. It
now cycles until the total is exact, and a test covers every feasible
combination of `n` and `min_per_stratum`.

## Key columns must be groupable

A data frame column can itself be a list or a matrix, and neither can serve as a
grouping or ordering key. Both failure modes are now refused by name:

* A **matrix column** used to fail silently. `design_cluster()` on one returned a
  single row from a twenty-row frame, because `unique()` on a matrix works
  element-wise; `design_stratified()` produced a meaningless stratification.
* A **list column** used to fail from inside base R with "unimplemented type
  'list' in 'orderVector1'".

`strata`, `clusters` and `order_by` now check the column and say which one is at
fault. Factors, character, numeric and Date keys are unaffected.

## Performance

The rewrite made the two hot spots vectorised rather than row-at-a-time:

* Reservoir sampling from a data frame: **6.5 s to 0.002 s** at 200,000 rows,
  because a data frame is not a stream and no longer pretends to be. A genuine
  stream of 200,000 items now runs Algorithm L in 0.24 s.
* Temporal sampling: **3.11 s to 0.036 s** at 50,000 rows across 833 buckets,
  assigning buckets arithmetically instead of re-filtering per interval.

At one million rows, every design draws in under a quarter of a second.

## Inclusion probabilities and design weights

New `inclusion_prob(data, design)` and `design_weight(data, design)` return the
first-order inclusion probability of every population row under a design, and
its reciprocal. `draw(..., weights = TRUE)` prepends them to a sample as `.prob`
and `.weight`, so `sum(y * .weight)` is an unbiased Horvitz-Thompson total.

Nine designs have exact closed forms, verified against 20,000-draw simulations.
Four do not, and the package reports that rather than approximating:

* `design_cluster(balanced = TRUE)` — the per-cluster take is the smallest
  *selected* cluster's size, which is random. Clusters of 2/4/6/8 give realised
  rates of 0.50, 0.41, 0.34, 0.25, against the flat 0.50 a naive
  `n_clusters / N_clusters` would claim.
* `design_multistage(allocation = "proportional")` — stage-two allocation
  depends on which clusters were drawn. Clusters of 3/5/7/9 give 0.17, 0.20,
  0.17, 0.15, against a naive 0.33, 0.20, 0.14, 0.11.
* `design_weighted(method = "successive")` — successive sampling has no closed
  form.
* `design_bootstrap()` — not a probability sample of a finite population.

`simulate = TRUE` estimates any of them by Monte Carlo.

## Probability-proportional-to-size sampling

`design_weighted()` gains `method`. The default `"successive"` is unchanged and
is what `base::sample(prob=)` does. Two new methods give inclusion probabilities
that really are proportional to size:

* `"systematic"` — systematic PPS with a random start. `pi_i = n * p_i` exactly,
  fixed sample size. Units heavy enough that `n * p_i > 1` are taken with
  certainty and the remainder rescaled, iterating until every probability is
  valid.
* `"poisson"` — independent inclusion at `pi_i = n * p_i`. Exactly proportional,
  random sample size averaging `n`.

Horvitz-Thompson totals under both were confirmed unbiased by simulation, as
were those under the simple, stratified, cluster, equal-allocation multi-stage
and systematic designs.

## Renamed

* The package is now `drawn`. `sampleR` differs from CRAN's `sampler` only
  in the case of its final letter, which *Writing R Extensions* asks authors to
  avoid.

## The API is now a design object and one verb

The ten `*_sampling()` functions are gone. Build a design, then draw from it:

```r
# was
stratified_sampling(df, strata_column = "site", sample_size = 100, seed = 1)

# now
draw(df, design_stratified("site", n = 100), seed = 1)
```

A design is an ordinary value: build it once, print it, store it, reuse it
across data sets. Every design shares one contract for validation, seeding,
missing keys and return shape, which is what the ten separate signatures could
not.

Argument names were settled at the same time:

| Was | Now |
|---|---|
| `sample_size`, `stage_two_sample_size` | `n` — always a total, never per-group |
| `sample_size` in `temporal_sampling()` | `per_interval` — named for what it is |
| `equal_samples`, `proportional_stage_two` | `allocation = c("proportional", "equal")` |
| `strata_column`, `cluster_column`, `weights_column`, `time_column` | `strata`, `clusters`, `weights`, `time` |
| `sort_column` | `order_by` |
| `start_time`, `end_time` | `from`, `to` |
| `time_zone` | `tz` |
| `num_clusters`, `num_samples` | `n_clusters`, `n_replicates` |
| `latitude_column`, `longitude_column` | `coords = c(x, y)` — longitude first, matching `sf` |
| `data_stream` | `data` |

Removed:

* **`max_rows`** everywhere. In simple random sampling it was provably identical
  to asking for fewer rows; everywhere else it meant "draw a proper sample, then
  throw part of it away", which is never what a design should do. Every design
  already has a size argument.
* **`normalization`** from weighted sampling. `prob=` normalises by the sum, so
  rescaling never changed anything a user could want; both offered methods
  subtracted the minimum, which silently re-weighted every row.
* **`complex_region`** from spatial sampling. A list of geometries is now unioned
  because it is a list, with no flag to set.
* **`set_seed()`, `get_random_state()`, `setup_logging()`**. Placeholders: the
  first two wrapped one line of base R, the third did nothing.

Also new: `n` is the total across selected clusters in `design_multistage()`
under both allocations, matching `design_stratified()`. It previously meant
per-cluster under one setting and total under the other.

## Correctness fixes

Each of these returned a plausible-looking result that was wrong.

* Cluster and multi-stage designs no longer mistake a single numeric cluster id
  for a range. `sample(x, n)` reinterprets a length-one numeric `x` as
  `seq_len(x)`, so a frame whose only cluster id was `5` silently returned no
  rows at all.
* Capping a result no longer discards whole strata. `max_rows` was `head()`, and
  stratified, cluster and temporal results arrive sorted by group, so a cap kept
  only the first group. `n` is now allocated across groups rather than applied
  afterwards.
* Multi-stage proportional allocation divides by the rows in the *selected*
  clusters rather than the whole population, and no longer floors every cluster
  at one row. The old behaviour both under- and over-delivered: 10 rows for a
  request of 20, and 20 for a request of 5.
* Weighted sampling returns the caller's weights unmodified. It used to hand back
  the internally rescaled values.
* No rescaling can drive a row's probability to zero. `min-max` and `z-score`
  both subtracted the minimum and added `.Machine$double.eps`, leaving the
  lightest row at ~2.2e-16 and effectively unselectable.
* Balanced cluster sampling works when the cluster column is a factor with
  unused levels. It used to return zero rows.
* `design_bootstrap(method = "block")` is a moving-block bootstrap. It used to
  pick one block and resample inside it, destroying the serial dependence the
  method exists to preserve. `block_length` is a separate argument from `n`.
* Temporal sampling uses `POSIXct` and `Date` columns as they are. Running
  `ymd_hms()` over an already-parsed column dropped midnight values to `NA`.
* A region given as a list of geometries is unioned whatever its length.
  `do.call(st_union, region)` mapped the third and later elements onto
  `st_union()`'s own arguments.
* Count arguments reject negative, non-numeric, fractional, `NA` and length-two
  values at construction time. `head(x, -3)` silently trimmed from the end, and a
  length-two value silently dropped columns.
* Single-column data frames survive every design intact.
* Reservoir sampling trims its result instead of padding with `NULL`s.
* An empty time window returns a zero-row frame with the input's columns, not a
  0x0 frame.
* Multi-stage matches cluster labels with `%in%`, not `==`. Because `NA == "a"`
  is `NA` and indexing with `NA` yields an all-`NA` row, missing labels used to
  manufacture rows that were never in the data — and inflated the row count
  enough to defeat the `replace = FALSE` guard.
* Spatial sampling warns when a region edge spans more than 180 degrees of
  longitude, which under spherical geometry is drawn the short way across the
  antimeridian and collapses a "whole world" rectangle into a 2-degree strip.

## Statistical accuracy

* Allocation uses the largest-remainder method and hits `n` exactly. Independent
  per-group rounding used to drift, and banker's rounding took four strata of 25
  down to 8 rows for a requested 10.
* `min_per_stratum` and `min_per_cluster` (default `0`) make the rare-group floor
  opt-in rather than mandatory.
* `unit = "months"` and `"years"` step by calendar boundaries. They previously
  used fixed 30.4375-day durations, so buckets drifted off month boundaries
  immediately. `"seconds"`, `"minutes"` and `"years"` are new.
* `tz` is applied when the column *and* both window bounds are parsed, so it
  actually shifts the window. It previously relabelled the output only.
* `?design_weighted` documents that sampling without replacement gives
  successive-sampling inclusion probabilities, not πPS.

## Performance

* Reservoir sampling takes a vectorised path for data frames and Algorithm L for
  real streams. It used to `split()` the entire frame into one-row data frames
  first: about 7 seconds for 200,000 rows against under a millisecond.
* Temporal sampling assigns buckets with `findInterval()` in one pass instead of
  re-filtering the whole frame per interval: 0.04s against 3.1s for 50,000 rows
  across 833 buckets.

## Packaging

* Proper package: `DESCRIPTION`, roxygen-generated `NAMESPACE` and `man/`,
  runnable examples throughout.
* `sf` moved to `Suggests` and required at run time, so installing the package no
  longer pulls in GDAL, GEOS and PROJ.
* `dplyr` and `rlang` are no longer needed; the internals are base R. `lubridate`
  is the only non-base import.
* Tests are discoverable: `tests/testthat.R` plus `test-*.R` files, testthat 3rd
  edition, with one regression test per correctness fix above. Each was confirmed
  to fail against the original implementation before its fix landed.
* `R CMD check --as-cran` is clean.
