# Describe a drawn sample

What you actually got, against what the design asked for. Worth a look
before analysing: it surfaces strata that came up short, weights that
vary more than you expected, and rows the design could never have
reached.

## Usage

``` r
sample_summary(sample)
```

## Arguments

- sample:

  A data frame returned by
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with
  `weights = TRUE`.

## Value

A list with a [`print()`](https://rdrr.io/r/base/print.html) method,
holding:

- `design`:

  The design's type, as a string.

- `n`, `N`:

  Rows drawn, and rows in the frame.

- `N_reachable`:

  Rows of the frame the design gives a probability above 0, the
  population it samples. For weights conditional on the clusters
  selected, the rows of the clusters not selected count, since another
  selection could have drawn them. The printed sampling fraction is
  `n / N_reachable`. `NA` only when the frame's probabilities cannot be
  worked out, for an empty sample saved by drawn 0.2.0, which records
  none, from a design with no closed-form inclusion probability. The
  fraction is then `n / N`.

- `expected_n`:

  The sample size the design expects, the sum of the inclusion
  probabilities over the frame. It is `n` for a design of fixed size,
  and printed beside `n` when the two differ, as they can under Poisson
  sampling, a
  [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md)
  `interval` that does not divide the frame, or clusters of unequal
  size. `NA` where `N_reachable` is, and for weights conditional on the
  clusters selected.

- `weight_range`:

  The smallest and largest design weight.

- `weight_cv`:

  Their coefficient of variation. Large values mean a few rows carry
  most of the estimate.

- `n_kish`:

  Kish's effective sample size, `sum(w)^2 / sum(w^2)` over the weights
  `w`, which is `n / (1 + cv^2)` with the coefficient of variation
  computed with divisor `n`: the size of a simple random sample that
  would estimate a mean as precisely if unequal weights were the only
  cost. Clustering and stratification are not in it;
  [`deff()`](https://elkronos.github.io/dRawn/reference/deff.md)
  measures those.

- `n_certain`:

  Sampled rows taken with certainty, at a probability of 1.

- `certain_by`:

  Those rows by why they are certain, a named count: `threshold`, above
  the threshold of
  [`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md);
  `size`, put at probability 1 by a size measure, in
  [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
  or in
  [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
  or
  [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
  with `size`, which includes the rows of a cluster put at 1 by its size
  and then taken whole; and `whole`, in a stratum, interval or cluster
  taken whole, or a census, with a size measure or without. Under a size
  measure a group counts as taken whole when the design gives every row
  of it probability 1, or every cluster of a multistage design a
  first-stage probability of 1, so a Poisson sample that holds every row
  of a stratum by chance files its certain rows under `size`.

- `conditional`:

  The number of clusters selected, when the weights are conditional on
  them, as for
  [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)
  with `balanced = TRUE` (see "Design weights" in
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md));
  otherwise `NA`.

- `unreachable`:

  Frame rows with inclusion probability 0, which the design could never
  have selected, such as rows `na_rm` drops. `NA` where `N_reachable`
  is.

- `by_group`:

  A data frame of `group`, `drawn`, `in_frame` and `rate` per stratum,
  interval or cluster, or `NULL` for a design with no grouping. The
  strata are those of
  [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
  and of
  [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
  and
  [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
  when they have `strata`, and every stratum of the frame is listed, one
  allocated no rows included. A stratum's `in_frame` counts the rows the
  design can draw from it, as `N` of
  [`allocation_table()`](https://elkronos.github.io/dRawn/reference/allocation_table.md)
  does, so rows `na_rm` drops are left out, and a stratum it empties is
  not listed. The intervals of
  [`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md)
  are listed in time order and named by their start, such as
  `"2024-01-03 00:00"`. For
  [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)
  and
  [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
  only the clusters drawn are listed. For
  [`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md),
  the groups are those of its `rest` over the rows below the threshold,
  after a `"(certain)"` row for the rows above it.

- `group_col`:

  The column(s) `by_group` is keyed on.

- `drawn`:

  How the sample was drawn, from its `drawn_meta` attribute: `seed`,
  `rng_kind`, `version`, `time`, `frame_rows` and `made_by` (see
  "Seeding" in
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)).
  `NULL` for a sample from an earlier version of drawn, which recorded
  none. The print shows the seed, the version and the date.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md),
[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md)

## Examples

``` r
set.seed(1)
pop <- data.frame(
  id = 1:400,
  site = rep(c("a", "b", "c", "d"), times = c(200, 100, 60, 40))
)
s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
sample_summary(s)
#> Sample of 40 from 400  (stratified design)
#>   sampling fraction  0.1
#>   design weights     10 to 10   (cv 0)
#>   effective n        40  (Kish, from the weights alone)
#>   drawn              with seed 1, by drawn 0.2.0.9000 on 2026-11-17
#> 
#>   by site:
#>     group  drawn  in frame   rate
#>     a         20      200  0.100
#>     b         10      100  0.100
#>     c          6       60  0.100
#>     d          4       40  0.100
```
