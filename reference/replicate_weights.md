# Replicate weights for a sample

Builds replicate weights from a drawn sample for
[`survey::svrepdesign()`](https://rdrr.io/pkg/survey/man/svrepdesign.html),
from the design the sample was drawn with. Each replicate is the
sample's `.weight` column times a factor, and the variance of a total
over the replicates is that of one of
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)'s
estimators:

## Usage

``` r
replicate_weights(sample, type = c("JKn", "JK1"))
```

## Arguments

- sample:

  A data frame returned by
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with
  `weights = TRUE`.

- type:

  The form of the jackknife. `"JKn"`, the default, deletes units within
  strata and keeps each stratum's correction in `rscales`. `"JK1"`
  deletes across the whole sample, with the correction in `scale`, and
  is given only where the jackknife has one stratum; with more, it would
  count the differences between strata as sampling variance, and is
  refused. A systematic sample with one start, and strata walked in
  order, get successive-difference replicates of type `"other"`
  whichever is given.

## Value

A numeric matrix with a row for each row of `sample`, in its order, and
a column for each replicate, holding replicate weights. Its attributes
`type`, `scale` and `rscales` are the arguments of the same names of
[`survey::svrepdesign()`](https://rdrr.io/pkg/survey/man/svrepdesign.html),
which with `combined.weights = TRUE` computes the variance as
`scale * sum(rscales * (t_r - t)^2)` over the replicate estimates `t_r`,
with `t` their mean. For a total that mean is the sample's estimate.

## Details

- For most designs, the delete-a-unit jackknife of
  `ht_total(variance = "jackknife")`. Each replicate deletes one primary
  sampling unit, a cluster where the design has them and a row
  otherwise, and reweights the other units of its stratum; strata are
  the design's strata, the intervals of
  [`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md)
  and the strata of
  [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
  and
  [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md).
  `rscales` holds `(m - 1) / m` for a stratum of `m` sampled units,
  times its finite population correction. Units taken with certainty
  keep their weight in every replicate, as the jackknife holds them
  fixed. Where the default of
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  is another estimator, such as Deville's approximation for
  [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
  with `method = "systematic"`, the local mean approximation for
  [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
  or the two-stage variance for
  [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
  and balanced
  [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
  the replicates give the jackknife's variance, not that default.

- For
  [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md)
  with one start, and
  [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md)
  with `within`, successive-difference replicates in the order the
  design walked, whose variance is the successive-difference
  approximation
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  gives by default, with its term for a random sample size. Each pair of
  neighbours in the walk shares a row of a Hadamard matrix, which raises
  one of the two and lowers the other. The first and last rows of a walk
  are not paired with each other: the circular form of
  successive-difference replication pairs them, which adds the squared
  difference between the two ends of the sort order. On a frame of 240
  rows sorted on a size measure related to `y`, at an interval of 8,
  that more than doubled the variance in each of the 8 possible samples.

- For
  [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md)
  with `n_starts` of 2 or more, one replicate per start, each deleting
  the rows its start drew. The variance is the random group variance
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  gives by default.

Rows above the threshold of
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md)
keep their weight in every replicate, and the replicates for the rows
below it are those of `rest`. A sample whose every row was taken with
certainty gets one replicate equal to its weights, and a variance of 0.

Poisson sampling, alone or as the `rest` of
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md),
is refused: its sample size is random, and deleting rows from it
understates the variance.
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
hands it to `survey`, which models it. A sample the jackknife declines,
such as one with a single unit in a stratum, is refused with the
jackknife's reason.

Every replicate is a multiple of `.weight`, so an adjustment written
into `.weight` after the draw, such as trimming, is carried into each
one unchanged. One that depends on which rows responded is more
accurately made again within each replicate.

The matrix of weights takes `8 * n * R` bytes for `n` sampled rows and
`R` replicates, and building it about twice that. The jackknife has a
replicate for each sampled unit, and successive differences have the
first power of 2 above the number of pairs of neighbours in the walks,
counting one more for each walk of random size, so the size grows with
the square of the sample. Past about 1 GB, `replicate_weights()`
refuses: from about 11,200 rows for a jackknife over rows, and about
8,200 for successive differences. Such a sample can still be estimated
with
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
or handed to `survey` with
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md).

## See also

[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md),
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)

## Examples

``` r
pop <- data.frame(
  id = 1:400,
  site = rep(c("a", "b", "c", "d"), times = c(200, 100, 60, 40)),
  spend = 10 + (1:400 * 37) %% 490
)
s <- draw(pop, design_stratified("site", n = 60), seed = 1, weights = TRUE)

rw <- replicate_weights(s)
des <- survey::svrepdesign(
  data = s, repweights = rw, weights = ~.weight,
  type = attr(rw, "type"), scale = attr(rw, "scale"),
  rscales = attr(rw, "rscales"), combined.weights = TRUE
)
survey::svytotal(~spend, des)
#>       total     SE
#> spend 93473 7191.9

# The same standard error as the jackknife here
ht_total(s, "spend", variance = "jackknife")
#> Horvitz-Thompson total  (stratified design, n = 60)
#>   estimate 93,473.33
#>   se       7,191.859  (jackknife)
#>   95% CI  79,066.31 to 107,880.4  (t, 56 df)
#>   deff     1.04  (about the same as simple random sampling)
#> 
#>   Jackknife over 60 primary sampling unit(s) in 4 strata, with a finite
#>   population correction.

# A systematic sample on a sorted frame: the replicates follow the walk
sys <- draw(pop, design_systematic(interval = 8, order_by = "spend"),
            seed = 1, weights = TRUE)
rw <- replicate_weights(sys)
des <- survey::svrepdesign(
  data = sys, repweights = rw, weights = ~.weight,
  type = attr(rw, "type"), scale = attr(rw, "scale"),
  rscales = attr(rw, "rscales"), combined.weights = TRUE
)
survey::SE(survey::svytotal(~spend, des))
#> [1] 367.5401
ht_total(sys, "spend")$se
#> [1] 367.5401
```
