# Bootstrap sampling

Generates bootstrap replicates. `"simple"` resamples rows independently;
`"block"` is a moving-block bootstrap, which preserves the serial
dependence in ordered data by concatenating `ceiling(n / block_length)`
independently chosen blocks per replicate.

## Usage

``` r
design_bootstrap(
  n_replicates = 1000L,
  n = NULL,
  method = c("simple", "block"),
  block_length = NULL
)
```

## Arguments

- n_replicates:

  Number of replicates to generate.

- n:

  Rows per replicate. `NULL` uses `nrow(data)`, the standard
  nonparametric bootstrap.

- method:

  `"simple"` or `"block"`.

- block_length:

  Block length for `method = "block"`, and an error for
  `method = "simple"`, where it would have no effect. `NULL` uses
  `floor(nrow(data)^(1/3))`, a common rule of thumb.

## Value

A design object, for use with
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

## Details

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) returns
all replicates in one data frame with a leading `.replicate` column.
Split them with `split(out, out$.replicate)`.

## Not a survey design

The bootstrap resamples the data you already have to approximate the
sampling distribution of a statistic (Efron 1979). It is not a
probability sample of a finite population, so it has no inclusion
probabilities and cannot be passed to
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md).
For the design-based variance of a total or mean, use the estimators
there; for a statistic they do not cover, such as a median, this is the
tool. Blocks preserve short-range dependence in ordered data (Künsch
1989).

## References

Efron, B. (1979). Bootstrap methods: another look at the jackknife.
*Annals of Statistics*, 7, 1–26.

Künsch, H. R. (1989). The jackknife and the bootstrap for general
stationary observations. *Annals of Statistics*, 17, 1217–1241.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)

Other designs:
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md),
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md),
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md),
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
[`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md),
[`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md),
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)

## Examples

``` r
df <- data.frame(id = 1:100, value = (1:100) / 10)

reps <- draw(df, design_bootstrap(n_replicates = 5, n = 20), seed = 1)
table(reps$.replicate)
#> 
#>  1  2  3  4  5 
#> 20 20 20 20 20 

# A statistic per replicate
vapply(split(reps, reps$.replicate), function(r) mean(r$value), numeric(1))
#>     1     2     3     4     5 
#> 5.450 5.285 4.395 4.735 5.570 
```
