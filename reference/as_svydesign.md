# Hand a sample to the survey package

Builds a
[`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html)
object from a drawn sample, so the analysis this package does not do —
regression, calibration, quantiles with proper standard errors,
replicate weights — can be done by the package that does.

## Usage

``` r
as_svydesign(sample, ...)
```

## Arguments

- sample:

  A data frame returned by
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with
  `weights = TRUE`.

- ...:

  Passed to
  [`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html).
  Anything named here overrides what the mapping above would have
  supplied, so `nest = TRUE` or a replacement `fpc` is yours to set.

## Value

A `survey.design` object.

## Details

The two packages compute variance from different starting points. This
one uses the design's inclusion probabilities; `survey` reconstructs the
variance from the design's *shape*. So the job here is to express each
design in `survey`'s own terms rather than hand over a weight column and
hope.

## What maps to what

|  |  |  |
|----|----|----|
| **Design** | **Expressed as** | **Standard errors** |
| [`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md), [`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md) | `ids = ~1` with `fpc` the rows the design can reach | identical |
| [`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md) | `ids = ~1` with `fpc` the rows the stream reached | identical |
| [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md) | `strata` from the strata columns, `fpc` each stratum's size | identical |
| [`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md) | `strata` from the sampling intervals, `fpc` each interval's size | identical |
| [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md) | `ids` the cluster column, `fpc` the number of clusters | identical |
| [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md) | two stages, `ids = ~cluster + row`, `fpc` the number of clusters and each cluster's size | identical |
| [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md), `"poisson"` | [`survey::poisson_sampling()`](https://rdrr.io/pkg/survey/man/poisson_sampling.html), which models the random size | identical |
| [`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md) | the certainty rows as their own stratum, taken whole | as `rest` |
| [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md), `"systematic"` | `pps = "brewer"` with the inclusion probabilities | within about 0.2 percent |
| [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md) | as its first-order design: `ids = ~1`, or `pps = "brewer"` with `size` | **`survey`'s is larger** |
| [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md) | `ids = ~1` with `fpc` the frame size | **differ** |

"Identical" means to floating point, and is checked by this package's
tests against
[`survey::svytotal()`](https://rdrr.io/pkg/survey/man/surveysummary.html),
as are domain estimates against
[`survey::svyby()`](https://rdrr.io/pkg/survey/man/svyby.html) and
degrees of freedom against
[`survey::degf()`](https://rdrr.io/pkg/survey/man/svychisq.html). The
exceptions are real and worth knowing:

- **Systematic PPS.** Neither package has its joint probabilities. This
  one uses Deville's approximation and `survey` uses Brewer's; the two
  are close relatives and agree to a fraction of a percent.

- **Spread.** `survey` has no estimator for a spatially balanced sample,
  so it is told only the first-order design and computes a variance that
  ignores the spreading. That errs conservative, often substantially.

- **Systematic.** `survey`, told only `ids = ~1`, computes the
  simple-random variance;
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  uses the successive-difference approximation, which can see a trend
  along the sort order. Neither is design-unbiased, because no such
  estimator exists for one systematic sample.

Certainty rows arrive in a stratum where `n == N`, so `survey`'s own
finite population correction zeroes them out, matching this package's
treatment.

Two compositions have no single `survey` design and are refused rather
than approximated:
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md)
over a cluster, multistage or Poisson `rest`, where the certainty rows
and the rest are different kinds of sampling unit.
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md)
is refused outright — it resamples the sample, so there is no finite
population for `svydesign()` to represent.

A stratum or interval with a single sampled row makes `survey` stop with
"Stratum has only one PSU" unless `options(survey.lonely.psu)` says
otherwise;
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
declines in the same situation. Draw at least two per stratum.

## See also

[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)

## Examples

``` r
set.seed(1)
pop <- data.frame(
  id = 1:400,
  site = rep(c("a", "b", "c", "d"), times = c(200, 100, 60, 40)),
  spend = round(stats::runif(400, 10, 500))
)
s <- draw(pop, design_stratified("site", n = 60), seed = 1, weights = TRUE)

des <- as_svydesign(s)
survey::svytotal(~spend, des)
#>        total     SE
#> spend 108333 6758.6

# The same total, and the same standard error
ht_total(s, "spend")
#> Horvitz-Thompson total  (stratified design, n = 60)
#>   estimate 108,333.3
#>   se       6,758.57  (analytic)
#>   95% CI  94,794.29 to 121,872.4  (t, 56 df)
#>   deff     1.02  (about the same as simple random sampling)

# Now the analysis this package does not do
survey::svyquantile(~spend, des, quantiles = 0.5)
#> $spend
#>     quantile ci.2.5 ci.97.5       se
#> 0.5      252    198     341 35.69217
#> 
#> attr(,"hasci")
#> [1] TRUE
#> attr(,"class")
#> [1] "newsvyquantile"
```
