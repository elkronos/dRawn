# Hand a sample to the survey package

Builds a
[`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html)
object from a drawn sample, so the analysis this package does not do —
subpopulation estimates, regression, calibration, quantiles with proper
standard errors — can be done by the package that does.

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
one uses the design's joint inclusion probabilities; `survey`
reconstructs the variance from the design's *shape*. So the job here is
to express each design in `survey`'s own terms rather than hand over a
weight column and hope.

## What maps to what

|  |  |  |
|----|----|----|
| **Design** | **Expressed as** | **Standard errors** |
| [`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md), [`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md), [`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md) | `ids = ~1` with `fpc` the frame size | identical |
| [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md) | `strata` from the strata columns, `fpc` each stratum's size | identical |
| [`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md) | `strata` from the sampling intervals, `fpc` each interval's size | identical |
| [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md) | `ids` the cluster column, `fpc` the number of clusters | identical |
| [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md), `"systematic"` | `ids = ~1` with `fpc` the frame size | identical |
| [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md), `"poisson"` | [`survey::poisson_sampling()`](https://rdrr.io/pkg/survey/man/poisson_sampling.html), which models the random size | identical |
| [`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md) | the certainty rows as their own stratum, taken whole | identical |
| [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md) | `ids` the cluster column | **differ by a few percent** |
| [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md) | `ids = ~1` with `fpc` the frame size | **`survey` returns one; this package declines** |

"Identical" means to floating point, and is checked by this package's
tests against
[`survey::svytotal()`](https://rdrr.io/pkg/survey/man/surveysummary.html).
The two exceptions are real and worth knowing:

- **Multistage.** `survey` uses the ultimate-cluster approximation,
  which attributes all the variance to the first stage and ignores
  sampling within clusters.
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  uses the exact two-stage form. `survey`'s is the smaller of the two,
  by around 5–10% on a typical frame.

- **Systematic.** Most pairs of rows can never co-occur, so no
  design-unbiased variance exists and
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  returns `NA` with a note. `survey`, having only been told `ids = ~1`,
  computes the simple-random variance — which is the conservative
  substitute, not the design's own.

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
population for
[`svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html) to
represent.

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
#>   95% CI  95,086.78 to 121,579.9
#>   deff     1.02  (about the same as simple random sampling)

# Now the analysis this package does not do
survey::svyby(~spend, ~site, des, survey::svymean)
#>   site    spend       se
#> a    a 286.0000 21.97916
#> b    b 247.1333 35.82923
#> c    c 297.6667 50.04617
#> d    d 214.0000 53.03870
```
