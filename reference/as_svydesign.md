# Hand a sample to the survey package

Builds a
[`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html)
object from a drawn sample, so that the survey package can do the
analysis this package does not, such as regression, calibration and
analysis with replicate weights.

## Usage

``` r
as_svydesign(sample, ...)
```

## Arguments

- sample:

  A data frame returned by
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with
  `weights = TRUE`. See "Weights and design columns".

- ...:

  Passed to
  [`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html)
  by name, and an unnamed argument is refused. Anything named here
  overrides what the mapping above would have supplied, so `nest = TRUE`
  or a replacement `fpc` is yours to set.

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
| [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md) | `ids` a `.cluster` column numbering the clusters, `fpc` the number of clusters | identical |
| [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md), and [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md) with `balanced = TRUE` | two stages, `ids` the `.cluster` numbers then the rows, `fpc` the number of clusters and each cluster's size; with one sampled row from a selected cluster of several, one stage, `ids` the `.cluster` numbers and no `fpc` | identical |
| [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md), `"poisson"` | [`survey::poisson_sampling()`](https://rdrr.io/pkg/survey/man/poisson_sampling.html), which models the random size; with `strata`, the strata too, which change only the degrees of freedom | identical |
| [`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md) | the certainty rows as one unit in a stratum of their own, taken whole, beside `rest` as it maps alone; over a Poisson `rest`, [`survey::poisson_sampling()`](https://rdrr.io/pkg/survey/man/poisson_sampling.html) from the whole sample's `.prob` | as `rest` |
| [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md), `"systematic"` | `pps = "brewer"` with the inclusion probabilities; with `strata`, `strata` from the strata columns and Brewer's within them | close for totals, less so for means |
| [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md) with `size` | two stages, `ids` the `.cluster` numbers then the rows, `pps = "brewer"` with each cluster's stage-one probability and `m / N_i` at the second stage; each cluster taken with certainty a stratum of its own, whose rows are its units; with one sampled row from a selected cluster of several, one stage, `ids` the `.cluster` numbers and no `fpc` | close; identical in one stage |
| [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md) | as its first-order design: `ids = ~1`, or `pps = "brewer"` with `size`; with `strata`, `strata` from the strata columns and `fpc` each stratum's size, or Brewer's within them | **differ** |
| [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md) with `n_starts` of 2 or more | `ids` a `.cluster` column numbering the replicates, `fpc` the `n_starts * interval` possible starts | identical |
| [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md) | `ids = ~1` with `fpc` the rows the walk can reach | **differ** |
| [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md) with `within` | `strata` from the strata columns, `fpc` each stratum's rows walked | **differ** |

"Identical" means to floating point under `survey`'s default options,
and is checked by this package's tests against
[`survey::svytotal()`](https://rdrr.io/pkg/survey/man/surveysummary.html),
as are domain estimates against
[`survey::svyby()`](https://rdrr.io/pkg/survey/man/svyby.html) and
degrees of freedom against
[`survey::degf()`](https://rdrr.io/pkg/survey/man/svychisq.html). The
exceptions are real and worth knowing:

- **Systematic PPS.** Neither package has its joint probabilities. This
  one uses Deville's approximation and `survey` uses Brewer's. In
  simulations on frames of 100 to 300 rows with gamma and lognormal
  sizes, the standard errors of totals usually agreed within 0.5
  percent, and within 3 percent in every sample. Those of means differed
  more, and more as the sampling fraction grew: typically by under 2
  percent at a fraction of a tenth or a fifth, 2 to 4 percent at a third
  and 3 to 7 percent at a half, where single samples differed by up to
  12 percent. On a frame whose units below certainty had probabilities
  from 0.05 to 0.95, with about 17 of them in each sample, totals
  differed by about 1 percent typically and by up to 9 percent, and
  means by 7 to 9 percent typically and by up to 17 percent. The degrees
  of freedom agree. With `strata`, each stratum is handed over the same
  way.
  [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
  with `size` meets the same two approximations at the first stage,
  while the variance within the selected clusters is the same in both
  packages. Where the clusters below certainty had stage-one
  probabilities of 0.12 to 0.41, the standard errors of totals agreed
  within 0.3 percent over 200 samples. With a cluster at 0.93 among
  them, `survey`'s came out 0.97 to 1.48 times this package's: Brewer's
  approximation centres on the plain mean of the clusters' totals over
  their probabilities, which that cluster pulls, where Deville's weights
  each by `1 - pi`.

- **Spread.** `survey` has no estimator for a spatially balanced sample,
  so it is told only the first-order design and computes a variance that
  ignores the spreading. Relative to the true variance that errs
  conservative, often substantially. Relative to
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)'s
  local mean approximation it goes either way. In simulations drawing 30
  of 300 points spread over a square, `survey`'s standard error was the
  larger in every sample when `y` varied smoothly over the square, three
  times as large at the median. When `y` had no spatial pattern it was
  the smaller in 67 to 79 percent of samples, by 4 to 6 percent at the
  median, and drawing 150 of the 300 it was the smaller in every sample,
  by 30 percent at the median, because the local mean has no finite
  population correction.

- **Systematic.** `survey`, told only `ids = ~1`, computes the
  simple-random variance;
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  uses the successive-difference approximation, which can see a trend
  along the sort order. Neither is design-unbiased, and at a step of 2
  rows or more no such estimator exists for one systematic sample.
  Within strata, `survey` computes the stratified simple-random variance
  and
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  adds up the strata's successive-difference approximations. With
  `n_starts` of 2 or more the replicates are clusters and the two agree.
  A replicated sample with a start that drew no rows, which happens only
  when `n_starts * interval` exceeds the rows walked, is refused,
  because `survey` would count one replicate fewer.

Units taken with certainty are in every sample and add nothing to the
variance: the rows above a
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md)
threshold, strata and intervals taken whole, the units
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
and
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
with `size` put at probability 1, and the clusters
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
and
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)
with `balanced = TRUE` take whole when they select every cluster. Each
group of them goes to `survey` as one primary sampling unit in a stratum
of its own, with an `fpc` of 1, so `survey` gives them no variance and
[`survey::degf()`](https://rdrr.io/pkg/survey/man/svychisq.html) leaves
them out of the degrees of freedom, as
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
does. For the PPS designs this also keeps them out of Brewer's
approximation, which left among the other units inflated `survey`'s
standard error, by a factor of 23 on a frame with four dominant units.
Beside a `rest` of clusters, or of the replicates of several systematic
starts, the certainty rows form one more primary sampling unit. In two
stages, such a unit is taken whole at the second stage too.

A sample in which every row was taken with certainty, a census, goes to
`survey` as `ids = ~1` with `fpc` the number of rows. Both packages then
report a variance of 0, and
[`survey::degf()`](https://rdrr.io/pkg/survey/man/svychisq.html) gives
`n - 1` where
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
gives `Inf`. A sample with no rows is refused, and so is a certainty
sample whose `rest` drew none of the rows below the threshold, whose
total `survey` would report as exact where
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
reports no variance.

[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md)
over a Poisson `rest` goes to `survey` as
[`survey::poisson_sampling()`](https://rdrr.io/pkg/survey/man/poisson_sampling.html)
built from the whole sample's `.prob`, so the certainty rows enter at
probability 1 and add no variance.
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md)
is refused outright — it resamples the sample, so there is no finite
population for `svydesign()` to represent.

A stratum or interval with a single sampled row makes `survey` stop with
"Stratum has only one PSU" unless `options(survey.lonely.psu)` says
otherwise;
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
declines in the same situation. Draw at least two per stratum. Under
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
with `size`, a cluster taken with certainty is a stratum of its sampled
rows, so one that gives a single row of several is such a stratum.
`options(survey.ultimate.cluster = TRUE)` makes `survey` leave out the
second stage of every two-stage design, so its standard errors then
differ from
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)'s.
The agreement the table states holds under the default, `FALSE`.

A sample weighted conditional on the clusters selected (see "Design
weights" in
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)) goes
over in two stages too. `survey` works out each cluster's second-stage
sampling fraction from the rows sampled in it, which is the take given
that selection, so its estimate and standard error are those of
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md).

One sampled row from a selected cluster of several leaves the variation
within clusters unmeasured.
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
then falls back to the delete-a-cluster jackknife, which is the
ultimate-cluster variance with no finite population correction, and
`survey` gets that design: one stage, `ids` the `.cluster` numbers and
no `fpc`. So does a certainty design's `rest` of balanced clusters,
which [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)
weights without regard to the clusters selected, leaving no two-stage
estimator. Under
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
with `size` the jackknife runs over the clusters below certainty, and a
cluster taken whole is held fixed in a stratum of its own in both
packages.

## Weights and design columns

`survey` gets the sample's `.weight` column as its weights. Weights
adjusted after the draw, for nonresponse or by trimming, and written
into `.weight` are therefore estimated here, while
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
and [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
keep to the design probabilities in `.prob` and warn that `.weight`
differs. `.prob` itself must stay as
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) set it:
it supplies the probabilities of the PPS mappings above, and a sample
whose `.prob` was changed is refused.

Strata, intervals, cluster labels and population counts come from the
copy of the frame stored with the sample, not from the sample's own
columns, and are passed in the columns `.stratum`, `.cluster`, `.psu`,
`.fpc`, `.fpc2` and `.unit`. Renaming or recoding a design column in the
sample after the draw does not change the survey design. Those six names
are reserved: a sample with a column of one of them is refused rather
than overwritten. Rename the column in the sample, which leaves the
design as drawn.

## Replicate weights and calibration

[`survey::as.svrepdesign()`](https://rdrr.io/pkg/survey/man/as.svrepdesign.html)
turns the design returned here into one with replicate weights. For
samples from
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md)
without `within`,
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)
without `balanced`,
[`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md),
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md)
and
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
without `size`, its jackknife gives the standard errors of
`ht_total(variance = "jackknife")`. It stops with an error on the
designs handed over with Brewer's approximation or as Poisson sampling,
those of
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md),
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
and
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
with `size`, and
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md)
over them, because `survey` builds no replicates for its own
probability-proportional-to-size and Poisson designs. That is a
limitation of `survey`. For a systematic sample, alone or within strata,
it builds a jackknife that ignores the order of the walk. For
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)
with `balanced = TRUE` and
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
without `size` it keeps the finite population correction of the first
stage, which `ht_total(variance = "jackknife")` leaves out, and warns
that it dropped the corrections after the first stage. Its variance is
then the jackknife's times `1 - n / N`, for `n` of the `N` clusters
selected.

[`replicate_weights()`](https://elkronos.github.io/dRawn/reference/replicate_weights.md)
builds replicate weights from the sample itself, for
[`survey::svrepdesign()`](https://rdrr.io/pkg/survey/man/svrepdesign.html):
the jackknife of `ht_total(variance = "jackknife")`, for balanced
clusters and the designs handed over with Brewer's approximation too,
and successive differences along the walk for systematic samples.
Poisson sampling, which the jackknife understates, has none, and
`survey` estimates it from the design returned here.

[`frame_totals()`](https://elkronos.github.io/dRawn/reference/frame_totals.md)
gives the population totals that
[`survey::calibrate()`](https://rdrr.io/pkg/survey/man/calibrate.html)
and
[`survey::postStratify()`](https://rdrr.io/pkg/survey/man/postStratify.html)
take, from the frame stored with the sample, over the rows the design
could reach. Calibration is left to `survey`, on the design returned
here.

## See also

[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md),
[`replicate_weights()`](https://elkronos.github.io/dRawn/reference/replicate_weights.md),
[`frame_totals()`](https://elkronos.github.io/dRawn/reference/frame_totals.md)

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
coef(summary(survey::svyglm(spend ~ site, des)))
#>              Estimate Std. Error    t value     Pr(>|t|)
#> (Intercept) 286.00000   21.97916 13.0123281 3.915560e-18
#> siteb       -38.86667   42.03352 -0.9246588 3.593344e-01
#> sitec        11.66667   54.65988  0.2134411 8.318019e-01
#> sited       -72.00000   57.41243 -1.2540839 2.153146e-01
```
