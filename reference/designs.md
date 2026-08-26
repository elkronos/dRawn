# Sampling designs

A design describes *how* to sample, independently of *what* you are
sampling. Build one with a `design_*()` constructor, then hand it to
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with the
data. The same design can be reused across data sets, stored, printed
and passed around.

## The designs

- [`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md):

  Rows uniformly at random.

- [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md):

  A share of each stratum.

- [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md):

  Every *k*-th row.

- [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md):

  Whole clusters.

- [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md):

  Clusters, then rows within them.

- [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md):

  Rows with probability governed by a weight.

- [`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md):

  Everything above a threshold, plus a sample of the rest.

- [`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md):

  A fixed-size sample from a stream.

- [`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md):

  Resampled replicates.

- [`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md):

  A share of each time interval.

- [`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md):

  Rows inside a region.

## What every design guarantees

Arguments mean the same thing everywhere they appear:

- `n` is always the **total** number of rows drawn, never a per-group
  figure.
  [`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md)
  uses `per_interval` precisely because that one *is* per group, and
  [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)
  has no `n` at all because the row count follows from which clusters
  were selected.

- `allocation` always says how a total is split across groups, either
  `"proportional"` or `"equal"`. Splitting is exact: the
  largest-remainder method is used, so a request for 60 rows returns 60
  rather than drifting with per-group rounding.

- `na_rm` always decides whether rows with a missing key are dropped or
  raise an error. It is never silently assumed either way.

- `replace` always means sampling with replacement within whatever group
  the design works on, and always rules out `draw(weights = TRUE)`: an
  inclusion probability describes distinct units, and a sample holding
  duplicates cannot be weighted by one.

- [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) always
  restores the caller's random number stream before returning, and
  always gives back a data frame with the input's class and column
  order.

- Rows come back in frame order for every design that selects a set of
  rows. The two exceptions are the ones where draw order is meaningful:
  [`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md)
  and
  [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
  return rows in the order they were drawn, and
  [`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md)
  returns replicates in order with a leading `.replicate` column.

## Estimating from a sample

Because a design is a value, it can report its own inclusion
probabilities before any sampling happens.
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
gives first-order probabilities and
[`sampling_weight()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
their reciprocals — the number of population rows each sampled row
stands for.
[`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)
gives second-order probabilities, and
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
and [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
combine them into a population total or mean with a standard error.
[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md) reports
what the design cost in precision against simple random sampling, and
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
reports what was actually drawn against what was in the frame.

Going the other way,
[`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md)
solves for the sample size a given margin of error requires — the step
*before* choosing a design. For analysis this package does not do,
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
hands a sample to the `survey` package.

Not every design has a closed form for these — see
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
for which, why, and what to do instead. Those that do not say so rather
than returning an approximation.

## Choosing one

- Every unit equally likely:

  [`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
  or
  [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md)
  when the frame has a useful order, or
  [`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md)
  when it does not fit in memory.

- Guaranteed coverage of subgroups:

  [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
  with `min_per_stratum` if rare groups must appear.

- Fieldwork cost matters more than efficiency:

  [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)
  or
  [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
  — visiting five sites is cheaper than visiting fifty, at the cost of
  precision.

- Large units matter more:

  [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md),
  with `method = "systematic"` or `"poisson"` if you intend to estimate.

- A few units dominate the total:

  [`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md)
  — take those with certainty and sample the tail, which removes them
  from the variance entirely.

- Coverage across time or space:

  [`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md),
  [`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md).

- Uncertainty of a statistic, not a population total:

  [`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md).

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md),
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md),
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md),
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
