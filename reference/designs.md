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

  A share of each stratum, at random or walked in order.

- [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md):

  Every *k*-th row, or `n` rows at an even step, from one random start
  or several.

- [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md):

  Whole clusters.

- [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md):

  Clusters, at random or in proportion to size, then rows within them.

- [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md):

  Rows with probability governed by a weight, within strata if asked.

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

- [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md):

  A sample spread evenly across a map or across auxiliary variables,
  within strata if asked.

## What every design guarantees

Arguments mean the same thing everywhere they appear:

- `n` is always the **total** number of rows drawn, never a per-group
  figure.
  [`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md)
  uses `per_interval` precisely because that one *is* per group, and
  [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)
  has no `n` at all because the row count follows from which clusters
  were selected. A plan from
  [`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md)
  can be given as `n`, and its `n` is used.

- `allocation` always says how a total is split across groups:
  `"proportional"`, `"equal"` or, in
  [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md)
  only, `"neyman"`.
  [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
  [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
  and
  [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
  also take the count for each stratum in place of a rule, and `n` is
  then their total. Splitting is exact: the largest-remainder method is
  used, so a request for 60 rows returns 60 rather than drifting with
  per-group rounding.

- `na_rm` always decides whether rows with a missing key are dropped or
  raise an error. It is never silently assumed either way.

- Groups always come in key order, whether strata, `by` domains or the
  rows of
  [`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md):
  a factor's levels, or the sorted values of any other column, with text
  sorted byte by byte as in the C locale. A seed therefore draws the
  same rows in every locale, and an `order_by` column of text is walked
  in the same order everywhere.

- `replace` always means sampling with replacement within whatever group
  the design works on, and always rules out `draw(weights = TRUE)`: an
  inclusion probability describes distinct units, and a sample holding
  duplicates cannot be weighted by one.

- With a `seed`,
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) always
  restores the caller's random number stream and generator before
  returning, and a seed draws the same sample whatever generator the
  session has selected. With `seed = NULL` it draws from, and advances,
  the caller's stream. See "Seeding" in
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

- [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) always
  gives back a data frame with the input's class and column order, and
  the frame's own row names on the rows drawn.
  [`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md)
  with `output = "indices"` or `"weights"` gives a matrix instead, and
  [`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md)
  a list when `data` is not a data frame, unless `batches = TRUE`.

- Rows come back in frame order for every design that selects a set of
  rows, including
  [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md)
  with `order_by` and
  [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
  with `method = "systematic"` or `"poisson"`. The exceptions are the
  ones where draw order is meaningful:
  [`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md)
  and `design_weighted(method = "successive")` return rows in the order
  they were drawn, and
  [`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md)
  returns replicates in order with a leading `.replicate` column. A
  design with `reserve` returns its main sample and then its reserves,
  in the order they are used, with a leading `.reserve` column.

## Estimating from a sample

Because a design is a value, it can report its own inclusion
probabilities before any sampling happens.
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
gives first-order probabilities and
[`sampling_weight()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
their reciprocals — the number of population rows each sampled row
stands for, and
[`certainty_set()`](https://elkronos.github.io/dRawn/reference/certainty_set.md)
the rows a design that selects in proportion to size takes with
certainty.
[`allocation_table()`](https://elkronos.github.io/dRawn/reference/allocation_table.md)
shows how a design with strata splits `n` over the strata of a frame,
and which strata it holds at the floor, takes whole or leaves empty.
[`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)
gives second-order probabilities, and
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
and [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
turn a sample into a population total, mean or proportion, overall or by
domain and for one variable or several, with a standard error from the
variance estimator suited to the design and a confidence interval on the
design's degrees of freedom. A factor `y` gives the count or share of
each level.
[`ht_ratio()`](https://elkronos.github.io/dRawn/reference/ht_ratio.md)
estimates a ratio of two totals the same way, `ht_total(ratio =)` the
ratio estimator of a total, and
[`ht_quantile()`](https://elkronos.github.io/dRawn/reference/ht_quantile.md)
medians and other quantiles with Woodruff intervals.
[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md) reports
what the design cost in precision against simple random sampling, and
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
reports what was actually drawn against what was in the frame.
[`add_fieldwork()`](https://elkronos.github.io/dRawn/reference/add_fieldwork.md)
joins fieldwork results to a sample by an id,
[`use_reserves()`](https://elkronos.github.io/dRawn/reference/use_reserves.md)
replaces the units of a list from
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md)
or
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
with `reserve` that were found unavailable,
[`as_drawn()`](https://elkronos.github.io/dRawn/reference/as_drawn.md)
rebuilds a sample read back from a file or a database, and
[`compact_sample()`](https://elkronos.github.io/dRawn/reference/compact_sample.md)
drops the frame columns its design does not read before it is saved or
shared.

Going the other way,
[`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md)
solves for the sample size a given margin of error requires — the step
*before* choosing a design. It plans a sample overall or split over
strata, and
[`deff_cluster()`](https://elkronos.github.io/dRawn/reference/deff_cluster.md)
gives the design effect of a cluster sample to plan with. For analysis
this package does not do,
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
hands a sample to the `survey` package,
[`replicate_weights()`](https://elkronos.github.io/dRawn/reference/replicate_weights.md)
builds replicate weights for it from the design, and
[`frame_totals()`](https://elkronos.github.io/dRawn/reference/frame_totals.md)
gives the frame's totals to calibrate it to.

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
  precision. With `size`,
  [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
  selects the sites in proportion to size and takes the same number of
  rows from each.

- Large units matter more:

  [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md),
  with `method = "systematic"` or `"poisson"` if you intend to estimate,
  and `strata` to sample in proportion to size within each stratum.

- A few units dominate the total:

  [`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md)
  — take those with certainty and sample the tail, which removes them
  from the variance entirely.

- Coverage across time or space:

  [`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md)
  for time;
  [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
  for an even spread over a map, which is usually far more precise than
  a simple random sample when what you measure varies smoothly over
  space;
  [`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md)
  to restrict sampling to a region.

- Balance on known covariates without choosing strata:

  [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
  across those covariates.

- Uncertainty of a statistic, not a population total:

  [`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md).

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md),
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md),
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md),
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
