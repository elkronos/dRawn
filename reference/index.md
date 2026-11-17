# Package index

## Drawing a sample

One verb, and the designs it takes.

- [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) : Draw
  a sample
- [`designs`](https://elkronos.github.io/dRawn/reference/designs.md) :
  Sampling designs
- [`plot(`*`<drawn_design>`*`)`](https://elkronos.github.io/dRawn/reference/plot.drawn_design.md)
  : Picture a sampling design
- [`is_design()`](https://elkronos.github.io/dRawn/reference/is_design.md)
  : Is this a sampling design?

## What a design knows before you draw

The probabilities that make a sample estimable.

- [`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
  [`sampling_weight()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
  : Inclusion probabilities and design weights
- [`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)
  : Joint inclusion probabilities
- [`certainty_set()`](https://elkronos.github.io/dRawn/reference/certainty_set.md)
  : Rows a size-proportional design takes with certainty
- [`allocation_table()`](https://elkronos.github.io/dRawn/reference/allocation_table.md)
  : How a design splits its sample over the strata of a frame

## Estimating from a sample

Totals, means, ratios, proportions and quantiles, overall or by domain,
with the variance estimator the design calls for and a design effect.

- [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  : Estimate a population total from a sample
- [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md) :
  Estimate a population mean
- [`ht_ratio()`](https://elkronos.github.io/dRawn/reference/ht_ratio.md)
  : Estimate a ratio of two population totals
- [`ht_quantile()`](https://elkronos.github.io/dRawn/reference/ht_quantile.md)
  : Estimate population quantiles from a sample
- [`deff()`](https://elkronos.github.io/dRawn/reference/deff.md) :
  Design effect

## From fieldwork to estimate

Join fieldwork results to a sample, replace units found unavailable from
a reserve list, rebuild a sample read back from a file, and drop frame
columns before saving or sharing one.

- [`add_fieldwork()`](https://elkronos.github.io/dRawn/reference/add_fieldwork.md)
  : Join fieldwork results to a weighted sample
- [`use_reserves()`](https://elkronos.github.io/dRawn/reference/use_reserves.md)
  : Replace unavailable units from a reserve list
- [`as_drawn()`](https://elkronos.github.io/dRawn/reference/as_drawn.md)
  : Rebuild a weighted sample from its rows
- [`compact_sample()`](https://elkronos.github.io/dRawn/reference/compact_sample.md)
  : Drop the frame columns a sample's design does not need

## Planning and checking

How many rows to draw, what you actually got, and where to go next.

- [`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md)
  : How large a sample do you need?
- [`deff_cluster()`](https://elkronos.github.io/dRawn/reference/deff_cluster.md)
  : Design effect of a cluster sample, from the intraclass correlation
- [`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
  : Describe a drawn sample
- [`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
  : Hand a sample to the survey package
- [`replicate_weights()`](https://elkronos.github.io/dRawn/reference/replicate_weights.md)
  : Replicate weights for a sample
- [`frame_totals()`](https://elkronos.github.io/dRawn/reference/frame_totals.md)
  : Population totals from the frame behind a sample

## Equal-probability designs

Every unit equally likely, or equally likely within a group.

- [`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md)
  : Simple random sampling
- [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md)
  : Systematic sampling
- [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md)
  : Stratified sampling
- [`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md)
  : Reservoir sampling

## Cluster designs

Sample groups, then optionally sample within them.

- [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)
  : Cluster sampling
- [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
  : Multi-stage sampling

## Unequal-probability designs

Selection driven by size, and units too big to leave to chance.

- [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
  : Weighted sampling
- [`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md)
  : Take some rows with certainty, sample the rest

## Space and time

A sample spread evenly over a map or over covariates, restricted to a
region, or taken across time intervals.

- [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
  : Spatially balanced sampling
- [`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md)
  : Sampling within a region
- [`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md)
  : Temporal sampling

## Resampling

- [`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md)
  : Bootstrap sampling

## Extending

A design class of your own, drawn by its own method.

- [`draw_design()`](https://elkronos.github.io/dRawn/reference/draw_design.md)
  : Draw method for a design class
