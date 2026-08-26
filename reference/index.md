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
- [`draw_design()`](https://elkronos.github.io/dRawn/reference/draw_design.md)
  : Design-specific draw method

## What a design knows before you draw

The probabilities that make a sample estimable.

- [`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
  [`sampling_weight()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
  : Inclusion probabilities and design weights
- [`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)
  : Joint inclusion probabilities

## Estimating from a sample

A total or a mean, with a standard error and the design’s price.

- [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  : Estimate a population total from a sample
- [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md) :
  Estimate a population mean
- [`deff()`](https://elkronos.github.io/dRawn/reference/deff.md) :
  Design effect

## Planning and checking

How many rows to draw, what you actually got, and where to go next.

- [`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md)
  : How large a sample do you need?
- [`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
  : Describe a drawn sample
- [`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
  : Hand a sample to the survey package

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

Selection driven by a weight, a time window, or a place.

- [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
  : Weighted sampling
- [`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md)
  : Temporal sampling
- [`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md)
  : Spatial sampling

## Resampling

- [`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md)
  : Bootstrap sampling

## Composite designs

- [`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md)
  : Take some rows with certainty, sample the rest
