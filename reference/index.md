# Package index

## Drawing a sample

One verb, and the probabilities that make a sample usable.

- [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) : Draw
  a sample
- [`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
  [`sampling_weight()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
  : Inclusion probabilities and design weights
- [`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)
  : Joint inclusion probabilities
- [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  : Estimate a population total from a sample
- [`designs`](https://elkronos.github.io/dRawn/reference/designs.md) :
  Sampling designs
- [`is_design()`](https://elkronos.github.io/dRawn/reference/is_design.md)
  : Is this a sampling design?

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
