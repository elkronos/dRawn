# Multi-stage sampling

Selects clusters at the first stage, then draws rows within them at the
second. `n` is the total across the selected clusters under both
allocations, matching
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md).

## Usage

``` r
design_multistage(
  clusters,
  n_clusters,
  n,
  allocation = c("equal", "proportional"),
  min_per_cluster = 0L,
  replace = FALSE,
  na_rm = FALSE,
  size = NULL
)
```

## Arguments

- clusters:

  Column naming each row's cluster.

- n_clusters:

  Number of clusters to select at stage one.

- n:

  Total rows to draw across the selected clusters, at least
  `n_clusters`, or a plan from
  [`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md).
  A multiple of `n_clusters`, no larger than `n_clusters` times the
  smallest cluster, gives every selected cluster the same take and every
  row a closed-form inclusion probability. Otherwise the take depends on
  which clusters were selected, and `draw(weights = TRUE)` weights the
  sample conditional on them. See "When the take depends on the clusters
  selected".

- allocation:

  `"equal"` splits `n` evenly across the selected clusters;
  `"proportional"` splits it in proportion to their size, and its
  weights are conditional on the clusters selected. Proportional
  allocation can give a small selected cluster no rows;
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) then
  warns with class `drawn_empty_cluster` and names it.

- min_per_cluster:

  Minimum rows from each selected cluster. Defaults to `0`, and
  `n_clusters * min_per_cluster` must not exceed `n`. A cluster smaller
  than the minimum is taken whole. `draw(weights = TRUE)` needs at least
  1 under proportional allocation. Take at least two rows per cluster if
  you intend to estimate a variance: with one, the variation within
  clusters cannot be measured and
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  falls back to the ultimate-cluster jackknife.

- replace:

  Sample with replacement within each cluster? A selected cluster then
  gives its full share of `n` however small it is, and a row's inclusion
  probability is the chance of being drawn at least once,
  `(n_clusters / N_clusters) * (1 - (1 - 1 / N_h)^m)` for a take of `m`
  from a cluster of `N_h` rows.

- na_rm:

  Drop rows whose cluster label, or `size`, is `NA` instead of raising
  an error.

- size:

  Optional column holding each cluster's size measure, the same on every
  row of a cluster: positive, finite and numeric. Clusters are then
  selected with probability proportional to it. See "Clusters in
  proportion to size".

## Value

A design object, for use with
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

## When the take depends on the clusters selected

Without `size`, the `n` rows are split over the clusters selected, so a
cluster's take is the same in every sample only under `"equal"`
allocation of an `n` that divides by `n_clusters`, with no cluster
smaller than `n / n_clusters`. Otherwise, under `"proportional"`
allocation, with a remainder that goes to the largest clusters selected,
or with a small cluster capped at its size and its shortfall dealt to
the others, a row's inclusion probability depends on which other
clusters were selected. It is the chance of drawing the row given the
clusters selected, averaged over every possible selection, which has no
closed form, so
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
refuses it.

`draw(weights = TRUE)` weights such a sample conditional on the clusters
it selected. With `a` of `A` clusters selected and `m_i` of the `N_i`
rows of cluster `i` taken, `.prob` is `(a / A) * m_i / N_i`, and the
estimate of a total is `A / a` times the sum of the selected clusters'
estimated totals. That is unbiased whichever clusters were selected, and
so over all of them.
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
and [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
estimate its variance with the two-stage formula, the variance between
the estimated cluster totals plus that within each cluster at its own
take, on `a - 1` degrees of freedom, and
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
hands it to `survey` in two stages, where the standard errors agree. The
variance falls back to the delete-a-cluster jackknife when a selected
cluster gives one row of several. Proportional allocation needs
`min_per_cluster` of at least 1 for weights, since a selected cluster
with no rows would count as a cluster total of 0.

## Clusters in proportion to size

With `size`, the clusters are selected with probability proportional to
a size measure, such as a count of households from a census, by
systematic PPS as in
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md):
cluster `i` gets `pi_i = n_clusters * M_i / sum(M)`, and a cluster large
enough that this passes 1 is taken with certainty, the rest of
`n_clusters` being spread over the others. Each selected cluster then
gives the same `m = n / n_clusters` rows by simple random sampling, so a
row of cluster `i`, of `N_i` rows, has inclusion probability
`pi_i * m / N_i`. When the measure is the cluster's row count and no
cluster is certain, every row gets `n / N`: the design is
self-weighting, with a fixed number of interviews per cluster.

The measure is a property of the cluster, so every row of a cluster must
carry the same value; a cluster with two values is refused. `n` must be
a multiple of `n_clusters`, and `allocation` must be `"equal"`.

[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
and [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
estimate the variance as Deville's approximation over the estimated
totals of the selected clusters below certainty, each scaled by its own
`1 - pi_i`, plus the variance within the selected clusters, each divided
by its `pi_i`. That is the two-stage estimator of `design_multistage()`
without `size` when every `pi_i` is equal. A cluster taken with
certainty is a stratum, and adds only the variance of the simple random
sample within it. The degrees of freedom count the clusters below
certainty less one, and in each certain cluster the sampled rows less
one. The variance is declined when one cluster below certainty is
selected. With one row taken per cluster, whose variation within
clusters cannot be measured,
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
falls back to the delete-a-cluster jackknife over the clusters below
certainty, with no finite population correction, and declines if a
cluster taken with certainty gives one row of several. In simulations of
2,000 draws each from frames of 78 or 80 clusters of 10 to 600 rows,
taking 2 or 5 rows from each of 12 or 20 clusters, with and without
clusters taken with certainty, the variance averaged 0.94 to 1.02 times
the true variance of the total, and 95% intervals covered 93% to 95% of
the time.

[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
hands the design to `survey` in two stages with Brewer's approximation
at the first. A sample with one row from a selected cluster of several
goes over in one stage instead, as the clusters of that jackknife, and
the standard errors agree. Joint probabilities have no closed form;
[`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)
estimates them with `simulate = TRUE`.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
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
df <- data.frame(id = 1:100, site = rep(paste0("s", 1:10), each = 10))
nrow(draw(df, design_multistage("site", n_clusters = 4, n = 12), seed = 1))
#> [1] 12

# Villages selected in proportion to their households, then 5 households
# from each: every household has the same chance, 20 / 260
hh <- data.frame(village = rep(paste0("v", 1:8),
                               times = c(10, 20, 30, 40, 20, 50, 60, 30)))
hh$households <- ave(seq_len(nrow(hh)), hh$village, FUN = length)
pps <- design_multistage("village", n_clusters = 4, n = 20,
                         size = "households")
unique(round(inclusion_prob(hh, pps), 4))
#> [1] 0.0769
```
