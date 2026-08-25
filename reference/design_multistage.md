# Multi-stage sampling

Selects clusters at the first stage, then draws rows within them at the
second. \`n\` is the total across the selected clusters under both
allocations, matching \[design_stratified()\].

## Usage

``` r
design_multistage(
  clusters,
  n_clusters,
  n,
  allocation = c("equal", "proportional"),
  min_per_cluster = 0L,
  replace = FALSE,
  na_rm = FALSE
)
```

## Arguments

- clusters:

  Column naming each row's cluster.

- n_clusters:

  Number of clusters to select at stage one.

- n:

  Total rows to draw across the selected clusters.

- allocation:

  \`"equal"\` splits \`n\` evenly across the selected clusters;
  \`"proportional"\` splits it in proportion to their size.

- min_per_cluster:

  Minimum rows from each selected cluster. Defaults to \`0\`, which
  leaves allocation unbiased.

- replace:

  Sample with replacement within each cluster?

- na_rm:

  Drop rows whose cluster label is \`NA\` instead of raising an error.

## Value

A design object, for use with \[draw()\].

## See also

\[draw()\]

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md),
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md),
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
[`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md),
[`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md),
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)

## Examples

``` r
df <- data.frame(id = 1:100, site = rep(paste0("s", 1:10), each = 10))
nrow(draw(df, design_multistage("site", n_clusters = 4, n = 12), seed = 1))
#> [1] 12
```
