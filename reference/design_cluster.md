# Cluster sampling

Selects whole clusters at random and returns all of their rows. The row
count follows from which clusters were picked, so there is no \`n\`; use
\[design_multistage()\] when you need to control it.

## Usage

``` r
design_cluster(clusters, n_clusters, balanced = FALSE, na_rm = FALSE)
```

## Arguments

- clusters:

  Column naming each row's cluster.

- n_clusters:

  Number of clusters to select.

- balanced:

  Take an equal number of rows from each selected cluster, equal to the
  smallest selected cluster's size.

- na_rm:

  Drop rows whose cluster label is \`NA\` instead of raising an error.
  When \`FALSE\`, missing labels are never treated as a cluster of their
  own.

## Value

A design object, for use with \[draw()\].

## See also

\[draw()\]

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
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
unique(draw(df, design_cluster("site", n_clusters = 3), seed = 1)$site)
#> [1] "s4" "s7" "s9"
```
