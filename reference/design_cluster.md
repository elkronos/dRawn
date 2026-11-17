# Cluster sampling

Selects whole clusters at random and returns all of their rows. The row
count follows from which clusters were picked, so there is no `n`; use
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
when you need to control it.

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
  smallest selected cluster's size. See "Balanced clusters".

- na_rm:

  Drop rows whose cluster label is `NA` instead of raising an error.
  When `FALSE`, missing labels are never treated as a cluster of their
  own.

## Value

A design object, for use with
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

## Balanced clusters

With `balanced = TRUE`, each of the `a` clusters selected out of `A`
gives `m` rows by simple random sampling, `m` being the size of the
smallest of them. A row of a cluster of `N_i` rows has inclusion
probability `(a / A) * E(m) / N_i`, where `E(m)` is the average take
over the selections that include its cluster.
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
computes it exactly: given that cluster `i` is selected, the take is at
least `s` rows when the other `a - 1` clusters selected all hold `s`
rows or more, which happens with probability
`choose(c_s - 1, a - 1) / choose(A - 1, a - 1)`, `c_s` being the number
of clusters of `s` rows or more, and `E(m)` adds these up for `s` from 1
to `N_i`. On clusters of 2, 4, 6 and 8 rows with two selected, that
gives 0.5, 0.417, 0.333 and 0.25.

`draw(weights = TRUE)` weights a balanced sample conditional on the
clusters it selected instead: `.prob` is `(a / A) * m / N_i` for the `m`
taken in this draw. The estimate of a total is then `A / a` times the
sum of the selected clusters' estimated totals, which is unbiased
whichever clusters were selected, and so over all of them.
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
and [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
use the two-stage variance of
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
with the variance between the estimated cluster totals and within each
cluster, and
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
hands the sample to `survey` in two stages. The smallest selected
cluster is taken whole and adds no variance within it. When it holds one
row, every cluster gives one, the variation within clusters cannot be
measured, and the variance falls back to the delete-a-cluster jackknife.
Joint probabilities have no closed form;
[`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)
estimates them with `simulate = TRUE`.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
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
unique(draw(df, design_cluster("site", n_clusters = 3), seed = 1)$site)
#> [1] "s4" "s7" "s9"

# Balanced: the take is the smallest selected cluster's size, so a row of
# a large cluster is less likely to be drawn
uneven <- data.frame(site = rep(c("a", "b", "c", "d"), c(2, 4, 6, 8)))
bal <- design_cluster("site", n_clusters = 2, balanced = TRUE)
round(tapply(inclusion_prob(uneven, bal), uneven$site, unique), 3)
#>     a     b     c     d 
#> 0.500 0.417 0.333 0.250 
```
