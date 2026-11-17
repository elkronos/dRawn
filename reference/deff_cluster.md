# Design effect of a cluster sample, from the intraclass correlation

The design effect of taking `m` rows from each cluster, for
[`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md)
before a pilot has measured one with
[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md):
`1 + (m - 1) * icc`, with `icc` the intraclass correlation of the
measure, the correlation between two rows of one cluster (Kish 1965).
Rows of one cluster tend to be alike, so each row taken from a cluster
already sampled adds less than a row from a new cluster would, and the
more alike they are and the more taken from each, the larger the design
effect.

## Usage

``` r
deff_cluster(m, icc)
```

## Arguments

- m:

  Rows taken from each cluster, 1 or more. With one row a cluster, the
  design effect is 1 whatever `icc` is.

- icc:

  Intraclass correlation, from a pilot or a published study of the same
  measure, between `-1 / (m - 1)` and 1. A negative one, rows of a
  cluster less alike than rows taken at random, gives a design effect
  below 1.

## Value

The design effect, a number. `m` and `icc` may be vectors, recycled to a
common length, and give one design effect per pair.

## Details

The formula assumes every cluster gives the same number of rows, drawn
by simple random sampling of clusters.
[`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md)
counts rows, so a plan of `n` rows in clusters of `m` takes
`ceiling(n / m)` clusters, and the "Where the formula comes from"
section of
[`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md)
explains why a few clusters call for more than the plan.

## References

Kish, L. (1965). *Survey Sampling*. Wiley.

## See also

[`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md),
[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md) to
measure the design effect of a past sample.

## Examples

``` r
deff_cluster(m = 10, icc = 0.05)
#> [1] 1.45

# Rows to draw, and the clusters of ten rows that hold them
p <- plan_size(margin = 5, sd = 40, N = 20000,
               deff = deff_cluster(m = 10, icc = 0.05))
p$n
#> [1] 351
ceiling(p$n / 10)
#> [1] 36

# Takes of 5, 10 and 20 rows a cluster, side by side
plan_size(margin = 5, sd = 40, N = 20000,
          deff = deff_cluster(m = c(5, 10, 20), icc = 0.05))
#>     n n_effective margin level     N deff response target spread  p capped
#> 1 291         291      5  0.95 20000 1.20        1   mean     40 NA  FALSE
#> 2 351         351      5  0.95 20000 1.45        1   mean     40 NA  FALSE
#> 3 469         469      5  0.95 20000 1.95        1   mean     40 NA  FALSE
#>   short raised n_needed
#> 1 FALSE  FALSE      291
#> 2 FALSE  FALSE      351
#> 3 FALSE  FALSE      469
```
