# Design effect

How much precision the design costs against simple random sampling of
the same size. `deff = 1` means the design is doing as well as a coin
flip over the frame; above 1 it is doing worse, which is the usual price
of clustering; below 1 it is doing better, which is what stratification
and probability-proportional-to-size buy you.

## Usage

``` r
deff(x)
```

## Arguments

- x:

  A result from
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  or
  [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md).

## Value

A single number, or `NA` when the design has no variance estimate.

## Details

Read it as an exchange rate on sample size: at `deff = 2`, a sample of
400 carries about as much information as 200 drawn at random.

## See also

[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)

## Examples

``` r
set.seed(1)
# Sites differ from each other, and rows within a cluster are alike --
# exactly the structure that makes stratifying pay and clustering cost.
pop <- data.frame(
  id = 1:400,
  site = rep(c("a", "b", "c", "d"), each = 100),
  cl = rep(paste0("c", 1:40), each = 10)
)
pop$y <- rep(c(20, 60, 120, 200), each = 100) + round(stats::rnorm(400, 0, 8))

# Stratifying on something that matters buys precision (deff below 1)
deff(ht_total(draw(pop, design_stratified("site", n = 40), seed = 1,
                   weights = TRUE), "y"))
#> [1] 0.01174153

# Clustering usually costs it (deff above 1)
deff(ht_total(draw(pop, design_cluster("cl", n_clusters = 4), seed = 1,
                   weights = TRUE), "y"))
#> [1] 12.91977
```
