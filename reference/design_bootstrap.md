# Bootstrap sampling

Generates bootstrap replicates. \`"simple"\` resamples rows
independently; \`"block"\` is a moving-block bootstrap, which preserves
the serial dependence in ordered data by concatenating \`ceiling(n /
block_length)\` independently chosen blocks per replicate.

## Usage

``` r
design_bootstrap(
  n_replicates = 1000L,
  n = NULL,
  method = c("simple", "block"),
  block_length = NULL
)
```

## Arguments

- n_replicates:

  Number of replicates to generate.

- n:

  Rows per replicate. \`NULL\` uses \`nrow(data)\`, the standard
  nonparametric bootstrap.

- method:

  \`"simple"\` or \`"block"\`.

- block_length:

  Block length for \`method = "block"\`. \`NULL\` uses
  \`floor(nrow(data)^(1/3))\`, a common rule of thumb.

## Value

A design object, for use with \[draw()\].

## Details

\[draw()\] returns all replicates in one data frame with a leading
\`.replicate\` column. Split them with \`split(out, out\$.replicate)\`.

## See also

\[draw()\]

Other designs:
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
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
df <- data.frame(id = 1:100, value = (1:100) / 10)

reps <- draw(df, design_bootstrap(n_replicates = 5, n = 20), seed = 1)
table(reps$.replicate)
#> 
#>  1  2  3  4  5 
#> 20 20 20 20 20 

# A statistic per replicate
vapply(split(reps, reps$.replicate), function(r) mean(r$value), numeric(1))
#>     1     2     3     4     5 
#> 5.450 5.285 4.395 4.735 5.570 
```
