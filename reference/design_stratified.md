# Stratified sampling

Draws a sample from each stratum. Allocation uses the largest-remainder
method, so the returned row count matches \`n\` exactly rather than
drifting with per-stratum rounding.

## Usage

``` r
design_stratified(
  strata,
  n,
  allocation = c("proportional", "equal", "neyman"),
  allocation_by = NULL,
  min_per_stratum = 0L,
  replace = FALSE,
  na_rm = FALSE
)
```

## Arguments

- strata:

  One or more column names defining the strata. Several columns are
  cross-classified.

- n:

  Total rows to draw across all strata, under either allocation.

- allocation:

  How \`n\` is split across strata. \`"proportional"\` gives each
  stratum a share of \`n\` in proportion to its size; \`"equal"\` splits
  \`n\` evenly; \`"neyman"\` gives shares proportional to \`size \*
  sd\`, using the column named by \`allocation_by\`. Neyman minimises
  the variance of a total for a fixed \`n\` by putting more rows where
  the values vary most, and is the right choice when you have a frame
  variable correlated with what you are measuring.

- allocation_by:

  Column whose within-stratum standard deviation drives \`allocation =
  "neyman"\`. Ignored otherwise.

- min_per_stratum:

  Minimum rows from each stratum. The default of \`0\` leaves allocation
  unbiased; \`1\` guarantees coverage of rare strata at the cost of
  over-representing them.

- replace:

  Sample with replacement within each stratum?

- na_rm:

  Drop rows whose stratum key is \`NA\` instead of raising an error.

## Value

A design object, for use with \[draw()\].

## See also

\[draw()\]

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md),
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md),
[`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md),
[`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md),
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)

## Examples

``` r
df <- data.frame(id = 1:100, site = rep(letters[1:4], each = 25))
table(draw(df, design_stratified("site", n = 20), seed = 1)$site)
#> 
#> a b c d 
#> 5 5 5 5 

# Rare strata are covered only if you ask
skewed <- data.frame(id = 1:1000, g = c(rep("common", 999), "rare"))
draw(skewed, design_stratified("g", n = 10, min_per_stratum = 1), seed = 1)$g
#>  [1] "common" "common" "common" "common" "common" "common" "common" "common"
#>  [9] "common" "rare"  
```
