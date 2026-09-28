# Stratified sampling

Draws a sample from each stratum. Allocation uses the largest-remainder
method, so the returned row count matches `n` exactly rather than
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

  How `n` is split across strata. `"proportional"` gives each stratum a
  share of `n` in proportion to its size; `"equal"` splits `n` evenly;
  `"neyman"` gives shares proportional to `size * sd`, using the column
  named by `allocation_by` (Neyman 1934). Neyman minimises the variance
  of a total for a fixed `n` by putting more rows where the values vary
  most, and is the right choice when you have a frame variable
  correlated with what you are measuring. Under any rule, a stratum that
  would be allocated more rows than it holds is taken whole and the rest
  of `n` re-split over the others by the same rule (Cochran 1977,
  section 5.9).

- allocation_by:

  Column whose within-stratum standard deviation drives
  `allocation = "neyman"`. Ignored otherwise.

- min_per_stratum:

  Minimum rows from each stratum. With the default of `0`, a stratum
  smaller than about `N / n` can be allocated no rows at all; its rows
  then have inclusion probability 0 and a total estimated from the
  sample silently leaves them out, so
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) warns
  when that happens. `1` makes every stratum reachable, which is what an
  unbiased Horvitz-Thompson total needs, and `2` also lets every stratum
  contribute to the variance estimate. Over-sampling a small stratum
  this way does *not* bias an estimate made with the design weights,
  which correct for it; it only moves precision around.

- replace:

  Sample with replacement within each stratum?

- na_rm:

  Drop rows whose stratum key is `NA` instead of raising an error.

## Value

A design object, for use with
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

## References

Neyman, J. (1934). On the two different aspects of the representative
method. *Journal of the Royal Statistical Society*, 97, 558–625.

Cochran, W. G. (1977). *Sampling Techniques*, 3rd ed. Wiley.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md),
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md),
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md),
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

# Rare strata are covered only if you ask -- without the floor, the one
# "rare" row gets no allocation, and draw() warns that it is unreachable
skewed <- data.frame(id = 1:1000, g = c(rep("common", 999), "rare"))
draw(skewed, design_stratified("g", n = 10, min_per_stratum = 1), seed = 1)$g
#>  [1] "common" "common" "common" "common" "common" "common" "common" "common"
#>  [9] "common" "rare"  
```
