# Systematic sampling

Takes every `interval`-th row, beginning at `start`. The number of rows
drawn follows from the interval and the size of the data, so there is no
`n`.

## Usage

``` r
design_systematic(interval, start = NULL, order_by = NULL, na_rm = FALSE)
```

## Arguments

- interval:

  Sampling interval. A whole number of 1 or more. Fractional intervals
  are rejected: they produce uneven gaps, which is not a systematic
  design.

- start:

  Starting row, in `1:interval`. Drawn at random from that range when
  `NULL`.

- order_by:

  Optional column to sort by before walking the data. It changes which
  rows can appear together, so
  [`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
  and
  [`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)
  compute against the sorted order too.

- na_rm:

  When `order_by` is given, drop rows whose sort key is `NA` instead of
  raising an error.

## Value

A design object, for use with
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

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
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
[`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md),
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)

## Examples

``` r
df <- data.frame(id = 1:100, value = (1:100) / 10)
draw(df, design_systematic(interval = 10, start = 3))
#>    id value
#> 1   3   0.3
#> 2  13   1.3
#> 3  23   2.3
#> 4  33   3.3
#> 5  43   4.3
#> 6  53   5.3
#> 7  63   6.3
#> 8  73   7.3
#> 9  83   8.3
#> 10 93   9.3
draw(df, design_systematic(interval = 25, order_by = "value"), seed = 1)
#>    id value
#> 1  25   2.5
#> 2  50   5.0
#> 3  75   7.5
#> 4 100  10.0
```
