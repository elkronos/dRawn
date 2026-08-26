# Describe a drawn sample

What you actually got, against what the design asked for. Worth a look
before analysing: it surfaces strata that came up short, weights that
vary more than you expected, and rows the design could never have
reached.

## Usage

``` r
sample_summary(sample)
```

## Arguments

- sample:

  A data frame returned by
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with
  `weights = TRUE`.

## Value

A list with a [`print()`](https://rdrr.io/r/base/print.html) method,
holding:

- `design`:

  The design's type, as a string.

- `n`, `N`:

  Rows drawn, and rows in the frame.

- `weight_range`:

  The smallest and largest design weight.

- `weight_cv`:

  Their coefficient of variation. Large values mean a few rows carry
  most of the estimate.

- `unreachable`:

  Frame rows with inclusion probability 0 — the design could never have
  selected them. `NA` if the design has no closed-form inclusion
  probability.

- `by_group`:

  A data frame of `group`, `drawn`, `in_frame` and `rate` per stratum or
  cluster, or `NULL` for a design with no grouping.

- `group_col`:

  The column(s) `by_group` is keyed on.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md),
[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md)

## Examples

``` r
set.seed(1)
pop <- data.frame(
  id = 1:400,
  site = rep(c("a", "b", "c", "d"), times = c(200, 100, 60, 40))
)
s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
sample_summary(s)
#> Sample of 40 from 400  (stratified design)
#>   sampling fraction  0.1
#>   design weights     10 to 10   (cv 0)
#> 
#>   by site:
#>     group  drawn  in frame   rate
#>     a         20      200  0.100
#>     b         10      100  0.100
#>     c          6       60  0.100
#>     d          4       40  0.100
```
