# Estimate a population mean

The Hajek estimator, `sum(y / pi) / sum(1 / pi)`, is the default. It
divides by the *estimated* population size rather than the known one, so
a sample that happens to over-represent heavy-weight rows inflates
numerator and denominator together and they partly cancel.
`estimator = "ht"` divides by the true `N` instead.

## Usage

``` r
ht_mean(
  sample,
  y,
  estimator = c("hajek", "ht"),
  variance = c("auto", "analytic", "jackknife", "none"),
  level = 0.95
)
```

## Arguments

- sample:

  A data frame returned by
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with
  `weights = TRUE`.

- y:

  The variable to average: a column name, or a numeric vector as long as
  `sample`.

- estimator:

  `"hajek"` or `"ht"`. See above.

- variance:

  Passed to the underlying total. See
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md).

- level:

  Confidence level for the interval.

## Value

A list with a [`print()`](https://rdrr.io/r/base/print.html) method,
holding `mean`, `variance`, `se`, `ci`, `level`, `n`, `design`, `deff`,
`method` and `note`. See
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
for what `method` and `note` say, and
[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md) for the
design effect.

## When the two differ, and which to use

They coincide **exactly** whenever the design weights of the rows you
drew sum to `N` — which covers every fixed-size equal-probability design
and every
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md)
without replacement, because each stratum contributes
`n_h * N_h / n_h = N_h`. There is nothing to choose between them there.

They differ when that sum is random or uneven:

- the sample size itself is random —
  [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
  with `method = "poisson"`,
  [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)
  over clusters of unequal size, or
  [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md)
  where the interval does not divide `N`;

- the weights vary within a fixed-size sample —
  probability-proportional-to- size selection.

Hajek is usually the steadier of the two and is the default for that
reason. The exception is worth knowing: when `y` is close to
proportional to the size measure that drove selection, `y / pi` is
nearly constant, the Horvitz-Thompson numerator barely moves, and
dividing it by the known `N` beats dividing by an estimate.

One caveat on `"ht"`. It is unbiased for the frame mean *provided every
row could have been selected*. Rows with inclusion probability 0 —
outside a time window, outside a region, zero weight — sit inside the
`N` it divides by but can never enter the numerator, so the estimate is
biased low by exactly their share of the frame.
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
reports how many such rows there are. The Hajek mean is unaffected,
because it estimates the mean of the part of the frame the design can
actually reach.

## See also

[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md),
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)

## Examples

``` r
set.seed(1)
pop <- data.frame(
  id = 1:200,
  site = rep(c("a", "b"), times = c(150, 50)),
  spend = round(stats::runif(200, 10, 500))
)
s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)

ht_mean(s, "spend")
#> Hajek mean  (stratified design, n = 40)
#>   estimate 270
#>   se       18.56849  (analytic)
#>   95% CI  233.6064 to 306.3936
#>   deff     1.03  (about the same as simple random sampling)
mean(pop$spend)   # the truth
#> [1] 263.66

# Stratified without replacement: the weights sum to N, so the two agree
ht_mean(s, "spend", estimator = "ht")$mean == ht_mean(s, "spend")$mean
#> [1] TRUE

# Poisson sampling has a random size, so they part company
p <- draw(pop, design_weighted("spend", n = 40, method = "poisson"),
          seed = 1, weights = TRUE)
c(hajek = ht_mean(p, "spend", variance = "none")$mean,
  ht    = ht_mean(p, "spend", "ht", variance = "none")$mean)
#> hajek    ht 
#>   NaN     0 
```
