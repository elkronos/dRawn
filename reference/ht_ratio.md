# Estimate a ratio of two population totals

`ht_ratio()` estimates `R = Y / X`, the ratio of the population totals
of `y` and `x`, by `sum(y / pi) / sum(x / pi)`: audited value per unit
of book value, spend per visit, or yield per hectare. A Hajek mean is
the ratio with `x` equal to 1 in every row.
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
with `ratio` multiplies the ratio by the frame's total of `x` to
estimate the total of `y`.

## Usage

``` r
ht_ratio(
  sample,
  y,
  x,
  by = NULL,
  variance = c("auto", "analytic", "jackknife", "none"),
  level = 0.95,
  df = NULL,
  alternative = c("two.sided", "less", "greater")
)
```

## Arguments

- sample:

  A data frame returned by
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with
  `weights = TRUE`. See "Which columns are read".

- y, x:

  The numerator and the denominator: each a column name of `sample`, or
  a numeric vector as long as `sample`.

- by:

  Optional column name(s) in `sample` defining domains. When given, the
  result is one row per domain. See "Domains". A column named like one
  of the table's own, such as `n` or `se`, is refused.

- variance:

  How to compute it. `"auto"` uses the design's own estimator (see
  "Variance") and falls back to the jackknife when that declines;
  `"analytic"` insists on the design's own estimator, returning `NA`
  with the reason in `note` rather than falling back; `"jackknife"` uses
  the jackknife, and returns `NA` with the reason in `note` for
  systematic and Poisson samples, where it does not apply; `"none"`
  skips it. The result reports which was used in `method`, and `"none"`
  when nothing could produce a figure.

- level:

  Confidence level for the interval.

- df:

  Degrees of freedom for the interval. `NULL` uses the design's (see
  "Confidence intervals"); `Inf` gives a normal interval; otherwise a
  single number of 1 or more.

- alternative:

  `"two.sided"` for a confidence interval, `"less"` for a one-sided
  upper bound or `"greater"` for a lower bound, each at `level`. See
  "Confidence intervals".

## Value

A list of class `drawn_ratio` with a
[`print()`](https://rdrr.io/r/base/print.html) method, holding `ratio`,
the estimate; `variable` and `x`, the names of `y` and `x` when they are
column names, else `NULL`; `total_y` and `total_x`, their estimated
totals; and `variance`, `se`, `ci`, `level`, `alternative`, `df`, `n`,
`n_random`, `ci_scale`, `design`, `deff`, `method`, `note` and
`unreachable` as for
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md).
With `by`, a `drawn_by` table with a row per domain and a `ratio`
column, as for
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md).
Both have [`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
`tidy()` methods; see "Working with the result" in
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md).

## Variance

The ratio is linearised. Its variance is that of the Horvitz-Thompson
total of the residuals `(y - R x) / X_hat`, with `X_hat = sum(x / pi)`
(Särndal, Swensson and Wretman 1992; Deville 1999), and is computed with
the design's own estimator, as for a total: see "Variance" and
"Confidence intervals" in
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
for the estimators, the fallback to the jackknife, the degrees of
freedom and the notes. For the designs
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
hands to `survey` with identical standard errors, the standard error
equals that of
[`survey::svyratio()`](https://rdrr.io/pkg/survey/man/svyratio.html) on
its result. For systematic, spread and systematic PPS samples it comes
from the estimator this package uses for their totals, which `survey`
does not have.

The residuals sum to zero over the sample, weighted by `1 / pi`, so a
ratio whose sampled rows lie in one primary sampling unit, one row or
one cluster, gets no standard error, as a Hajek mean does. When `y` is
proportional to `x` in every sampled row the residuals are 0, the
variance is set to exactly 0 rather than left to rounding error, and
`note` says so.

The ratio of two estimated totals is biased, by an amount of order
`1 / n`, which is small against its standard error except in small
samples. The design effect compares its variance with that of the same
ratio from a simple random sample of the same size, as
`survey::svyratio(deff = TRUE)` does.

## Domains

`by` gives a ratio for each domain, `R_d = Y_d / X_d` over the rows of
the domain, estimated from the whole sample. Each domain's residuals,
`(y - R_d x) / X_hat_d`, are 0 outside it, as for the domain totals and
means described under "Domains" in
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
and the standard errors are those of
`survey::svyby(~y, ~d, design, survey::svyratio, denominator = ~x)`. A
domain whose estimated total of `x` is 0 is refused by name.

A missing `y` or `x` is refused. `survey::svyratio(na.rm = TRUE)` on
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
estimates the ratio over the rows that have both.

## References

Deville, J.-C. (1999). Variance estimation for complex statistics and
estimators: linearization and residual techniques. *Survey Methodology*,
25, 193–203.

Särndal, C.-E., Swensson, B. and Wretman, J. (1992). *Model Assisted
Survey Sampling*. Springer.
[doi:10.1007/978-1-4612-4378-6](https://doi.org/10.1007/978-1-4612-4378-6)

## See also

[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md),
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)

## Examples

``` r
set.seed(1)
pop <- data.frame(
  id = 1:300,
  site = rep(c("a", "b", "c"), each = 100),
  book = round(stats::rlnorm(300, 6, 1))
)
pop$audited <- round(pop$book * stats::runif(300, 0.9, 1.05))
s <- draw(pop, design_stratified("site", n = 45), seed = 1, weights = TRUE)

# Audited value per unit of book value
ht_ratio(s, "audited", "book")
#> Ratio of Horvitz-Thompson totals  (stratified design, n = 45)
#>   ratio    `audited` / `book`
#>   totals   173,026.7 / 177,573.3
#>   estimate 0.9743956
#>   se       0.007899163  (analytic)
#>   95% CI  0.9584544 to 0.9903367  (t, 42 df)
#>   deff     1.01  (about the same as simple random sampling)
sum(pop$audited) / sum(pop$book)   # the truth
#> [1] 0.9718688

# For each site, from the whole sample
ht_ratio(s, "audited", "book", by = "site")
#> Ratio `audited` / `book` by domain  (stratified design, 95% CI, t, 42 df)
#>  site  n  ratio      se ci_lower ci_upper   method
#>     a 15 0.9874 0.01226   0.9627   1.0122 analytic
#>     b 15 0.9707 0.01336   0.9438   0.9977 analytic
#>     c 15 0.9629 0.01515   0.9323   0.9935 analytic

# The ratio estimator of the total uses the frame's total of `book`, and
# is far more precise than the Horvitz-Thompson total here
ht_total(s, "audited", ratio = "book")
#> Ratio estimate of the total  (stratified design, n = 45)
#>   ratio    0.9743956, `audited` / `book` in the sample
#>   x total  198,996, `book` over the frame
#>   estimate 193,900.8
#>   se       1,571.902  (analytic)
#>   95% CI  190,728.6 to 197,073  (t, 42 df)
#>   deff     1.01  (about the same as simple random sampling)
ht_total(s, "audited")
#> Horvitz-Thompson total  (stratified design, n = 45)
#>   estimate 173,026.7
#>   se       20,306.49  (analytic)
#>   95% CI  132,046.5 to 214,006.8  (t, 42 df)
#>   deff     1.03  (about the same as simple random sampling)
sum(pop$audited)
#> [1] 193398
```
