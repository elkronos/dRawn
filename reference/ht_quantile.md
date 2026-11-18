# Estimate population quantiles from a sample

`ht_quantile()` estimates quantiles of `y`, such as the median or the
quartiles, from the Hajek estimate of its distribution function,
`F(t) = sum(w * (y <= t)) / sum(w)` with `w = 1 / pi`, and gives each a
Woodruff confidence interval and a standard error from the design's own
variance estimator.

## Usage

``` r
ht_quantile(
  sample,
  y,
  probs = 0.5,
  by = NULL,
  variance = c("auto", "analytic", "jackknife", "none"),
  level = 0.95,
  df = NULL
)
```

## Arguments

- sample:

  A data frame returned by
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with
  `weights = TRUE`. See "Which columns are read".

- y:

  The variable: a numeric column name of `sample`, or a numeric vector
  as long as `sample`. A logical `y` is refused, since its quantiles are
  0 or 1;
  [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
  estimates the share that is `TRUE`.

- probs:

  The probabilities, each from 0 to 1, such as `0.5` for the median or
  `c(0.25, 0.5, 0.75)` for the quartiles.

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

## Value

For one probability without `by`, a list of class `drawn_quantile` with
a [`print()`](https://rdrr.io/r/base/print.html) method, holding
`quantile`, the estimate; `prob`; `variable`, the name of `y` when it is
a column name, else `NULL`; `cdf` and `cdf_se`, the estimated
distribution function at the estimate and its standard error; `steps`,
the estimated distribution function as the sorted sampled values `x` and
the share of the weight `F` at or below each, which
[`confint()`](https://rdrr.io/r/stats/confint.html) reads to form the
interval at another level; and `variance`, `se`, `ci`, `level`,
`alternative` (always `"two.sided"`), `df`, `n`, `n_random`, `design`,
`deff`, `method`, `note` and `unreachable` as for
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
with `ci_scale` `"woodruff"`.

For several probabilities, or with `by`, a `drawn_by` table with a row
for each domain and probability: the `by` columns, `prob`, `n`,
`quantile`, `variance`, `se`, `ci_lower`, `ci_upper` and `method`.

Both have [`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
`tidy()` methods, as described under "Working with the result" in
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md).
[`confint()`](https://rdrr.io/r/stats/confint.html) at another level
carries the interval for `F(q)` at that level back through the quantile
rule. [`vcov()`](https://rdrr.io/r/stats/vcov.html) holds the variances,
and `NA` for the covariances between quantiles, which are not estimated.

## The quantile rule

The estimate at `p` is the smallest sampled value of `y` whose estimated
distribution function is `p` or more, which inverts the step function
`F` without interpolating between sampled values. It is always a value
some sampled row has. This is the rule of
`survey::svyquantile(qrule = "math")`, survey's default, and with equal
weights it is `stats::quantile(type = 1)`, apart from the tolerance
below. `p = 0` gives the smallest sampled value and `p = 1` the largest.

A share of the weight equal to `p` can come out of the sums a rounding
error below `p`, so a share within `1e-9` below `p` counts as reaching
it. A `p` that is a rounding error above a share therefore gives the
value at that share. `seq(0, 1, 0.01)[71]`, which is 0.7 plus a rounding
error, gives the 21st of 30 equally weighted values, where
`stats::quantile(type = 1)` and
[`survey::svyquantile()`](https://rdrr.io/pkg/survey/man/svyquantile.html)
give the 22nd.

## Interval and standard error

The interval is Woodruff's (1952). At the estimate `q`, `F(q)` is a
Hajek mean of the indicator `y <= q`, and its variance is the variance
of the total of the linearised residuals `((y <= q) - F(q)) / N_hat`,
with `N_hat = sum(w)`, from the design's own estimator, as for
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md).
"Variance" and "Confidence intervals" in
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
describe the estimators, the fallback to the jackknife, the degrees of
freedom and the notes. The interval `F(q)` plus or minus `t` standard
errors, with `t` on the design's degrees of freedom, is carried back
through the quantile rule, so its ends are the estimated quantiles at
its two ends. The standard error is the interval's width divided by
`2 t`, and `variance` is its square. `cdf` and `cdf_se` in the result
hold `F(q)` and its standard error.

Where the interval for `F(q)` runs below 0 or above 1, the sample bounds
the quantile on one side only. That end of `ci` is then `NA`, and so are
`se` and `variance`, with a note.

At the largest sampled value `F(q)` is 1 and its residuals are 0, so its
estimated variance is 0 whatever the spread of `y`, and the interval
would have no width. A quantile there, which a `p` close to 1 or a `y`
with one value gives, gets no standard error, with a note, unless every
frame row the design can reach was taken with certainty, when it is
exact. A quantile whose sampled rows lie in one primary sampling unit
gets none either, as a Hajek mean does, because its residuals sum to
zero over that unit.

A quantile has no design effect, and `deff` is `NA`.

## Domains

`by` gives the quantiles of each domain, each from the distribution
function of the domain's sampled rows. The variance of `F(q)` is that of
the domain's Hajek mean of `y <= q`, from the whole sample, as described
under "Domains" in
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
and the degrees of freedom are the design's. The result is a `drawn_by`
table with a row for each domain and probability.

## Comparison with survey

For the designs
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
maps with identical standard errors, the estimate, the interval and the
standard error equal those of
`survey::svyquantile(ci = TRUE, interval.type = "mean")` on its result,
whose interval is also centred on `F(q)` rather than on `p`. They differ
in these places:

- At the largest sampled value survey gives an interval of no width and
  a standard error of 0, where `ht_quantile()` gives none.

- Where the interval for `F(q)` passes 0 or 1, survey gives `NaN` for
  that end and the standard error, where `ht_quantile()` gives `NA`.

- For a `p`, or an end of the interval for `F(q)`, below the share of
  the weight held by the smallest sampled value, survey's rule gives the
  second smallest value and `ht_quantile()` the smallest.

- For a `p` within `1e-9` above a share of the weight, `ht_quantile()`
  gives the value at that share and survey the next one, as "The
  quantile rule" describes.

- With `by`,
  [`survey::svyby()`](https://rdrr.io/pkg/survey/man/svyby.html) with
  [`survey::svyquantile()`](https://rdrr.io/pkg/survey/man/svyquantile.html)
  forms each interval on the degrees of freedom of the domain's own
  units, and `ht_quantile()` on the design's, as for every domain
  estimate here. With `df = Inf` in both, the standard errors agree.

For systematic, spread and systematic PPS samples the variance of `F(q)`
comes from the estimator this package uses for their totals, which
`survey` does not have.

A missing `y` is refused. `survey::svyquantile(na.rm = TRUE)` on
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
estimates the quantiles over the rows that have one.

## References

Woodruff, R. S. (1952). Confidence intervals for medians and other
position measures. *Journal of the American Statistical Association*,
47, 635–646.
[doi:10.1080/01621459.1952.10483443](https://doi.org/10.1080/01621459.1952.10483443)

Särndal, C.-E., Swensson, B. and Wretman, J. (1992). *Model Assisted
Survey Sampling*. Springer.
[doi:10.1007/978-1-4612-4378-6](https://doi.org/10.1007/978-1-4612-4378-6)

## See also

[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md),
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)

## Examples

``` r
set.seed(1)
pop <- data.frame(
  id = 1:300,
  site = rep(c("a", "b", "c"), each = 100),
  spend = round(stats::rlnorm(300, 4, 0.8))
)
s <- draw(pop, design_stratified("site", n = 60), seed = 1, weights = TRUE)

# The median, and its Woodruff interval
ht_quantile(s, "spend")
#> Quantile of `spend`  (stratified design, n = 60)
#>   prob     0.5
#>   estimate 55
#>   se       6.741689  (analytic)
#>   95% CI  49 to 76  (Woodruff, t, 57 df)
stats::median(pop$spend)   # the truth
#> [1] 53

# The quartiles of each site, from the whole sample
ht_quantile(s, "spend", probs = c(0.25, 0.5, 0.75), by = "site")
#> Quantiles of `spend` by domain  (stratified design, 95% CI, Woodruff, t, 57 df)
#>  site prob  n quantile     se ci_lower ci_upper   method
#>     a 0.25 20       33  9.488       20       58 analytic
#>     a 0.50 20       75 13.733       40       95 analytic
#>     a 0.75 20       95 13.483       86      140 analytic
#>     b 0.25 20       38  5.993       26       50 analytic
#>     b 0.50 20       51  9.239       39       76 analytic
#>     b 0.75 20       76 21.723       53      140 analytic
#>     c 0.25 20       32  8.989       18       54 analytic
#>     c 0.50 20       55 18.977       40      116 analytic
#>     c 0.75 20      116 38.453       62      216 analytic
```
