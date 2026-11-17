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
  level = 0.95,
  by = NULL,
  df = NULL,
  proportion = NULL,
  ci_method = c("logit", "beta"),
  alternative = c("two.sided", "less", "greater"),
  contrast = NULL,
  na_rm = FALSE
)
```

## Arguments

- sample:

  A data frame returned by
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with
  `weights = TRUE`. It is read as
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  reads it: `.prob` for the probabilities, and the stored frame for the
  design's columns. See "Which columns are read" there.

- y:

  The variable to average: a column name, several column names, a
  numeric vector as long as `sample`, or a single number, which stands
  for that value in every row. A factor or character column, or a factor
  as long as `sample`, gives the share of each level. See "Several
  variables" and "Categories" in
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md).

- estimator:

  `"hajek"` or `"ht"`. See above.

- variance, level, by, df, alternative, contrast:

  As for
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md).

- proportion:

  Whether `y` is a proportion, which decides the shape of the interval.
  `NULL` decides from `y`: a logical `y`, or a column of the frame whose
  values are all 0 or 1, is a proportion, and a numeric vector or a
  column added after the draw is one when it is 0 or 1 in every sampled
  row. `TRUE` treats `y` as one, and needs it to be 0 or 1 in every
  sampled row; `FALSE` gives the symmetric interval. See "Proportions".

- ci_method:

  The interval for a proportion: `"logit"`, or `"beta"` for Korn and
  Graubard's. A proportion of exactly 0 or 1 gets `"beta"` either way.
  See "Proportions".

- na_rm:

  `FALSE` refuses a missing `y`. `TRUE` estimates the mean over the rows
  with `y` present, as a domain of the sample. See "Missing values".

## Value

A list with a [`print()`](https://rdrr.io/r/base/print.html) method,
holding `mean`, `estimator`, `variable`, `variance`, `se`, `ci`,
`level`, `alternative`, `df`, `n`, `n_random`, `design`, `deff`,
`method`, `note`, `ci_scale` (`"linear"`, `"logit"` or `"beta"`),
`unreachable` and `n_missing`, the sampled rows with `y` missing that
`na_rm = TRUE` left out. See
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
for what `variable`, `n_random`, `method`, `df`, `note` and
`unreachable` say, and
[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md) for the
design effect. With `by`, a data frame with one row per domain, and with
`contrast` a `drawn_contrast`, as for
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md).
With several column names as `y`, a `drawn_by` table with a row per
variable, and with a factor or character `y` one with a row per level,
as for
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md).
Each has [`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
`tidy()` methods; see "Working with the result" in
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md).

## Proportions

For a 0/1 or logical variable either one estimates a proportion, and the
confidence interval is then formed on the logit scale and transformed
back (Korn and Graubard 1999, section 5.3). A symmetric interval around
a small proportion runs below zero; this one cannot, and matches
`survey::svyciprop(method = "xlogit")`.

Whether `y` is a proportion is read from what `y` is, not from the
values the sample happened to draw. A logical `y` is one, and so is a
column of the frame whose values are all 0 or 1. A count that came out 0
or 1 in every sampled row is not: a logit interval cannot pass 1, so it
would never cover a mean of 1.2. A numeric vector, or a column added to
the sample after the draw, such as a 0/1 result recorded in fieldwork,
has no column of the frame to read. It is taken for a proportion when it
is 0 or 1 in every sampled row, and `note` says the choice came from the
sample. `proportion = FALSE` asks for the symmetric interval whatever
`y` is, and is the choice for a count.

`ci_method = "beta"` gives the Korn and Graubard (1998) interval
instead, a Clopper-Pearson interval at an effective sample size of
`p(1 - p) / v`, with `p` the estimate and `v` its variance, capped at
`n` and then scaled by `(qt(a, n - 1) / qt(a, df))^2` for the design's
degrees of freedom. Here `n` is the number of sampled rows below
certainty behind the estimate and `a` the tail probability of the
interval. The scaled size is capped at `n` too, so a domain whose design
has more degrees of freedom than the domain has rows gets no narrower an
interval than Clopper-Pearson's over those rows. With no rows taken with
certainty it matches `survey::svyciprop(method = "beta")` whenever
neither cap binds. It ignores the finite population correction, so it
errs wide when the sample is a large share of the frame.

An estimate of exactly 0 or 1, such as no deviations found in 60 sampled
invoices, has a variance of 0, and a logit interval around it has no
width. It gets the beta interval whatever `ci_method` says, with the
effective sample size set to `n`, and `note` says so. For 0 of 60 from a
simple random sample, the 95% upper bound (`alternative = "less"`) is
4.87%. The interval needs no variance, so a design that leaves 0 degrees
of freedom, such as a Poisson sample of one row, gets Clopper-Pearson's
over the sampled rows. When every sampled row behind the estimate was
taken with certainty, the estimate is exact and the interval has no
width.

## When the two differ, and which to use

They coincide **exactly** whenever the design weights of the rows you
drew sum to `N` — which covers every fixed-size equal-probability design
that can reach every row of the frame, and every
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md)
without replacement that draws from every stratum, because each stratum
contributes `n_h * N_h / n_h = N_h`. There is nothing to choose between
them there. A
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md)
sample from a region, or a
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md)
sample whose stream ran past `max_items`, is a simple random sample of
the rows it can reach, and its weights sum to their number instead; see
the caveat on `"ht"` below.

They differ when that sum is random or uneven:

- the sample size itself is random —
  [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
  with `method = "poisson"`,
  [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)
  over clusters of unequal size, or
  [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md)
  with an `interval` that does not divide `N`;

- the weights vary within a fixed-size sample —
  probability-proportional-to- size selection.

Hajek is usually the steadier of the two and is the default for that
reason (Särndal, Swensson and Wretman 1992, section 5.7). The exception
is worth knowing: when `y` is close to proportional to the size measure
that drove selection, `y / pi` is nearly constant, the Horvitz-Thompson
numerator barely moves, and dividing it by the known `N` beats dividing
by an estimate. The Hajek mean's ratio bias, and a variance estimate
that runs 11 to 12% low, also lower the coverage of its interval there.
Drawing with `n = 30` from 240 rows by Poisson PPS, its 95% interval
covered 91% of the time and that of `"ht"` 94%. By systematic PPS both
covered about 90%. "Confidence intervals" in
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
gives the figures and the reasons.

One caveat on `"ht"`. It is unbiased for the frame mean *provided every
row could have been selected*. Rows with inclusion probability 0, such
as rows outside a time window or a region, past `max_items` or dropped
under `na_rm = TRUE`, sit inside the `N` it divides by but can never
enter the numerator, so the estimate is biased low by exactly their
share of the frame.
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
reports how many such rows there are. The Hajek mean is unaffected,
because it estimates the mean of the part of the frame the design can
actually reach.

## Variance and domains

The Hajek mean is a ratio, and its variance is the variance of the total
of the linearised residuals `(y - mean) / N_hat` (Deville 1999),
computed with whichever estimator suits the design — see
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
which also explains the confidence interval's degrees of freedom, how
`by` estimates domain means from the whole sample rather than a subset
of it, and what an empty sample gives.

The residuals sum to zero inside a primary sampling unit that holds
every sampled row of the domain, so a Hajek mean over one unit, such as
a domain with one sampled row or a Poisson sample of one row, gets no
standard error. The "Domains" section of
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
says more.

With `estimator = "ht"` each domain's total is divided by the domain's
number of rows in the frame, read from the frame's `by` columns. A frame
row with a `by` value missing could be in any domain, so such a frame is
refused.

`contrast` gives a weighted sum of domain means, such as the difference
between two, with its standard error; see "Contrasts" in
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md).
Its interval is symmetric whatever `proportion` and `ci_method` say.

When `y` has the same value in every sampled row, the residuals are 0
and the Hajek mean's variance is set to exactly 0, with a note, rather
than left to rounding error, which gave standard errors near `1e-16` and
design effects that changed from seed to seed. The design effect is then
`NA`.

## Shares and counts of categories

A factor or character column named as `y`, or a factor as long as
`sample`, gives the share of each of its levels: a `drawn_by` table led
by a `level` column, with a proportion's interval for each level and the
covariances between the shares in
[`vcov()`](https://rdrr.io/r/stats/vcov.html).
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
gives the estimated number of frame rows at each level the same way; see
"Categories" there. The share of one level is also the mean of its
indicator, `ht_mean(s, s$habitat == "forest")`, which takes `contrast`
as any `y` does.

A single number as `y` stands for that value in every row, so
`ht_total(s, 1)` estimates the number of frame rows the design can
reach, and `ht_total(s, 1, by = "habitat")` the number in each domain.

## Missing values

A missing `y` is refused unless `na_rm = TRUE`. Then the mean is the
Hajek mean over the rows with `y` present, treated as a domain of the
whole sample, with the design's own variance estimator, and equals
`survey::svymean(na.rm = TRUE)` on
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md).
It assumes the values are missing completely at random, and the printed
result says so. `n` counts the rows with `y` present, `n_missing` the
others, and the design effect is `NA`, as for a domain. With `by`, each
domain is the rows of its `by` values with `y` present.
`estimator = "ht"` is refused, since the number of frame rows that would
have a `y` is unknown. See "Missing values and nonresponse" in
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
for when the assumption fails and what to do then.

## References

Korn, E. L. and Graubard, B. I. (1998). Confidence intervals for
proportions with small expected number of positive counts estimated from
survey data. *Survey Methodology*, 24, 193–201.

Korn, E. L. and Graubard, B. I. (1999). *Analysis of Health Surveys*.
Wiley.

Hájek, J. (1971). Comment on "An essay on the logical foundations of
survey sampling, part one" by D. Basu. In V. P. Godambe and D. A. Sprott
(eds.), *Foundations of Statistical Inference*, p. 236. Holt, Rinehart
and Winston.

Särndal, C.-E., Swensson, B. and Wretman, J. (1992). *Model Assisted
Survey Sampling*. Springer.

Deville, J.-C. (1999). Variance estimation for complex statistics and
estimators: linearization and residual techniques. *Survey Methodology*,
25, 193–203.

## See also

[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`ht_ratio()`](https://elkronos.github.io/dRawn/reference/ht_ratio.md),
[`ht_quantile()`](https://elkronos.github.io/dRawn/reference/ht_quantile.md),
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
#>   95% CI  232.4101 to 307.5899  (t, 38 df)
#>   deff     1.03  (about the same as simple random sampling)
mean(pop$spend)   # the truth
#> [1] 263.66

# A mean for each site, estimated from the whole sample
ht_mean(s, "spend", by = "site")
#> Hajek mean by domain  (stratified design, 95% CI, t, 38 df)
#>  site  n  mean    se ci_lower ci_upper   method
#>     a 30 271.4 21.01    228.9    314.0 analytic
#>     b 10 265.7 39.27    186.2    345.2 analytic

# The difference between the two sites' means
ht_mean(s, "spend", by = "site", contrast = c(a = 1, b = -1))
#> Contrast of Hajek means  (stratified design, n = 40)
#>   contrast `a` - `b`
#>   domains  `a` 271.4333, `b` 265.7
#>   estimate 5.733333
#>   se       44.53818  (analytic)
#>   95% CI  -84.42949 to 95.89616  (t, 38 df)

# A proportion is the mean of a 0/1 variable
ht_mean(s, s$spend > 250)
#> Hajek mean  (stratified design, n = 40)
#>   estimate 0.5
#>   se       0.07207914  (analytic)
#>   95% CI  0.3580894 to 0.6419106  (logit, t, 38 df)
#>   deff     1.01  (about the same as simple random sampling)

# The share of each level of a factor, with the covariances between them
s$band <- cut(s$spend, c(0, 100, 250, 500))
shares <- ht_mean(s, "band")
shares
#> Hajek mean of `band` by level  (stratified design, 95% CI, logit, t, 38 df)
#>      level  n  mean      se ci_lower ci_upper   method
#>    (0,100] 40 0.075 0.03828  0.02585   0.1986 analytic
#>  (100,250] 40 0.425 0.07151  0.29016   0.5720 analytic
#>  (250,500] 40 0.500 0.07208  0.35809   0.6419 analytic
vcov(shares)
#>                 (0,100]     (100,250]     (250,500]
#> (0,100]    0.0014655172 -0.0006915709 -0.0007739464
#> (100,250] -0.0006915709  0.0051130268 -0.0044214559
#> (250,500] -0.0007739464 -0.0044214559  0.0051954023

# Stratified without replacement: the weights sum to N, so the two agree
all.equal(ht_mean(s, "spend", estimator = "ht")$mean, ht_mean(s, "spend")$mean)
#> [1] TRUE

# Poisson sampling has a random size, so they part company
p <- draw(pop, design_weighted("spend", n = 40, method = "poisson"),
          seed = 2, weights = TRUE)
c(hajek = ht_mean(p, "spend", variance = "none")$mean,
  ht    = ht_mean(p, "spend", "ht", variance = "none")$mean)
#>    hajek       ht 
#> 330.0012 329.5750 
```
