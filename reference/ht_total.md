# Estimate a population total from a sample

Forms the Horvitz-Thompson total `sum(y / pi)` with a standard error, a
confidence interval and a design effect, using a variance estimator
matched to how the design actually randomised.

## Usage

``` r
ht_total(
  sample,
  y,
  variance = c("auto", "analytic", "jackknife", "none"),
  level = 0.95,
  by = NULL,
  df = NULL
)
```

## Arguments

- sample:

  A data frame returned by
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with
  `weights = TRUE`.

- y:

  The variable to total: a column name, or a numeric vector as long as
  `sample`.

- variance:

  How to compute it. `"auto"` uses the design's own estimator (see
  "Variance") and falls back to the jackknife when that declines;
  `"analytic"` insists on the design's own estimator, returning `NA`
  with the reason in `note` rather than falling back; `"jackknife"`
  always resamples; `"none"` skips it. The result reports which was used
  in `method`, and `"none"` when nothing could produce a figure.

- level:

  Confidence level for the interval.

- by:

  Optional column name(s) in `sample` defining domains. When given, the
  result is one row per domain. See "Domains".

- df:

  Degrees of freedom for the interval. `NULL` uses the design's (see
  "Confidence intervals"); `Inf` gives a normal interval.

## Value

A list with a [`print()`](https://rdrr.io/r/base/print.html) method,
holding:

- `total`:

  The Horvitz-Thompson total, `sum(y / pi)`.

- `variance`, `se`, `ci`, `level`, `df`:

  Its estimated variance, standard error, confidence interval, and the
  degrees of freedom the interval used. `NA` where the design supports
  no variance.

- `n`, `design`:

  Rows used, and the design's type.

- `deff`:

  The design effect — see
  [`deff()`](https://elkronos.github.io/dRawn/reference/deff.md).

- `method`:

  Which estimator produced the variance. See "Variance".

- `note`:

  Why a variance is missing, what approximation was used, or which
  fallback was taken. `NULL` when an exact estimator applied cleanly.

With `by`, a data frame of class `drawn_by` instead: one row per domain,
with the domain columns, `n`, `total`, `se`, `ci_lower`, `ci_upper` and
`method`.

## Variance

Each design gets the estimator the survey-sampling literature recommends
for it, and `method` in the result says which one was used:

- `"analytic"`:

  Exact and design-unbiased. Fixed-size designs with known joint
  probabilities use the Sen-Yates-Grundy estimator; Poisson sampling,
  whose rows are independent, uses `sum((1 - pi) / pi^2 * y^2)`; cluster
  designs apply the estimator at the cluster level, because taking whole
  clusters makes the row count random; certainty designs hand the
  problem to `rest`, since certainty rows are in every sample and add
  nothing.

- `"deville"`:

  Systematic probability-proportional-to-size
  ([`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
  with `method = "systematic"`). Its joint probabilities have no closed
  form, so this uses Deville's (1999) approximation, which needs only
  the first-order probabilities and is the best-performing of the
  standard approximations for high-entropy designs (Matei and Tillé
  2005). Randomised systematic selection is close to high entropy
  (Hartley and Rao 1962). Unlike a plain jackknife, each row's
  contribution is scaled by its own `1 - pi`, so rows near certainty
  stop adding variance they do not have.

- `"successive difference"`:

  [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md).
  No design-unbiased estimator exists — the design has one random start,
  and most pairs of rows can never be drawn together — so this uses the
  successive-difference approximation, which compares each sampled row
  with its neighbour in the order the design walked (Wolter 2007). It
  removes a trend along that order rather than counting it as noise,
  which is why it beats treating the sample as simple random: on a frame
  sorted by a variable related to `y`, the simple random formula
  overstated the variance two- to thirtyfold in the package's
  simulations, where this one stayed near the truth. Two limits no
  estimator from a single systematic sample can escape: it is
  understated if the frame cycles with a period matching `interval`, and
  when a smooth trend so dominates `y` that the random start is almost
  the only source of variation, it runs low (0.3 to 0.7 of the truth
  when the trend spanned 25 noise standard deviations).

- `"local mean"`:

  [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md).
  A spatially balanced sample is well spread precisely because
  neighbours are rarely taken together, so this compares each sampled
  row with the local mean of its nearest sampled neighbours (Stevens and
  Olsen 2003; Grafström and Schelin 2014).

- `"jackknife"`:

  The delete-one-unit jackknife, by stratum (JKn) where the design has
  strata, and deleting whole clusters where it has them. Used when
  requested, and as the fallback under `variance = "auto"`.

The analytic estimator declines — returning `NA` with the reason in
`note` — rather than understate. Two cases matter in practice:

- **A stratum or interval with a single sampled row.** There is no
  variation inside it to measure, and simply leaving it out, as the
  Sen-Yates-Grundy sum does, reports a standard error that is too small
  — zero, if every stratum is like that. Draw at least two per stratum
  (`min_per_stratum = 2`, `per_interval = 2`).

- **One row per cluster in a multistage design.** The within-cluster
  variation is unmeasurable and the exact estimator understates by
  around a third. `variance = "auto"` falls back to the delete-a-cluster
  jackknife, which is the ultimate-cluster approximation standard in
  survey practice.

A Sen-Yates-Grundy estimate that comes out negative is treated the same
way: a failure of the estimator on an unlucky sample, not a variance.

## Confidence intervals

Intervals use Student's t with the design's degrees of freedom — the
number of primary sampling units minus the number of strata — rather
than the normal distribution (Korn and Graubard 1999). The difference is
negligible with hundreds of units and decisive with a handful: at four
clusters a normal interval covers about 87% where it claims 95%, and the
t interval restores most of that. The degrees of freedom used are
reported as `df`; pass `df = Inf` for the normal interval.

## Domains

`by` estimates for each level of one or more columns — sites, regions,
months — without subsetting the sample first. Subsetting a sample and
treating the piece as its own sample gets the variance wrong, because
the number of rows that happen to fall in a domain is itself random.
Here each domain's total is estimated from `y * (row in domain)` over
the whole sample, so the same design and the same variance machinery
apply (Särndal, Swensson and Wretman 1992, section 10.3).

## References

Horvitz, D. G. and Thompson, D. J. (1952). A generalization of sampling
without replacement from a finite universe. *Journal of the American
Statistical Association*, 47, 663–685.

Sen, A. R. (1953). On the estimate of the variance in sampling with
varying probabilities. *Journal of the Indian Society of Agricultural
Statistics*, 5, 119–127.

Yates, F. and Grundy, P. M. (1953). Selection without replacement from
within strata with probability proportional to size. *Journal of the
Royal Statistical Society B*, 15, 253–261.

Hartley, H. O. and Rao, J. N. K. (1962). Sampling with unequal
probabilities and without replacement. *Annals of Mathematical
Statistics*, 33, 350–374.

Deville, J.-C. (1999). Variance estimation for complex statistics and
estimators: linearization and residual techniques. *Survey Methodology*,
25, 193–203.

Matei, A. and Tillé, Y. (2005). Evaluation of variance approximations
and estimators in maximum entropy sampling with unequal probability and
fixed sample size. *Journal of Official Statistics*, 21, 543–570.

Wolter, K. M. (2007). *Introduction to Variance Estimation*, 2nd ed.
Springer.

Korn, E. L. and Graubard, B. I. (1999). *Analysis of Health Surveys*.
Wiley.

Särndal, C.-E., Swensson, B. and Wretman, J. (1992). *Model Assisted
Survey Sampling*. Springer.

## See also

[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md),
[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md),
[`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md),
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)

## Examples

``` r
set.seed(1)
pop <- data.frame(
  id = 1:200,
  site = rep(c("a", "b"), times = c(150, 50)),
  spend = round(stats::runif(200, 10, 500))
)

s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
ht_total(s, "spend")
#> Horvitz-Thompson total  (stratified design, n = 40)
#>   estimate 54,000
#>   se       3,713.698  (analytic)
#>   95% CI  46,482.01 to 61,517.99  (t, 38 df)
#>   deff     1.03  (about the same as simple random sampling)

sum(pop$spend)   # the truth
#> [1] 52732

# One estimate per site, from the same sample
ht_total(s, "spend", by = "site")
#> Horvitz-Thompson total by domain  (stratified design, 95% CI, t, 38 df)
#>  site  n total   se ci_lower ci_upper   method
#>     a 30 40715 3152    34334    47096 analytic
#>     b 10 13285 1963     9310    17260 analytic
tapply(pop$spend, pop$site, sum)
#>     a     b 
#> 39123 13609 
```
