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
  df = NULL,
  alternative = c("two.sided", "less", "greater"),
  contrast = NULL,
  ratio = NULL
)
```

## Arguments

- sample:

  A data frame returned by
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with
  `weights = TRUE`. See "Which columns are read".

- y:

  The variable to total: a column name, several column names, a numeric
  vector as long as `sample`, or a single number, which stands for that
  value in every row. A factor or character column, or a factor as long
  as `sample`, gives the count of each level. See "Several variables"
  and "Categories".

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

- by:

  Optional column name(s) in `sample` defining domains. When given, the
  result is one row per domain. See "Domains". A column named like one
  of the table's own, such as `n` or `se`, is refused.

- df:

  Degrees of freedom for the interval. `NULL` uses the design's (see
  "Confidence intervals"); `Inf` gives a normal interval; otherwise a
  single number of 1 or more.

- alternative:

  `"two.sided"` for a confidence interval, `"less"` for a one-sided
  upper bound or `"greater"` for a lower bound, each at `level`. See
  "Confidence intervals".

- contrast:

  Optional weights over the `by` domains, as a named numeric vector such
  as `c(a = 1, b = -1)` for domain `a` less domain `b`. A name is the
  domain's value of a single `by` column, or the values of several
  joined by `" . "`, as in `"a . east"`, the way the notes of a `by`
  table name domains. Needs `by`. See "Contrasts".

- ratio:

  Optional name of a numeric column `x` of the frame, for the ratio
  estimator of the total, `sum(y / pi) / sum(x / pi)` times the frame's
  total of `x`. Not with `by`. See "Ratio estimation".

## Value

A list with a [`print()`](https://rdrr.io/r/base/print.html) method,
holding:

- `total`:

  The Horvitz-Thompson total, `sum(y / pi)`.

- `variance`, `se`, `ci`, `level`, `alternative`, `df`:

  Its estimated variance, standard error, confidence interval or
  one-sided bound, and the degrees of freedom the interval used. `NA`
  where the design supports no variance.

- `ci_scale`:

  The scale the interval was formed on: `"linear"` for a total, and
  `"logit"` or `"beta"` for a proportion from
  [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md).

- `variable`:

  The name of `y` when it is a column name, and `NULL` when it is a
  vector.

- `n`, `n_random`, `design`:

  Rows used, those of them not taken with certainty, which the beta
  interval of a proportion counts as its sample size, and the design's
  type.

- `deff`:

  The design effect — see
  [`deff()`](https://elkronos.github.io/dRawn/reference/deff.md).

- `method`:

  Which estimator produced the variance. See "Variance".

- `note`:

  Why a variance is missing, what approximation was used, or which
  fallback was taken. `NULL` when an exact estimator applied cleanly.

- `unreachable`:

  Rows of the frame with inclusion probability 0, which the estimate
  does not cover; `NA` for a sample saved by drawn 0.2.0. See "Rows the
  design cannot reach".

- `unavailable`:

  For a sample from
  [`use_reserves()`](https://elkronos.github.io/dRawn/reference/use_reserves.md)
  only, the number of those rows found unavailable. The estimate takes
  every frame row nobody visited as available; see "What the estimates
  describe" in
  [`use_reserves()`](https://elkronos.github.io/dRawn/reference/use_reserves.md).

With `ratio`, the list also holds `ratio`, the estimated ratio of `y` to
`x`, `x`, the column's name, and `frame_total`, its total over the
frame.

With `by`, a data frame of class `drawn_by` instead: one row per domain,
with the domain columns, `n` (an integer), `total`, `variance`, `se`,
`ci_lower`, `ci_upper` and `method`, and the level, alternative, scale,
degrees of freedom, notes, unreachable count and `variable` as
attributes. It prints without `variance`, which is the square of `se`. A
`by` table has no design effect.

With `contrast`, a list of class `drawn_contrast` with a
[`print()`](https://rdrr.io/r/base/print.html) method, holding
`contrast`, the weighted sum; `weights` and `estimates`, the weight and
estimate of each domain it weighs; `what` (`"total"` or `"mean"`),
`estimator` for a mean, and `by`; and `variable`, `variance`, `se`,
`ci`, `level`, `alternative`, `df`, `n`, `n_random`, `ci_scale`,
`design`, `method`, `note` and `unreachable` as above. `deff` is `NA`.

With several column names as `y`, a `drawn_by` table with a row for each
variable, or each variable and domain, led by a `variable` column; see
"Several variables". With a factor or character `y`, the same table with
a row for each level, led by a `level` column; see "Categories".

Each has [`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
`tidy()` methods; see "Working with the result".

## Variance

Each design gets the estimator the survey-sampling literature recommends
for it, and `method` in the result says which one was used:

- `"analytic"`:

  Exact and design-unbiased. Fixed-size designs with known joint
  probabilities use the Sen-Yates-Grundy estimator. For simple random,
  reservoir and spatial samples, strata and intervals it reduces to the
  textbook formula built from the variance of `y` within each stratum or
  interval, and for an equal-take
  [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
  to the two-stage formula built from the variance between the estimated
  cluster totals and within clusters (Särndal, Swensson and Wretman
  1992). Those formulas are what is computed, in time and memory that
  grow with the number of sampled rows rather than with its square.
  Poisson sampling, whose rows are independent, uses
  `sum((1 - pi) / pi^2 * y^2)`; cluster designs apply the estimator at
  the cluster level, because taking whole clusters makes the row count
  random; certainty designs hand the problem to `rest`, since certainty
  rows are in every sample and add nothing. A
  [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md)
  sample with `n_starts` of 2 or more uses the variance between its
  replicates' totals, the random group estimator (Wolter 2007), on
  `n_starts - 1` degrees of freedom.

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
  stop adding variance they do not have. With `strata` it is added up
  over the strata.
  [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
  with `size` applies it to the estimated totals of the clusters, and
  adds the variance within them.

- `"successive difference"`:

  [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md).
  The design has one random start, and at a step of 2 rows or more
  neighbouring rows are never drawn together, so no design-unbiased
  estimator exists. This uses the successive-difference approximation,
  which compares each sampled row with its neighbour in the order the
  design walked (Wolter 2007), at every step, including a step under 2
  rows, where every pair of rows can be drawn together. It removes a
  trend along that order rather than counting it as noise, which is why
  it beats treating the sample as simple random: on a frame sorted by a
  variable related to `y`, the simple random formula overstated the
  variance two- to thirtyfold in the package's simulations. This one
  came closer, but its average over the random starts depends on the
  frame, and on the estimation article's example it was 1.6 times the
  true variance. No estimator from a single systematic sample escapes
  its limits. A frame that cycles with a period equal to the step, or
  one that divides it, such as 5 rows at a step of 10, puts every row of
  one sample at the same phase of the cycle, and the variance is
  understated, to under 0.001 of the truth in one simulation. A period
  of twice the step puts neighbouring sampled rows at opposite phases,
  and the variance is overstated, 29-fold in the same simulation. When a
  smooth trend so dominates `y` that the random start is almost the only
  source of variation, it runs low, and has no floor. It was 0.3 to 0.7
  of the truth when the trend spanned 25 noise standard deviations, 0.2
  at 100 and 0.06 at 1,000. With `y` a straight line along the walk and
  no noise, it is `6 * k / (n * (k + 1))` times the truth, about
  `6 / n`, for `n` rows sampled at a step of `k`. When `interval` does
  not divide the number of rows walked, the sample size is random, one
  of two values, and the total moves with it in a way successive
  differences cannot see. A term for that,
  `mean(y / pi)^2 * (r / k) * (1 - r / k)` with `k` the interval and `r`
  the remainder, is added, and the note says so. Without it a constant
  `y` had a standard error of 0, and a variable with a coefficient of
  variation of 0.02 an interval that never covered the truth. A design
  given `n` draws exactly `n` rows at a step of `N / n` and needs no
  such term.
  [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md)
  with `within = design_systematic()` adds up this estimator over its
  strata.

- `"local mean"`:

  [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md).
  A spatially balanced sample is well spread precisely because
  neighbours are rarely taken together, so this compares each sampled
  row with the local mean of its nearest sampled neighbours (Stevens and
  Olsen 2003; Grafström and Schelin 2014). With `strata` it runs within
  each stratum, and the strata's variances are added.

- `"jackknife"`:

  The delete-one-unit jackknife, by stratum (JKn) where the design has
  strata, and deleting whole clusters where it has them. Used when
  requested, and as the fallback under `variance = "auto"`.

Rows taken with certainty are in every possible sample and add nothing
to the variance. That covers rows above a
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md)
threshold, rows whose size puts them past certainty in
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
or
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md),
strata taken whole, clusters taken whole by a two-stage design that
selects every cluster, and a census. Every estimator leaves them out.
The jackknife holds them fixed instead of deleting them, and its finite
population correction counts only the other units. A sample made of
nothing else has a variance of 0 under the design's own estimator.

The analytic estimator declines — returning `NA` with the reason in
`note` — rather than understate. Only a decline sends
`variance = "auto"` to the jackknife. Any other error inside an
estimator, such as running out of memory, stops `ht_total()` with that
error. Two cases matter in practice:

- **A stratum or interval with a single sampled row.** There is no
  variation inside it to measure, and simply leaving it out, as the
  Sen-Yates-Grundy sum does, reports a standard error that is too small
  — zero, if every stratum is like that. Draw at least two per stratum
  (`min_per_stratum = 2`, `per_interval = 2`).

- **One row per cluster in a multistage design.** Two rows of one
  cluster are never drawn together, so the within-cluster variation is
  unmeasurable. With `f1 = n_clusters / N_clusters`, the share of
  clusters selected, the exact estimator understates by `f1` times the
  within-cluster component of the variance. `variance = "auto"` falls
  back to the delete-a-cluster jackknife, the ultimate-cluster
  approximation standard in survey practice. It has no first-stage
  finite population correction, so it overstates by `f1` times the
  between-cluster component, the variance the cluster totals would have
  if clusters were drawn with replacement. The excess is small when few
  clusters are selected. With 30 of 40 clusters whose totals differed
  widely, the jackknife averaged 3.3 times the true variance. Whenever
  `f1` is above 0.2, the note gives the factor by which the jackknife
  inflates the between-cluster part. When every cluster is selected, the
  clusters are strata and the variance between them is zero by design,
  so
  [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md)
  with the cluster column as `strata` and two rows in each gives an
  exact variance instead.

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

Units taken with certainty count in neither number. They are in every
sample and carry no variance, and counting them made the interval too
narrow. With eight certain units among ten, systematic PPS reported 9
degrees of freedom where the variance rests on two units, and its
interval covered 79%. A stratum or interval with no other unit drops out
of the count of strata. With no unit below certainty the variance is 0
and `df` is `Inf`.
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
hands certain units to `survey` in a form that makes
[`survey::degf()`](https://rdrr.io/pkg/survey/man/svychisq.html) give
the same count, except for a census, where it gives `n - 1`.

An estimate whose variance is 0 is exact and gets a zero-width interval
whatever `df` is. Otherwise a design that leaves 0 degrees of freedom,
such as Poisson sampling that drew one row, gives a standard error but
no interval, and `note` says why. `df` set by hand must be at least 1. A
proportion of 0 or 1 is the exception to the first rule; see
"Proportions" in
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md).

The t interval assumes the estimate is close to normally distributed
over repeated samples. When `y / pi` is skewed and the sample small it
is not, and the interval covers less often than `level` says, even where
the variance estimate is close to unbiased. In 2,000 draws with `n = 30`
from 240 rows in proportion to a size whose coefficient of variation was
0.66, with `y` close to proportional to size, the 95% interval for the
total covered 90% of the time under systematic PPS and 94% under Poisson
PPS. The mean with `estimator = "ht"`, which
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
recommends in that case, is the total divided by `N` and covered as
often. The default Hajek mean covered 91% under both methods. Skew is
one reason. The Hajek mean also has a ratio bias, here 1.5% of the mean
under systematic PPS and 1.8% under Poisson PPS, and its estimated
variance averaged 0.89 and 0.88 times the variance of the mean over the
draws. Skew lowers coverage under equal probabilities too. A simple
random sample of 30 from 600 lognormal values with a coefficient of
variation of 1.5 covered 87%. A larger sample brings the coverage closer
to `level`. The audit-sampling article checks a bound's coverage by
drawing the design many times, which is how to check it on a frame of
your own.

`alternative = "less"` gives a one-sided upper bound, the figure an
auditor signs off against, and `"greater"` a lower bound. The other end
of `ci` is then `-Inf` or `Inf`, or 0 or 1 for a proportion, and the
printout labels the bound. A one-sided bound at `level` is the end of a
two-sided interval at `1 - 2 * (1 - level)`, so the 95% upper bound is
the upper end of the 90% interval, which `level = 0.90` also gives.

## Rows the design cannot reach

Rows with inclusion probability 0, such as rows past `max_items`,
outside a region or time window, dropped by `na_rm` or in a stratum
allocated no rows, can never be drawn, so a total covers the other rows
only. [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)
warns about some of them by name. The estimate keeps their number in
`unreachable` and prints it on a line of its own, since a sample saved
and read back, or drawn with the warning muffled, would otherwise give
no sign of them. The design effect is measured against a simple random
sample of the rows the design can reach; see
[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md).

## An empty sample

A design of random size, such as Poisson sampling or a systematic sample
whose random start falls past the last row, can draw nothing.
`ht_total()` then returns a total of 0, which keeps the estimator
unbiased over repeated draws, with `df` 0, the variance, standard error
and interval `NA`, and a note that the sample is empty.
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
refuses the Hajek mean of an empty sample, which is 0 / 0;
`estimator = "ht"` gives 0. An empty sample from a design that gives
every row a probability of 0, such as one with `n = 0`, is refused, and
so is any empty sample in
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md).

## Domains

`by` estimates for each level of one or more columns — sites, regions,
months — without subsetting the sample first. Subsetting a sample and
treating the piece as its own sample gets the variance wrong, because
the number of rows that happen to fall in a domain is itself random.
Here each domain's total is estimated from `y * (row in domain)` over
the whole sample, so the same design and the same variance machinery
apply (Särndal, Swensson and Wretman 1992, section 10.3).

Domains come out in key order: a factor's levels, numbers and dates in
sequence, and text sorted byte by byte as in the C locale, the same in
every session. A domain with no sampled row is left out, as
[`survey::svyby()`](https://rdrr.io/pkg/survey/man/svyby.html) leaves it
out.

A stratum or interval with a single sampled row stops a domain's
variance, as it stops the whole sample's, unless it holds no unit of the
domain. That is read from the `by` columns of the frame stored with the
sample, so it needs them in the frame, with the same values as in the
sample. A domain made of strata taken whole is then exact, with a
zero-width interval. A lonely stratum whose one sampled row lies outside
the domain still stops it when other units of the stratum belong to the
domain, because their variation went unmeasured. A frame row whose `by`
value is missing may belong to any domain, so a lonely stratum holding
one stops every domain.

A domain's Hajek mean
([`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md))
needs sampled rows in at least two primary sampling units, rows or
clusters. Its linearised residuals sum to zero inside the unit that
holds the whole domain, so with one unit every estimator gives a
standard error of about 0, and `survey` does the same. Such a domain
gets no standard error, with a note. When stage one of a cluster design
takes every cluster, the analytic estimator's units are rows, since each
cluster is then in every sample and only its rows were drawn at random.
The jackknife still deletes whole clusters, so a domain inside one
cluster gets no standard error from it, whether `variance = "jackknife"`
asks for it or `"auto"` turns to it after the analytic estimator
declines. Rows taken with certainty count as one unit, since every
estimator holds them fixed and
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
hands them to `survey` together, so a domain whose sampled rows were all
taken with certainty gets no standard error either, unless every row of
it that the design could reach is certain, which makes its mean exact. A
domain total, and the `estimator = "ht"` mean, keep the standard error
of 0 their estimators give such a domain, since those estimators are
unbiased over samples, and a note says the interval has no width because
the sampled rows were all certain. A domain with a handful of sampled
rows has a standard error that is itself unreliable, and the interval
covers less often than `level` says. In one simulation of a skewed
variable, a 95% interval for the mean of a domain with two or three
sampled rows covered 65% of the time, and 85% with four or five.

A sample grouped with `dplyr::group_by()` is estimated as a whole, with
a warning of class `drawn_grouped_df`; the groups play no part. Use
`by`.

The strata, clusters and other pieces the variance reads from the frame
are worked out once for all the domains. For simple random, reservoir,
spatial, stratified and temporal samples the analytic variances of all
the domains come from one pass over the sample, in time that grows with
the number of sampled rows rather than with rows times domains. A domain
with units in a stratum or interval of one sampled row, and every domain
under `variance = "jackknife"`, is estimated on its own.

## Contrasts

`contrast` weighs the domain estimates of `by` and gives their weighted
sum, such as the difference between two sites' totals or means, with a
standard error, an interval and `df` as for any other estimate. Each
domain's estimate is the total, over the whole sample, of a variable
that is 0 outside the domain: `y` for a total, the linearised residual
for a Hajek mean and `y / N_d` for a Horvitz-Thompson mean. The weighted
sum of the estimates is the total of the same weighted sum of those
variables, so one run of the design's variance estimator gives its
variance, with the covariances between domains included. For simple
random, stratified, cluster and multistage samples the standard error
equals that of
[`survey::svycontrast()`](https://rdrr.io/pkg/survey/man/svycontrast.html)
applied to `survey::svyby(..., covmat = TRUE)`.

The interval is symmetric, for a proportion too, since a difference of
proportions can be negative. Domains left out of `contrast` have a
weight of 0. A name that is not a domain with a sampled row is refused,
and so is a contrast whose weights are all 0. When a Hajek mean that the
contrast weighs rests on one primary sampling unit, the contrast gets no
standard error, for the reason given under "Domains".

## Ratio estimation

`ratio` names a column `x` of the frame, and the total of `y` is then
estimated by the ratio estimator `sum(y / pi) / sum(x / pi) * X`, where
`X` is the total of `x` over the frame rows the design can reach, read
from the frame stored with the sample, as
[`frame_totals()`](https://elkronos.github.io/dRawn/reference/frame_totals.md)
reads it. When `y` is close to proportional to `x`, such as an audited
value and the book value it was drawn on, the residuals `y - R x` vary
far less than `y`, and the estimate is more precise than the
Horvitz-Thompson total. Its variance is `X^2` times that of the ratio
from
[`ht_ratio()`](https://elkronos.github.io/dRawn/reference/ht_ratio.md),
estimated from the residuals `X (y - R x) / X_hat` with the design's own
estimator, and it equals
`predict(survey::svyratio(~y, ~x, design), total = X)` for the designs
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
maps with identical standard errors. The sampled rows' values of `x`
come from the frame too. A value of `x` missing from a frame row the
design can reach is refused, since `X` is then unknown, and so is
`ratio` with `by`. The result holds the ratio as `ratio`, the column as
`x` and `X` as `frame_total`, and prints them.

The ratio estimator is biased, by an amount of order `1 / n`, and its
design effect compares it with the same estimator from a simple random
sample of the same size. Calibration in `survey`, with the totals from
[`frame_totals()`](https://elkronos.github.io/dRawn/reference/frame_totals.md),
covers several auxiliary columns at once.

## Several variables

`y` may name several columns, such as `c("spend", "visits")`. The result
is then a `drawn_by` table with a row for each variable, or for each
variable and domain with `by`, led by a `variable` column. Each row
holds what a call for that variable alone gives: the same estimate,
standard error, interval and method, with the same refusals for a
missing or infinite value, named by its column. Notes that hold for
every variable are printed once, and the others after the variable's
name.

[`vcov()`](https://rdrr.io/r/stats/vcov.html) of the table gives the
covariances between the variables' estimates as well as their variances.
Every variance estimator here is a quadratic form in the variable it
totals, so the covariance of two estimates whose variables are `z_a` and
`z_b` is `(V(z_a + z_b) - V(z_a) - V(z_b)) / 2`, one more run of the
same estimator. For a total these variables are the `y` columns, and for
a Hajek mean their linearised residuals. For the designs
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
maps with identical standard errors, the matrix equals
`vcov(survey::svytotal(~a + b, design))`, or that of
[`survey::svymean()`](https://rdrr.io/pkg/survey/man/surveysummary.html)
for means. With `by`, the covariances between the variables within each
domain are estimated, and those between domains are not, as for a `by`
table of one variable. A covariance is `NA` when either variance is
missing, or when the two variances, or the variance of their sum, did
not come from one estimator, as when the design's own declined for one
and the jackknife stood in.

A `by` column named `variable` is refused, and so is `contrast`. The
table has no design effect.

## Categories

A factor or character column named as `y`, or a factor as long as
`sample`, is split into a 0/1 indicator for each of its levels:
`ht_total()` gives the estimated number of frame rows at each level, and
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md) the
share of each, with a proportion's interval. The result is the table of
"Several variables" with a row for each level, or each level and domain
with `by`, led by a `level` column in place of `variable`, and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) holds the covariances
between levels. The levels are a factor's levels, unused ones included,
and the values of a character column sorted byte by byte. A missing
value is a missing `y`. For the designs
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
maps with identical standard errors, the estimates and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) equal those of
`survey::svytotal(~f, design)` and `survey::svymean(~f, design)`. A
character vector given as `y` is read as column names, so pass
categories held in a vector as a factor. A `by` column named `level` is
refused, and so is `contrast`; the indicator of one level, such as
`sample$f == "a"`, takes it.

A single number as `y` stands for that value in every row:
`ht_total(sample, 1)` estimates the number of frame rows the design can
reach, and with `by` the number in each domain.

## Working with the result

[`coef()`](https://rdrr.io/r/stats/coef.html) gives the estimate, named
after `y` when `y` is a column name,
[`vcov()`](https://rdrr.io/r/stats/vcov.html) its variance as a 1 by 1
matrix, and [`confint()`](https://rdrr.io/r/stats/confint.html) its
interval as a matrix with a column for each end. At the estimate's own
`level`, [`confint()`](https://rdrr.io/r/stats/confint.html) returns the
interval the estimate holds. At another level it forms the interval
again on the design's degrees of freedom and on the estimate's scale
(`ci_scale`), so it gives what `level =` in `ht_total()` or
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
would have given.
[`survey::SE()`](https://rdrr.io/pkg/survey/man/SE.html) and
[`survey::cv()`](https://rdrr.io/pkg/survey/man/surveysummary.html) work
through [`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html).

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives one
row of `variable`, `n`, the estimate, `variance`, `se`, `ci_lower`,
`ci_upper`, `level`, `df`, `deff` and `method`, so estimates of several
variables can be bound with
[`rbind()`](https://rdrr.io/r/base/cbind.html). With the generics
package installed,
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html),
and so `broom::tidy()`, gives `term`, `estimate` and `std.error`, with
`conf.low` and `conf.high` when `conf.int = TRUE`, at `conf.level`,
which is the estimate's own level unless given.

A `by` table has the same methods.
[`coef()`](https://rdrr.io/r/stats/coef.html) names each domain as
`contrast` names it, [`confint()`](https://rdrr.io/r/stats/confint.html)
and `tidy()` give a row per domain, `tidy()` with the `by` columns, and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) drops the
table's class and attributes.
[`vcov()`](https://rdrr.io/r/stats/vcov.html) holds the domains'
variances. The covariances between domains are not estimated, and its
entries for them are `NA`, except 0 for two domains of a
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md)
sample whose sampled rows share no stratum, since strata are sampled
independently. The standard error of a difference between domains comes
from `contrast`, not from [`vcov()`](https://rdrr.io/r/stats/vcov.html).
A table of several variables names each row by its variable, then its
domain, as in `"spend . a"`, and its
[`vcov()`](https://rdrr.io/r/stats/vcov.html) holds the covariances
between variables described under "Several variables". Selecting columns
of the table with `[` or
[`subset()`](https://rdrr.io/r/base/subset.html) drops the attributes
these methods read, and the result prints as a plain data frame.
Selecting rows keeps them. A `drawn_contrast` has the same methods, with
the estimate named by its weights, as in `"a - b"`, and so does a ratio
from
[`ht_ratio()`](https://elkronos.github.io/dRawn/reference/ht_ratio.md),
named as `"y/x"`.

## Missing values and nonresponse

`ht_total()` refuses a missing `y`, from a nonresponse or a measurement
that failed. A total over the rows with `y` present counts the others as
0, so it is low by about the share of the weight they hold, the weighted
nonresponse rate, and the message gives that share.

`ht_mean(na_rm = TRUE)` estimates the mean over the rows with `y`
present. They are a domain of the whole sample, as a `by` domain is: the
estimate is the Hajek mean over them, and its variance comes from the
linearised residuals over every sampled row, 0 for the rows left out,
with the design's own estimator and degrees of freedom.
`survey::svymean(na.rm = TRUE)` on
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
gives the same estimate and standard error. With `by`, each domain's
mean is over its rows with `y` present, and a domain whose sampled rows
all have `y` missing has no mean, no row in the table and no place in
`contrast`; the note names it.
[`add_fieldwork()`](https://elkronos.github.io/dRawn/reference/add_fieldwork.md)
joins results keyed by an id and leaves `y` missing for a sampled row
with none.

That mean is unbiased for the frame's mean only if the values are
missing completely at random, with no relation between whether a row has
a `y` and its `y`. When response rates differ between groups whose `y`
differs, such as sites, the rows with `y` present over-represent the
groups that responded more. A weighting-class adjustment corrects for
that, under the weaker assumption that response is unrelated to `y`
within each group. Multiply the `.weight` of each row with `y` present
by its group's sum of `.weight` over all its sampled rows, divided by
the same sum over its rows with `y` present, and estimate totals and
means from
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
with `survey`'s `na.rm = TRUE`. `survey` treats the adjusted weights as
fixed, so its standard errors leave out that the adjustment factors were
estimated from the sample. `ht_total()` and
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
keep to `.prob` and do not see the adjustment. The estimation article
shows the steps.

A unit that cannot be reached at all, rather than one that does not
answer, can be replaced from a reserve list drawn with the sample:
`reserve` in
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md)
and
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md),
then
[`use_reserves()`](https://elkronos.github.io/dRawn/reference/use_reserves.md).
The estimates then describe the units that can be reached.

## Which columns are read

`ht_total()` and
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
divide by `.prob`, the inclusion probabilities
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) attached,
find each sampled row's frame row through `.row`, and read strata,
clusters, time stamps and the other design columns from the copy of the
frame stored with the sample. Of the sample's own columns they read
`.prob`, `.row`, `y` and `by`, and `.weight` only to warn when it was
edited, so renaming, recoding or dropping a design column in the sample
changes nothing.

`.weight` is not used in the estimate. A weight adjusted after the draw,
for nonresponse or by trimming, and written into `.weight` is ignored
with a warning of class `drawn_weight_edited`;
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
hands `.weight` to `survey`, which is where adjusted weights are
estimated from. A sample whose `.prob` was changed or removed is
refused. So is one whose design gives the frame different inclusion
probabilities in the current session than it did at draw time, which a
change of
[`sf::sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html)
does to a
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md)
sample: the variance is recomputed from the design and would describe
another one. A sample saved by drawn 0.2.0 has no `.row` and records no
probabilities for the frame, so its rows are matched by position and it
is refused when the design gives a sampled row a probability other than
its `.prob`.

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
[`ht_ratio()`](https://elkronos.github.io/dRawn/reference/ht_ratio.md),
[`ht_quantile()`](https://elkronos.github.io/dRawn/reference/ht_quantile.md),
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

# How much more site a spent than site b
ht_total(s, "spend", by = "site", contrast = c(a = 1, b = -1))
#> Contrast of Horvitz-Thompson totals  (stratified design, n = 40)
#>   contrast `a` - `b`
#>   domains  `a` 40,715, `b` 13,285
#>   estimate 27,430
#>   se       3,713.698  (analytic)
#>   95% CI  19,912.01 to 34,947.99  (t, 38 df)

# Two variables at once, with the covariance between their totals
s$visits <- round(s$spend / 50) + 1
both <- ht_total(s, c("spend", "visits"))
both
#> Horvitz-Thompson total by variable  (stratified design, 95% CI, t, 38 df)
#>  variable  n total     se ci_lower ci_upper   method
#>     spend 40 54000 3713.7    46482    61518 analytic
#>    visits 40  1285   75.1     1133     1437 analytic
vcov(both)
#>             spend    visits
#> spend  13791552.3 277139.69
#> visits   277139.7   5640.23
```
