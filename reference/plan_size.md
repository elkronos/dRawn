# How large a sample do you need?

Solves for `n` given the precision you want, rather than asking you to
guess it. This is the step before
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md): you
rarely know the sample size, you know the margin of error you can live
with.

## Usage

``` r
plan_size(
  margin,
  sd = NULL,
  p = 0.5,
  N = Inf,
  level = 0.95,
  deff = 1,
  response = 1,
  target = c("mean", "proportion", "total"),
  strata = NULL,
  allocation = c("proportional", "equal", "neyman"),
  min_per_stratum = 2L
)
```

## Arguments

- margin:

  Half-width of the confidence interval you want, the half-width the
  plan expects (see "Where the formula comes from"). For a proportion it
  is a proportion, and a `margin` of 1 or more is refused. Several
  margins give several plans (see "Several plans at once").

- sd:

  Standard deviation of the measure, for `target = "mean"` or `"total"`,
  as [`stats::sd()`](https://rdrr.io/r/stats/sd.html) computes it. Not
  used with `strata`.

- p:

  Expected proportion, for `target = "proportion"`. Not used with
  `strata`.

- N:

  Population size, at least 1. `Inf` for an effectively unbounded frame.
  An estimated `N` need not be a whole number; `n` is at most
  `floor(N)`. Not used with `strata`, whose sizes add up to it.

- level:

  Confidence level.

- deff:

  Design effect to inflate by: a positive number, or the result of
  [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
  for a mean or a proportion, or of
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  for a total, on a pilot sample. `1` assumes simple random sampling. A
  number should come from
  [`deff()`](https://elkronos.github.io/dRawn/reference/deff.md) on the
  estimator you will report, rounded up, or from
  [`deff_cluster()`](https://elkronos.github.io/dRawn/reference/deff_cluster.md)
  (see "The corrections"). With `strata` it inflates the variance of the
  stratified sample, as for clusters sampled within strata.

- response:

  Expected response rate, between 0 and 1.

- target:

  `"mean"`, `"proportion"` or `"total"`. For `"total"` the margin is on
  the population total and `N` must be finite.

- strata:

  `NULL`, or a data frame with a row per stratum, to plan a stratified
  sample (see "Planning a stratified sample"). Its columns are those
  that name the strata; `N`, the rows of the frame in each stratum;
  `sd`, the standard deviation of the measure within each stratum, or
  for a proportion `p`, the expected proportion within it; and, with
  `allocation = "neyman"` only, an optional `cost`, the cost of one row
  in each stratum.

- allocation:

  How `n` is split over `strata`: `"proportional"`, `"equal"` or
  `"neyman"`, as in
  [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md).
  Neyman shares are proportional to `N_h * S_h`, or
  `N_h * S_h / sqrt(cost_h)` with a `cost` column.

- min_per_stratum:

  The fewest rows a stratum is allocated, as in
  [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
  1 or more. A stratum smaller than this is taken whole. The default of
  2 is the fewest that give each stratum a variance, which
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  and
  [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
  need from every stratum not taken whole. 0, the default of
  [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
  is refused, since a split without a floor can leave a small stratum
  without rows, and the estimate then leaves it out.

## Value

A plan, a list of class `drawn_plan` with a
[`print()`](https://rdrr.io/r/base/print.html) method, or for several
plans a data frame (see "Several plans at once"). A plan holds:

- `n`:

  Draw this many. Never more than `N`, and without `strata` at least two
  when the frame has two rows or more.

- `n_effective`:

  What you expect to analyse once `response` has taken its share.

- `capped`:

  `TRUE` when only a census meets the margin, so `n` is every row of the
  frame. A census meets any margin when everyone responds.

- `short`:

  `TRUE` when a sample smaller than the frame would meet the margin if
  everyone responded, but not at this `response` rate, so `n` is the
  whole frame and the margin achieved will be wider than the one asked
  for. It is never `TRUE` with `capped`.

- `raised`:

  `TRUE` when the margin alone needs one row, and `n_effective` was
  raised to two, the fewest that give a standard error and an interval.

- `margin`, `level`, `N`, `deff`, `response`, `target`:

  The inputs, returned so the assumptions travel with the number. `deff`
  is the number used, also when an estimate was passed.

- `p`:

  The expected proportion for a proportion, and `NA` for a mean or a
  total.

- `n_needed`:

  Rows that would reach the margin if everyone responded, two at least
  as for `n_effective`. Equal to `n_effective` unless the whole frame is
  drawn at a `response` below 1.

- `spread`:

  The standard deviation used, `sd`, or `sqrt(p * (1 - p))` for a
  proportion. `NA` with `strata`.

With `strata` there are four more, and `raised` is always `FALSE`, since
`min_per_stratum` sets the floor:

- `allocation`:

  The rows to draw from each stratum, a data frame of the strata columns
  and `n`, in the order of `strata`. Pass it to
  [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md)
  as `allocation` to draw that split.

- `strata`:

  The strata as given, checked, with `n` added.

- `rule`, `min_per_stratum`:

  The `allocation` rule and the floor used.

A plan can be passed as `n` to any design that takes `n`, such as
`design_simple(n = plan)`, which draws `plan$n` rows.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) turns it
into a data frame of one row, with a column for each field above, except
that a stratified plan has the number of `strata` in place of the two
data frames. Several plans at once give a data frame of those rows, one
per plan.

## What you need to supply

A margin of error, and some idea of how much the thing you are measuring
varies:

- For a **mean or total**, an `sd` — from a pilot, last year's data, or
  a range divided by four as a rough stand-in.

- For a **proportion**, a `p`. Leave it at `0.5`, the most pessimistic
  value, unless you have a better guess; that is the honest default
  because it maximises the required size. The margin is a proportion
  too, `0.03` for 3 percentage points.

`sd` given for a proportion, or `p` for a mean or a total, is not used,
and `plan_size()` warns that it is ignored.

For a stratified sample, give the size and the spread of each stratum in
`strata` instead (see "Planning a stratified sample").

## The corrections

`N` applies a finite population correction: sampling 400 from 500 is
very different from 400 from 500,000, and past a certain point a bigger
frame stops mattering. `response` inflates for non-response: at `0.6`
you draw enough to end up with what you need.

`deff` inflates for the design, since a clustered design of 400 may
carry the information of 100. Before any pilot,
[`deff_cluster()`](https://elkronos.github.io/dRawn/reference/deff_cluster.md)
gives the design effect of taking `m` rows from each cluster from their
intraclass correlation. Measure it with
[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md) on a
comparable past sample, using the estimator you will report:
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md) for
a mean or a proportion, and
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
for a total. The two design effects differ for the same sample, and can
differ many times over. In one cluster sample with clusters of 2 to 40
rows, the total's design effect was 28 and the mean's 0.66, and a total
planned with the mean's figure came out a small fraction of the size it
needed. A total reported as `N` times a mean takes the mean's design
effect, and is planned as a mean with `margin / N`. You can pass the
result of
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md) or
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
itself as `deff`, and `plan_size()` refuses one of the wrong kind for
`target`.

A design effect from one pilot is an estimate, and a noisy one, which
comes out below the truth more often than above it. Round it up.
[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md) gives
`NA` for a pilot without a variance estimate, such as one taken with
`variance = "none"`, a census, or one with a single cluster, and for one
in which `y` took a single value. It is `NA` too for a mean from
`ht_mean(na_rm = TRUE)` with `y` missing in some sampled rows, since the
rows with `y` are a domain, whose design effect is not estimated. Such a
pilot cannot supply a design effect, and `plan_size()` refuses an `NA`.
It refuses an estimate from `survey` too, such as
`survey::svymean(deff = TRUE)`, whose value is the estimate itself; pass
[`survey::deff()`](https://rdrr.io/pkg/survey/man/surveysummary.html) of
it instead.

## Planning a stratified sample

With `strata`, `plan_size()` plans a stratified simple random sample, as
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md)
draws one. For each total `n` it splits `n` over the strata as that
design would, by the same `allocation` rule, the floor `min_per_stratum`
and largest-remainder rounding, with no stratum given more rows than it
has, and works out the variance of the stratified mean,
`deff * sum(W_h^2 * (1 / n_h - 1 / N_h) * S_h^2)` with `W_h = N_h / N`
(Cochran 1977, chapter 5). For a proportion `S_h^2` is
`N_h * p_h * (1 - p_h) / (N_h - 1)`. The plan is the smallest `n` whose
variance is at most `(margin / z)^2`, with `margin / N` for a total.
`min_per_stratum` must be 1 or more, since an estimate leaves a stratum
without rows out. Without a floor, proportional allocation gives a
stratum of one row in a frame of 200,000 its first row at `n` of about
100,000, and the plan would be that large. The rounded split can take a
row from one stratum as `n` grows, so the variance does not always fall
with `n`, and the search checks each `n` from the first that could meet
the margin.

`response` divides `n` as without strata, and the rows to draw are split
by the same rule. The split is returned as `allocation`, which
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md)
takes as the count for each stratum. Passing the plan as `n` instead
lets the design split `n` again by its own rule and `min_per_stratum`,
whose default there is 0, which can leave a small stratum without rows.

## Several plans at once

`margin`, `sd`, `p`, `deff` and `response` each take one value, or
several for several plans. Those with one value are used for every plan,
and those with several must have the same number, one per plan, so
`plan_size(margin = c(3, 5, 10), sd = 40)` plans three margins at the
same spread. The result is then a data frame with a row per plan, as
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives for
each, without the split over strata or the printed notes.

## Where the formula comes from

`n0 = (z * spread / margin)^2 * deff` is the size for an unbounded
frame, with `z` the normal quantile for `level` and, for a total,
`margin / N` as the margin on the mean. For a mean or a total the finite
population correction gives `n0 / (1 + n0 / N)`, the smallest `n` with
`(1 / n - 1 / N) * deff * sd^2 <= (margin / z)^2`. For a proportion,
whose variance `p * (1 - p)` has the divisor `N` rather than `N - 1`, it
gives `n0 / (1 + (n0 - 1) / N)`. The result is rounded up, raised to two
when it is one, and divided by `response` (Cochran 1977, chapter 4;
Valliant, Dever and Kreuter 2018, chapter 3).

`margin` is the half-width the plan expects, the usual convention for
sizing a sample: the estimate from a sample of size `n` lands within
`margin` of the truth with probability about `level`. The interval a
sample reports is built from that sample's own estimate of the variance,
which comes out above the planned variance about as often as below it,
so the reported half-width is wider than `margin` about half the time.

That interval also uses t on the design's degrees of freedom rather than
`z`, which matters when the primary sampling units are few. With 3, 4
and 6 clusters in one stratum, which give 2, 3 and 5 degrees of freedom,
the 95% t quantile is 2.2, 1.6 and 1.3 times `z`. `plan_size()` counts
rows and knows nothing of clusters, so when you turn `n` into a few
clusters, plan more of them.

For a proportion,
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
reports an interval on the logit scale, which is asymmetric. When `p` is
small its upper arm is longer than its lower one, and exceeds `margin`
in most samples even when the estimate is as precise as planned.

## References

Cochran, W. G. (1977). *Sampling Techniques*, 3rd ed. Wiley.

Kish, L. (1965). *Survey Sampling*. Wiley.

Valliant, R., Dever, J. A. and Kreuter, F. (2018). *Practical Tools for
Designing and Weighting Survey Samples*, 2nd ed. Springer.
[doi:10.1007/978-3-319-93632-1](https://doi.org/10.1007/978-3-319-93632-1)

## See also

[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md) to
measure the design effect of a past sample,
[`deff_cluster()`](https://elkronos.github.io/dRawn/reference/deff_cluster.md)
for one from an intraclass correlation,
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) to take
the sample.

## Examples

``` r
# A proportion, no prior guess, 20,000 in the frame
plan_size(margin = 0.03, N = 20000, target = "proportion")
#> Sample size for a proportion
#>   draw           1,014
#>   margin         +/- 0.03 at 95% confidence
#>   assuming       p 0.5, N 20,000

# A mean, when a pilot put the spread near 40
plan_size(margin = 5, sd = 40, N = 20000)
#> Sample size for a mean
#>   draw           243
#>   margin         +/- 5 at 95% confidence
#>   assuming       sd 40, N 20,000

# The same, in a clustered design with 70% response
plan_size(margin = 5, sd = 40, N = 20000, deff = 2.5, response = 0.7)
#> Sample size for a mean
#>   draw           853
#>   to analyse     597  (after 70% response)
#>   margin         +/- 5 at 95% confidence
#>   assuming       sd 40, deff 2.5, N 20,000

# A total, to within 100,000 across a 20,000-row frame
plan_size(margin = 1e5, sd = 40, N = 20000, target = "total")
#> Sample size for a total
#>   draw           243
#>   margin         +/- 100,000 at 95% confidence
#>   assuming       sd 40, N 20,000

# Clusters of ten rows, whose rows correlate at 0.05, before any pilot
plan_size(margin = 5, sd = 40, N = 20000,
          deff = deff_cluster(m = 10, icc = 0.05))
#> Sample size for a mean
#>   draw           351
#>   margin         +/- 5 at 95% confidence
#>   assuming       sd 40, deff 1.45, N 20,000

# Three margins at once, one plan per row
plan_size(margin = c(3, 5, 10), sd = 40, N = 20000)
#>     n n_effective margin level     N deff response target spread  p capped
#> 1 661         661      3  0.95 20000    1        1   mean     40 NA  FALSE
#> 2 243         243      5  0.95 20000    1        1   mean     40 NA  FALSE
#> 3  62          62     10  0.95 20000    1        1   mean     40 NA  FALSE
#>   short raised n_needed
#> 1 FALSE  FALSE      661
#> 2 FALSE  FALSE      243
#> 3 FALSE  FALSE       62

# A plan is a valid `n`
design_simple(n = plan_size(margin = 5, sd = 40, N = 20000))
#> <sampling design: simple>
#>   n        243
#>   replace  FALSE
#>   reserve  0

# Three strata, with the spread of each from an earlier survey
bands <- data.frame(band = c("small", "medium", "large"),
                    N = c(1200, 600, 200), sd = c(10, 25, 80))
plan <- plan_size(margin = 2, strata = bands, allocation = "neyman")
plan
#> Sample size for a mean, over 3 strata
#>   draw           312
#>   margin         +/- 2 at 95% confidence
#>   assuming       Neyman allocation, at least 2 rows per stratum, N 2,000
#> 
#>        band    N sd   n
#>       small 1200 10  87
#>      medium  600 25 109
#>       large  200 80 116
frame <- data.frame(id = 1:2000,
                    band = rep(bands$band, times = bands$N))
s <- draw(frame, design_stratified("band", allocation = plan$allocation),
          seed = 1)
table(s$band)
#> 
#>  large medium  small 
#>    116    109     87 

# The design effect of a cluster pilot, passed as the estimate itself
pop <- data.frame(cl = rep(paste0("c", 1:50), each = 8))
pop$y <- rep(seq(10, 59), each = 8) + rep(c(-3, 3), 200)
pilot <- draw(pop, design_cluster("cl", n_clusters = 6), seed = 1,
              weights = TRUE)
plan_size(margin = 3, sd = sd(pilot$y), N = nrow(pop),
          deff = ht_mean(pilot, "y"))
#> Sample size for a mean
#>   draw           294
#>   margin         +/- 3 at 95% confidence
#>   assuming       sd 16.86, deff 9.1, N 400
```
