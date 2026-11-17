# Systematic sampling

Walks the rows at a fixed step. Give `interval` to take every
`interval`-th row, beginning at `start`, which is drawn at random unless
given; the number of rows drawn then follows from the interval and the
size of the data. Give `n` instead to draw exactly `n` rows from a
random start: the step is the number of rows walked divided by `n`, and
need not be a whole number. With `n_starts`, several random starts are
walked at a longer step, which gives the sample a design-unbiased
variance. Given neither `interval` nor `n`, the design walks each
stratum of
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md)
as its `within` argument.

## Usage

``` r
design_systematic(
  interval = NULL,
  start = NULL,
  order_by = NULL,
  na_rm = FALSE,
  n = NULL,
  n_starts = 1L
)
```

## Arguments

- interval:

  Sampling interval, a whole number of 1 or more: every `interval`-th
  row of the walk is taken. Give `interval` or `n`, not both. For a step
  that is not a whole number, give `n`.

- start:

  Position in the walk of the first row taken, in `1:interval`. The walk
  is the frame's rows in their order, or sorted by `order_by` when
  given, so with `order_by` a `start` of 3 takes the row with the third
  smallest value, not the frame's third row. Drawn at random from that
  range when `NULL`. A fixed start makes the design deterministic: each
  row is selected with probability 0 or 1, so
  [`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md),
  [`sampling_weight()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md),
  [`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md),
  `draw(weights = TRUE)`,
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  and
  [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
  refuse it, while `inclusion_prob(simulate = TRUE)` reports those 0s
  and 1s. Not used with `n` or `n_starts`, whose starts are always
  random.

- order_by:

  Optional column to sort by before walking the data. It changes which
  rows can appear together, so
  [`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
  and
  [`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)
  compute against the sorted order too. Text is sorted byte by byte, as
  in the C locale, so the walk is the same in every session; ties keep
  frame order, so where the column has ties a seed is sure to draw the
  same rows only from a frame with its rows in the same order. Without
  ties the walk does not depend on frame order. One column: to walk by
  several, add a column of ranks built from them and give its name, as
  in the examples.

- na_rm:

  When `order_by` is given, drop rows whose sort key is `NA` instead of
  raising an error.

- n:

  Number of rows to draw, in place of `interval`. With `N` rows walked,
  each row covers a stretch of length `n / N` on a line from 0 to `n`, a
  start `u` is drawn uniformly from (0, 1), and the rows under `u`,
  `u + 1`, ..., `u + n - 1` are taken. The step is `N / n` rows, every
  row has probability `n / N`, and the sample always holds `n` rows.
  When `n` divides `N` the probabilities are those of
  `interval = N / n`. A plan from
  [`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md)
  gives its `n`.

- n_starts:

  Number of random starts, 1 by default. With `n_starts = m` and
  `interval = k`, `m` different starts are drawn from `1:(m * k)` and
  each is walked at a step of `m * k` rows. Every row still has
  probability `1 / k`, and the expected sample size is that of one start
  at `interval`. Each start gives a replicate, a systematic sample of
  its own, and two or more of them give a design-unbiased variance (see
  "Variance"). Needs `interval`, and a random `start`.

## Value

A design object, for use with
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

## Variance

A systematic sample has a single random start. At a step of 2 rows or
more, which is every `interval` above 1 and every `n` up to half the
rows walked, neighbouring rows are never drawn together, so no
design-unbiased variance estimator exists.
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
uses the successive-difference approximation (Wolter 2007), which
compares neighbouring sampled rows in the order the design walked. With
`n` above half the rows walked, every pair of rows can be drawn together
(see
[`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)),
and
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
uses the same approximation. Sorting on a variable related to what you
measure (`order_by`) is what makes systematic sampling efficient, and
this estimator is the one that can see it; see
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
for its limits. When `interval` does not divide the number of rows
walked, the sample size is random, and a term for that is added to the
variance. A design given `n` has a fixed size and needs no such term.

With `n_starts = m` of 2 or more, the `m` starts are a simple random
sample of the `m * k` possible ones, and each replicate `r` gives a
total of its own, `t_r`, which is `m` times the sum of `y / pi` over its
rows. The sample's total is their mean, and
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
uses the random group variance (Wolter 2007),
`(1 - 1 / k) * var(t_r) / m`, on `m - 1` degrees of freedom. It is
design-unbiased whether or not `m * k` divides the rows walked, a
replicate that drew no rows counting with a total of 0.
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
gives the same standard error with each replicate as a primary sampling
unit, and refuses a sample with a replicate that drew no rows, which
`survey` would leave out of its count of replicates. The price is
precision: `m` starts spread the sample over the walk less evenly than
one start at `interval` does, and few replicates give few degrees of
freedom and a wide interval.

Within strata, as `design_stratified(within = design_systematic())`,
each stratum has a start of its own and the variance is the sum of the
strata's successive-difference approximations; see
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md).

## References

Madow, W. G. and Madow, L. H. (1944). On the theory of systematic
sampling, I. *Annals of Mathematical Statistics*, 15, 1–24.

Wolter, K. M. (2007). *Introduction to Variance Estimation*, 2nd ed.
Springer.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md),
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md)

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md),
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md),
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md),
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
[`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md),
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)

## Examples

``` r
df <- data.frame(id = 1:100, value = (1:100) / 10,
                 site = rep(c("b", "a"), each = 50))
draw(df, design_systematic(interval = 10, start = 3))
#>    id value site
#> 1   3   0.3    b
#> 2  13   1.3    b
#> 3  23   2.3    b
#> 4  33   3.3    b
#> 5  43   4.3    b
#> 6  53   5.3    a
#> 7  63   6.3    a
#> 8  73   7.3    a
#> 9  83   8.3    a
#> 10 93   9.3    a
draw(df, design_systematic(interval = 25, order_by = "value"), seed = 1)
#>    id value site
#> 1  25   2.5    b
#> 2  50   5.0    b
#> 3  75   7.5    a
#> 4 100  10.0    a

# 30 rows from 100, a step of 3.33
draw(df, design_systematic(n = 30), seed = 1)$id
#>  [1]  1  5  8 11 15 18 21 25 28 31 35 38 41 45 48 51 55 58 61 65 68 71 75 78 81
#> [26] 85 88 91 95 98

# Four starts from 1:40, each walked at a step of 40: ten rows on
# average, as with one start, and a variance from the spread of the four
# replicates
s <- draw(df, design_systematic(interval = 10, n_starts = 4), seed = 1,
          weights = TRUE)
ht_total(s, "value")
#> Horvitz-Thompson total  (systematic design, n = 10)
#>   estimate 481
#>   se       18.99737  (analytic)
#>   95% CI  420.5419 to 541.4581  (t, 3 df)
#>   deff     0.0424  (better than simple random sampling)

# Walking by two columns: rank the rows by both, then walk the ranks
df$walk <- order(order(df$site, -df$value, method = "radix"))
draw(df, design_systematic(interval = 20, order_by = "walk"), seed = 1)
#>   id value site walk
#> 1 17   1.7    b   84
#> 2 37   3.7    b   64
#> 3 57   5.7    a   44
#> 4 77   7.7    a   24
#> 5 97   9.7    a    4
```
