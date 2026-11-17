# Weighted sampling

Draws rows with probability governed by a weights column. The weights
are used as supplied: `prob=` normalises by their sum, so rescaling them
changes nothing and the design does not offer it.

## Usage

``` r
design_weighted(
  weights,
  n,
  replace = FALSE,
  method = c("successive", "systematic", "poisson"),
  na_rm = FALSE,
  strata = NULL,
  allocation = c("proportional", "equal")
)
```

## Arguments

- weights:

  Column of positive, finite, numeric weights. A weight of 0 is refused
  under every method, not given a probability of 0; drop such rows from
  the frame first, since a design that selects in proportion to size
  could never select them.

- n:

  Number of rows to draw, the total across `strata` when given, or a
  plan from
  [`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md).
  Under `method = "poisson"` this is the expected number, not a
  guarantee. With counts in `allocation` it can be left out and is their
  sum; when given, it must equal that sum.

- replace:

  Sample with replacement? Only available for `method = "successive"`.

- method:

  One of `"successive"`, `"systematic"` or `"poisson"`. See "Choosing a
  method".

- na_rm:

  Drop rows whose weight or stratum key is `NA` instead of raising an
  error.

- strata:

  Optional column name(s) defining strata, to sample in proportion to
  size within each one. Several columns are cross-classified. Needs
  `method = "systematic"` or `"poisson"`. See "Strata".

- allocation:

  How `n` is split across `strata`: `"proportional"` to the number of
  rows in each stratum, or `"equal"`, by the largest-remainder method of
  [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
  with a stratum that would get more rows than it holds taken whole; or
  the count for each stratum, as a named vector such as
  `c(north = 10, south = 5)` when `strata` is one column, or a data
  frame of the `strata` columns and `n`. Counts follow the rules of
  [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md):
  every stratum of the frame needs one, a count for a stratum the frame
  does not have is refused, and no count may exceed the rows in its
  stratum. Only with `strata`.

## Value

A design object, for use with
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

## Choosing a method

The three methods differ in what the weights actually control, which
decides whether the sample can be used for estimation.

- `"successive"`:

  The default, and what
  [`base::sample()`](https://rdrr.io/r/base/sample.html) does. The
  weights govern each sequential draw, not the probability that a row
  ends up in the sample, and the realised inclusion probabilities are
  **not** proportional to the weights. How far off depends on the
  weights: with a single dominant unit the two nearly coincide (about 1%
  apart), but with a moderate spread the gap is large — weights `1:10`
  at `n = 5` give inclusion probabilities up to about 36% away from
  proportional, by exact enumeration of the ordered draws. Treating them
  as if they were proportional biases a Horvitz-Thompson total by around
  1% at `n = 15` out of 30 in simulation, small in size but unmistakable
  in sign. There is no closed form for `pi`, so
  [`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
  refuses to give one. Fine when you want a weighted selection; wrong as
  the basis for an estimate.

- `"systematic"`:

  Randomised systematic probability-proportional-to- size (Madow 1949;
  Hartley and Rao 1962). Shuffles the rows, then walks the cumulative
  weights with a fixed step from a random start, giving `pi_i = n * p_i`
  exactly. Rows heavy enough that `n * p_i > 1` are taken with certainty
  and the rest rescaled, repeatedly, until every probability is valid.
  Fixed sample size. The shuffle means every pair of rows can appear
  together, but their joint probabilities have no closed form, so
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
  uses Deville's approximation for the variance, which needs only the
  first-order probabilities.

- `"poisson"`:

  Each row is included independently with probability `pi_i = n * p_i`.
  Rows heavy enough that `n * p_i > 1` are taken with certainty and the
  rest rescaled, repeatedly, as for `"systematic"`, so the probabilities
  still sum to `n`. Inclusion probabilities are exactly proportional to
  size among the rows below certainty, and every pair can co-occur, at
  the cost of a **random sample size** averaging `n`.

`replace = TRUE` is with-replacement PPS and applies only to
`method = "successive"`.

## Strata

With `strata`, `n` is split over the strata by `allocation`, as
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md)
splits it over the rows in each stratum, and each stratum is then
sampled in proportion to size on its own: a row of stratum `h` gets
`n_h * w_i / sum(w)`, the sum running over the stratum, with rows past
certainty taken with certainty within it. Under `"systematic"` each
stratum is shuffled and walked from a start of its own, in key order;
under `"poisson"` every row is included independently with the
probability its stratum gives it, so the size of each stratum's sample
is random. `"successive"` has no closed-form inclusion probability and
is refused with `strata`. To split `n` in proportion to each stratum's
total size rather than its rows, work out the counts and give them as
`allocation`.

The strata are sampled independently.
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
and [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
add up Deville's approximation over the strata for `"systematic"`, and
the Poisson variance needs no change. The degrees of freedom are the
sampled rows below certainty less the strata that hold them. A stratum
with one sampled row below certainty has nothing to compare it with, and
the variance is refused with the stratum named; draw at least two in
each.

In simulations of 2,000 draws each from frames of 600 rows in four
strata, drawing 24 or 60 rows by `"systematic"`, the variance averaged
0.94 to 1.04 times the true variance of the total, and 95% intervals
covered 89% to 96% of the time. The low coverage came where `y / pi` was
skewed, so that samples with a low total also had a low estimated
variance; the same frames without strata covered 92% to 94%.

[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
hands the strata to `survey`, with Brewer's approximation within each,
[`allocation_table()`](https://elkronos.github.io/dRawn/reference/allocation_table.md)
shows the split before anything is drawn, and
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
counts the rows drawn from each stratum.

Drawing each stratum with a design of its own and binding the samples
with [`rbind()`](https://rdrr.io/r/base/cbind.html) does not give a
sample to estimate from. The bound sample keeps the first draw's frame
and design, so
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
refuses it; use `strata` instead.

## References

Madow, W. G. (1949). On the theory of systematic sampling, II. *Annals
of Mathematical Statistics*, 20, 333–354.

Hartley, H. O. and Rao, J. N. K. (1962). Sampling with unequal
probabilities and without replacement. *Annals of Mathematical
Statistics*, 33, 350–374.

Brewer, K. R. W. and Hanif, M. (1983). *Sampling with Unequal
Probabilities*. Springer.

Tillé, Y. (2006). *Sampling Algorithms*. Springer.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md),
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md),
[`certainty_set()`](https://elkronos.github.io/dRawn/reference/certainty_set.md)
for the rows taken with certainty,
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
to select in proportion to size *and* spread the sample across a map.
The `sampling` package implements further fixed-size unequal-probability
algorithms (Brewer, Sampford, Tillé, maximum entropy) as vector
functions.

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
[`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md),
[`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md)

## Examples

``` r
df <- data.frame(id = 1:20, w = 1:20)
draw(df, design_weighted("w", n = 5), seed = 1)
#>   id  w
#> 1 18 18
#> 2 16 16
#> 3 12 12
#> 4  6  6
#> 5 19 19

# Inclusion probabilities proportional to size, and checkable
d <- design_weighted("w", n = 5, method = "systematic")
round(inclusion_prob(df, d), 3)
#>  [1] 0.024 0.048 0.071 0.095 0.119 0.143 0.167 0.190 0.214 0.238 0.262 0.286
#> [13] 0.310 0.333 0.357 0.381 0.405 0.429 0.452 0.476

# Which is what makes an unbiased total possible
s <- draw(df, d, seed = 1, weights = TRUE)
sum(s$w * s$.weight)   # estimates sum(df$w) = 210
#> [1] 210

# In proportion to size within each of two regions
df$region <- rep(c("east", "west"), each = 10)
by_region <- design_weighted("w", n = 8, method = "systematic",
                             strata = "region")
allocation_table(df, by_region)
#>   region  N n prob sd cost status
#> 1   east 10 4   NA NA   NA   rule
#> 2   west 10 4   NA NA   NA   rule
round(inclusion_prob(df, by_region), 3)
#>  [1] 0.073 0.145 0.218 0.291 0.364 0.436 0.509 0.582 0.655 0.727 0.284 0.310
#> [13] 0.335 0.361 0.387 0.413 0.439 0.465 0.490 0.516
```
