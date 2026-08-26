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
  na_rm = FALSE
)
```

## Arguments

- weights:

  Column of positive, finite weights.

- n:

  Number of rows to draw. Under `method = "poisson"` this is the
  expected number, not a guarantee.

- replace:

  Sample with replacement? Only available for `method = "successive"`.

- method:

  One of `"successive"`, `"systematic"` or `"poisson"`. See "Choosing a
  method".

- na_rm:

  Drop rows whose weight is `NA` instead of raising an error.

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
  at `n = 5` give inclusion probabilities up to 35% away from
  proportional. Treating them as if they were proportional biases a
  Horvitz-Thompson total by around 1% at `n = 15` out of 30 in
  simulation, small in size but unmistakable in sign. There is no closed
  form for `pi`, so
  [`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
  refuses to give one. Fine when you want a weighted selection; wrong as
  the basis for an estimate.

- `"systematic"`:

  Systematic probability-proportional-to-size. Walks the cumulative
  weights with a fixed step from a random start, giving `pi_i = n * p_i`
  exactly. Rows heavy enough that `n * p_i > 1` are taken with certainty
  and the rest rescaled, repeatedly, until every probability is valid.
  Fixed sample size. Some pairs of rows can never appear together, so
  joint inclusion probabilities are zero for them and variance
  estimation needs care.

- `"poisson"`:

  Each row is included independently with probability `pi_i = n * p_i`,
  capped at 1. Inclusion probabilities are exactly proportional to size
  and every pair can co-occur, at the cost of a **random sample size**
  averaging `n`.

`replace = TRUE` is with-replacement PPS and applies only to
`method = "successive"`.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md),
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md),
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md),
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
```
