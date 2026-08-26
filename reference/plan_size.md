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
  target = c("mean", "proportion", "total")
)
```

## Arguments

- margin:

  Half-width of the confidence interval you want.

- sd:

  Standard deviation of the measure, for `target = "mean"` or `"total"`.

- p:

  Expected proportion, for `target = "proportion"`.

- N:

  Population size. `Inf` for an effectively unbounded frame.

- level:

  Confidence level.

- deff:

  Design effect to inflate by. `1` assumes simple random sampling.

- response:

  Expected response rate, between 0 and 1.

- target:

  `"mean"`, `"proportion"` or `"total"`. For `"total"` the margin is on
  the population total and `N` must be finite.

## Value

A list with a [`print()`](https://rdrr.io/r/base/print.html) method,
holding:

- `n`:

  Draw this many. Never more than `N`.

- `n_effective`:

  What you expect to analyse once `response` has taken its share.

- `capped`:

  `TRUE` when the frame is not large enough to reach the margin at all,
  however you sample it.

- `short`:

  `TRUE` when the frame could reach the margin but not at this
  `response` rate, so `n` is the whole frame and the margin achieved
  will be wider than the one asked for.

- `margin`, `level`, `N`, `deff`, `response`, `target`:

  The inputs, returned so the assumptions travel with the number.

- `n_needed`:

  Rows that would reach the margin if everyone responded. Equal to
  `n_effective` unless `short` is `TRUE`.

- `spread`:

  The standard deviation used — `sd`, or `sqrt(p * (1 - p))` for a
  proportion.

## What you need to supply

A margin of error, and some idea of how much the thing you are measuring
varies:

- For a **mean or total**, an `sd` — from a pilot, last year's data, or
  a range divided by four as a rough stand-in.

- For a **proportion**, a `p`. Leave it at `0.5`, the most pessimistic
  value, unless you have a better guess; that is the honest default
  because it maximises the required size.

## The corrections

`N` applies a finite population correction: sampling 400 from 500 is
very different from 400 from 500,000, and past a certain point a bigger
frame stops mattering. `deff` inflates for the design — pass the
[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md) from a
comparable past sample, since a clustered design of 400 may carry the
information of 100. `response` inflates for non-response: at `0.6` you
draw enough to end up with what you need.

## See also

[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md) to
measure the design effect of a past sample,
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) to take
the sample.

## Examples

``` r
# A proportion, no prior guess, 20,000 in the frame
plan_size(margin = 0.03, N = 20000, target = "proportion")
#> Sample size for a proportion
#>   draw           1,014
#>   margin         +/- 0.03 at 95% confidence
#>   assuming       sd 0.5, N 20,000

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
#>   margin         +/- 1e+05 at 95% confidence
#>   assuming       sd 40, N 20,000
```
