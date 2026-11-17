# Rows a size-proportional design takes with certainty

A design whose inclusion probabilities are proportional to size gives a
row whose share of `n` reaches 1 a probability of exactly 1, and spreads
what is left of `n` over the other rows, repeating until no share passes

1.  `certainty_set()` says which rows of `data` those are, and where the
    cut-off falls, before any sample is drawn.

## Usage

``` r
certainty_set(data, design)
```

## Arguments

- data:

  A data frame.

- design:

  [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
  with `method = "systematic"` or `"poisson"`, or
  [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
  with a `size` column. Other designs are refused, since none of them
  has inclusion probabilities proportional to size.

## Value

An integer vector of the rows of `data` whose inclusion probability is
1, in frame order, the same rows as
`which(inclusion_prob(data, design) == 1)`. Its attribute `threshold` is
the smallest size among them, or `NA` when there are none. Of the rows
the design can select, those whose size is at least `threshold` are the
set. With `strata`,
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
and
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
put the cut-off in each stratum at its own size, so there `threshold` is
`NA`. For
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
without `strata`,
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md)
with that `threshold` on the weights column therefore takes the same
rows with certainty, and a `rest` of the same method with `n` reduced by
the number of rows in the set gives the other rows the same
probabilities too.

## See also

[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md),
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md)

## Examples

``` r
df <- data.frame(id = 1:10, w = c(1, 2, 3, 4, 5, 6, 7, 8, 60, 90))
d <- design_weighted("w", n = 4, method = "systematic")
cs <- certainty_set(df, d)
df[cs, ]
#>    id  w
#> 9   9 60
#> 10 10 90
attr(cs, "threshold")
#> [1] 60

# The same probabilities, with the cut-off made explicit
same <- design_certainty("w", attr(cs, "threshold"),
                         rest = design_weighted("w", n = 4 - length(cs),
                                                method = "systematic"))
all.equal(inclusion_prob(df, same), inclusion_prob(df, d))
#> [1] TRUE
```
