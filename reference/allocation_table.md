# How a design splits its sample over the strata of a frame

The rows [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)
takes from each stratum, the inclusion probability that gives each row
there, and the reason for any stratum held at the floor, capped at its
size or left with no rows, all before anything is drawn. Under Neyman
allocation it also shows the standard deviation of each stratum the
split used, and the cost per row when the design has costs. The numbers
come from the same computation as
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) and
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md),
so they always agree with both.

## Usage

``` r
allocation_table(data, design)
```

## Arguments

- data:

  A data frame.

- design:

  [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
  [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
  or
  [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
  with `strata`, or
  [`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md)
  whose `rest` is one of these. For the last the table describes `rest`
  over the rows below the threshold, and the rows at or above it, which
  are taken with certainty, are not counted in any stratum. Other
  designs are refused.

## Value

A data frame with one row per stratum of `data`, in key order (see
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md)):
the `strata` columns, then

- `N`:

  Rows in the stratum.

- `n`:

  Rows drawn from it.

- `prob`:

  The inclusion probability of each of its rows, as
  [`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
  gives it. `NA` for
  [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
  with `size` and for
  [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md),
  whose rows in one stratum have different probabilities; there
  [`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
  gives them row by row.

- `sd`:

  Under Neyman allocation, the stratum's standard deviation that the
  split used: of the `allocation_by` column within the stratum, or as
  given. `NA` under other allocations.

- `cost`:

  The cost of one row, when the design has `cost`, and `NA` otherwise.

- `status`:

  Why the stratum got `n`: `"rule"` for its share under the allocation
  rule, rounded; `"count"` for a count given in `allocation`; `"floor"`
  for a stratum raised to `min_per_stratum` from a smaller share;
  `"capped"` for one the rule asked more rows of than it holds, so it is
  taken whole and the rows it could not give go to the other strata;
  `"take-all"` for one taken whole for another reason, such as a size
  below `min_per_stratum` or a count equal to its size; and `"empty"`
  for one allocated no rows, whose rows can never be drawn. With
  `replace = TRUE` no stratum is capped or taken whole.

## See also

[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md),
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)

## Examples

``` r
df <- data.frame(id = 1:300,
                 site = rep(c("a", "b", "c", "d"), c(150, 100, 45, 5)),
                 v = c(1:150 %% 7, 1:100 %% 19, 1:45 %% 3, 1:5 * 40))
d <- design_stratified("site", n = 30, allocation = "neyman",
                       allocation_by = "v", min_per_stratum = 2)
allocation_table(df, d)
#>   site   N  n       prob         sd cost status
#> 1    a 150  8 0.05333333  1.9948592   NA   rule
#> 2    b 100 15 0.15000000  5.5331873   NA   rule
#> 3    c  45  2 0.04444444  0.8257228   NA  floor
#> 4    d   5  5 1.00000000 63.2455532   NA capped

# The same counts as a draw, and the same probabilities as inclusion_prob()
table(draw(df, d, seed = 1)$site)
#> 
#>  a  b  c  d 
#>  8 15  2  5 
tapply(inclusion_prob(df, d), df$site, unique)
#>          a          b          c          d 
#> 0.05333333 0.15000000 0.04444444 1.00000000 
```
