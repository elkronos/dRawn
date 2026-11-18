# Stratified sampling

Draws a sample from each stratum. Allocation uses the largest-remainder
method, so the returned row count matches `n` exactly rather than
drifting with per-stratum rounding.

## Usage

``` r
design_stratified(
  strata,
  n,
  allocation = c("proportional", "equal", "neyman"),
  allocation_by = NULL,
  min_per_stratum = 0L,
  replace = FALSE,
  na_rm = FALSE,
  cost = NULL,
  within = NULL
)
```

## Arguments

- strata:

  One or more column names defining the strata. Several columns are
  cross-classified.

- n:

  Total rows to draw across all strata, under any allocation, or a plan
  from
  [`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md).
  With counts in `allocation` it can be left out and is their sum; when
  given, it must equal that sum.

- allocation:

  How `n` is split across strata: a rule, or the count for each stratum.
  Of the rules, `"proportional"` gives each stratum a share of `n` in
  proportion to its size; `"equal"` splits `n` evenly; `"neyman"` gives
  shares proportional to `size * sd`, with the standard deviation of
  each stratum from `allocation_by` (Neyman 1934). Neyman minimises the
  variance of a total for a fixed `n` by putting more rows where the
  values vary most, and is the right choice when you have a frame
  variable correlated with what you are measuring, or estimates of each
  stratum's standard deviation from a pilot or an earlier survey. Under
  any rule, a stratum that would be allocated more rows than it holds is
  taken whole, one that would get fewer than `min_per_stratum` gets that
  many, and the rest of `n` is re-split over the others by the same rule
  (Cochran 1977, section 5.9).

  Counts are a named vector, such as
  `c(north = 10, south = 5, east = 5)`, when `strata` is one column, its
  names matched to that column's values as text; or a data frame holding
  the `strata` columns and a column `n`, matched on the values
  themselves. Every stratum of the frame needs a count, 0 included, and
  a count for a stratum the frame does not have is refused, so a
  misspelt name cannot quietly leave a stratum out. Without replacement
  no count may exceed the rows in its stratum. The design then draws
  exactly those counts, and the inclusion probabilities, variances and
  [`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
  follow them.

- allocation_by:

  The standard deviation of each stratum, for `allocation = "neyman"`.
  Ignored otherwise. Either a column of the frame, whose standard
  deviation within each stratum is used, and which needs a finite value
  on every row allocated over; or the standard deviations themselves, so
  the frame needs no such column. Those are a named vector, such as
  `c(north = 40, south = 15)`, when `strata` is one column, its names
  matched to that column's values as text; or a data frame holding the
  `strata` columns and a column `sd`, matched on the values themselves.
  Every stratum of the frame needs one, and each must be finite and 0 or
  more. A value for a stratum the frame does not have is ignored, so the
  same values serve a part of the frame, such as the `rest` of
  [`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md)
  when every row of a stratum is above the threshold. A stratum with a
  standard deviation of 0 gets no share.

- min_per_stratum:

  Minimum rows from each stratum, under an allocation rule. With the
  default of `0`, a stratum smaller than about `N / n` can be allocated
  no rows at all; its rows then have inclusion probability 0 and a total
  estimated from the sample silently leaves them out, so
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) warns
  when that happens. `1` makes every stratum reachable, which is what an
  unbiased Horvitz-Thompson total needs, and `2` also lets every stratum
  contribute to the variance estimate. Over-sampling a small stratum
  this way does *not* bias an estimate made with the design weights,
  which correct for it; it only moves precision around. Without
  replacement a stratum smaller than the floor is taken whole, so `n`
  must cover `min_per_stratum` rows from each stratum, or all of a
  smaller one. With counts in `allocation`, give each stratum its count
  there instead.

- replace:

  Sample with replacement within each stratum?

- na_rm:

  Drop rows whose stratum key is `NA` instead of raising an error.

- cost:

  The cost of sampling one row in each stratum, for
  `allocation = "neyman"` only. In the same forms as standard deviations
  given in `allocation_by`, with the column of a data frame named
  `cost`, so one data frame of the `strata` columns, `sd` and `cost` can
  serve both. Each must be finite and above 0. The shares become
  proportional to `size * sd / sqrt(cost)`, which gives the smallest
  variance among the splits of the same total cost, `sum(cost * n_h)`
  over strata (Cochran 1977, section 5.5). Before any stratum is held at
  its floor or its size, the `n` that spends a budget `C`, net of any
  fixed cost, is
  `C * sum(N_h * S_h / sqrt(c_h)) / sum(N_h * S_h * sqrt(c_h))`, rounded
  down. When no stratum has any spread, the shares are proportional to
  `size / sqrt(cost)`, as for equal spreads, and without `cost` to
  `size`.

- within:

  How rows are drawn inside each stratum: `NULL`, the default, for a
  simple random sample, or
  [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md)
  without `interval` or `n`, to walk each stratum in the order of its
  `order_by` column, or in frame order without one. Each stratum gets
  the rows the allocation gives it, at an even step from a random start
  of its own; see "Systematic sampling within strata". Needs
  `replace = FALSE`.

## Value

A design object, for use with
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).
[`allocation_table()`](https://elkronos.github.io/dRawn/reference/allocation_table.md)
shows the rows it takes from each stratum of a frame before anything is
drawn.

## Rounding, ties and the order of strata

Each stratum first gets the whole part of its share of `n`. The rows
left over go one at a time to the strata with the largest fractional
parts. When two strata tie, the one with more rows not yet allocated
goes first, and after that the one that comes first in key order.
`min_per_stratum` and the size of each stratum bound the shares before
they are rounded, so the rounding never takes a stratum below its floor
or past its size.

Key order is the order of the levels for a factor column and the sorted
values for any other column, with text sorted byte by byte as in the C
locale (`"B"` before `"a"`). With several `strata` columns the first
varies slowest. Strata are also drawn from in this order, so a seed
selects the same rows in every locale. Renaming a stratum, or reordering
a factor's levels, can move a tied row from one stratum to another; the
inclusion probabilities always follow the allocation used.

An empty string is a stratum like any other, shown as `(blank)` in
messages.

## Systematic sampling within strata

With `within = design_systematic(order_by = "x")`, the strata are
allocated as without it, each stratum is sorted on `x`, and the `n_h`
rows allocated to a stratum of `N_h` are taken as
`design_systematic(n = n_h)` takes them from a frame of `N_h` rows: a
start drawn from (0, 1) and a step of `N_h / n_h` rows, which need not
be whole. Strata draw their starts in key order. Every row has
probability `n_h / N_h`, as under simple random sampling within strata,
so the allocation,
[`allocation_table()`](https://elkronos.github.io/dRawn/reference/allocation_table.md)
and the weights are the same; what changes is which rows can be drawn
together. Two rows of one stratum have the joint probability of the
systematic design there, which is 0 for rows at most `N_h / n_h - 1`
rows apart, and two rows of different strata the product of their
probabilities. Rows whose `x` is missing are refused, or with the
`na_rm = TRUE` of the systematic design left out of their stratum before
it is allocated.

Each stratum has one random start, so when any stratum is walked at a
step of 2 rows or more no design-unbiased variance exists, as for
[`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md)
on its own.
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
and [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
add up the strata's successive-difference approximations,
`(1 - n_h / N_h) * n_h / (2 * (n_h - 1)) * sum(diff(y / pi)^2)` over
each stratum's sampled rows in the order it was walked (Wolter 2007),
with the degrees of freedom of a stratified sample, the sampled rows
less the strata. A stratum with one sampled row from several is refused,
as without `within`. So is the jackknife: a stratum's rows all come from
one start, and deleting them one at a time does not mirror the design.
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
hands `survey` the strata, whose simple-random variance differs.

Only
[`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md)
is taken as `within`, with one start and no `interval` or `n`: the
allocation sets the size of each stratum. To walk a stratum by several
columns, give `order_by` a column of ranks, as in the examples of
[`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md).

## References

Neyman, J. (1934). On the two different aspects of the representative
method: the method of stratified sampling and the method of purposive
selection. *Journal of the Royal Statistical Society*, 97, 558–625.
[doi:10.2307/2342192](https://doi.org/10.2307/2342192)

Cochran, W. G. (1977). *Sampling Techniques*, 3rd ed. Wiley.

Wolter, K. M. (2007). *Introduction to Variance Estimation*, 2nd ed.
Springer.
[doi:10.1007/978-0-387-35099-8](https://doi.org/10.1007/978-0-387-35099-8)

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md),
[`allocation_table()`](https://elkronos.github.io/dRawn/reference/allocation_table.md)

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md),
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md),
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md),
[`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md),
[`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md),
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)

## Examples

``` r
df <- data.frame(id = 1:100, site = rep(letters[1:4], each = 25))
table(draw(df, design_stratified("site", n = 20), seed = 1)$site)
#> 
#> a b c d 
#> 5 5 5 5 

# Rare strata are covered only if you ask -- without the floor, the one
# "rare" row gets no allocation, and draw() warns that it is unreachable
skewed <- data.frame(id = 1:1000, g = c(rep("common", 999), "rare"))
draw(skewed, design_stratified("g", n = 10, min_per_stratum = 1), seed = 1)$g
#>  [1] "common" "common" "common" "common" "common" "common" "common" "common"
#>  [9] "common" "rare"  

# A count for each stratum, drawn as given; `n` is their sum
fixed <- design_stratified("site", allocation = c(a = 8, b = 4, c = 4, d = 4))
table(draw(df, fixed, seed = 1)$site)
#> 
#> a b c d 
#> 8 4 4 4 

# Neyman from standard deviations known from an earlier survey, with the
# cost of reaching a row in each stratum
pilot <- design_stratified("site", n = 20, allocation = "neyman",
                           allocation_by = c(a = 4, b = 1, c = 2, d = 8),
                           cost = c(a = 1, b = 1, c = 4, d = 4))
allocation_table(df, pilot)
#>   site  N n prob sd cost status
#> 1    a 25 8 0.32  4    1   rule
#> 2    b 25 2 0.08  1    1   rule
#> 3    c 25 2 0.08  2    4   rule
#> 4    d 25 8 0.32  8    4   rule

# Each stratum walked in order of `size`, from a start of its own
df$size <- (df$id * 37) %% 101
walked <- design_stratified("site", n = 20,
                            within = design_systematic(order_by = "size"))
draw(df, walked, seed = 1)
#>    id site size
#> 1   5    a   84
#> 2  10    a   67
#> 3  17    a   23
#> 4  22    a    6
#> 5  23    a   43
#> 6  27    b   90
#> 7  33    b    9
#> 8  39    b   29
#> 9  40    b   66
#> 10 45    b   49
#> 11 51    c   69
#> 12 57    c   89
#> 13 63    c    8
#> 14 69    c   28
#> 15 75    c   48
#> 16 77    d   21
#> 17 83    d   41
#> 18 89    d   61
#> 19 90    d   98
#> 20 95    d   81
```
