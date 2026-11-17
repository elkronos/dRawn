# Replace unavailable units from a reserve list

[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md)
and
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
with `reserve` give a list: the `n` units of the main sample, then the
reserves, with a `.reserve` column that numbers them 1, 2, ... in the
order they are used. Visit every unit of the main sample, and for each
one found unavailable, visit the next reserve, until the units found
available make up `n`. `use_reserves()` takes the list and what the
visits found, and returns the sample to estimate from: the units of the
main sample found available, and as many reserves, taken in order, as
replace the others.

## Usage

``` r
use_reserves(sample, unavailable)
```

## Arguments

- sample:

  A list from
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with
  `weights = TRUE` and a design with `reserve`, or from
  [`as_drawn()`](https://elkronos.github.io/dRawn/reference/as_drawn.md).

- unavailable:

  A logical vector with one value for each row of `sample`, in its
  order, or the name of such a column of `sample`. See "Marking the
  visits".

## Value

A weighted sample, ready for
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md),
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
and
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md).

## What the estimates describe

Replacing units changes the population being estimated. The sample that
comes back stands for the units of the frame that can be reached, and
its totals and means are theirs, not the whole frame's. It is weighted
as a sample of `n` drawn at equal probability from the frame without the
`u` listed units found unavailable: each unit used has probability
`n / (N - u)`, and the units found unavailable have probability 0, so
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
and [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
count them among the frame rows the estimate does not cover.

The probabilities assume that every unit of the frame nobody visited is
available. A frame that holds more unavailable units than were found has
fewer than `N - u` units that can be reached, and the totals are too
high by about the share of the frame those units make up. The Hajek mean
of [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
does not depend on `N - u` and is not affected. An estimate from the
sample keeps `u` in `unavailable`, and its printed `unreachable` line
adds that the units nobody visited are taken as available.

For
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md)
the weighting is exact under that assumption, when the reserves suffice.
The list is one random draw in random order, so its first `n` available
units are a simple random sample of `n` from the available units of the
frame, and the variance is that of a simple random sample from `N - u`.
For
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
the reserves are a second spread sample from the rows the main sample
left, in random order. The units used are spread over the space, and
under the same assumption their probabilities average `n / (N - u)`, but
each unit's own probability is not known exactly, so the weights are an
approximation. The variance is the local mean estimator over the units
used.

Reserves are for units that cannot be reached or measured at all. A unit
that is reached but gives no answer is nonresponse: keep it in the
sample with `y` missing, as "Missing values and nonresponse" in
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
describes.

## Marking the visits

`unavailable` holds one value for each row of `sample`: `TRUE` for a
unit found unavailable, `FALSE` for one found available, and `NA` for
one not visited. Every unit of the main sample needs `TRUE` or `FALSE`,
and so does each reserve up to the one that completes the sample; a
missing value there is refused, with the frame row named. Reserves
visited past that one are not used, and those of them found unavailable
count among the `u` units outside the population. When the reserves run
out first, the sample holds fewer than `n` units and is weighted as a
sample of that many, with a warning.

## What comes back

The units used, in the order of `sample`, with `.prob`, `.weight`,
`.row` and `.reserve` first. `.reserve` is 0 for a unit of the main
sample and the reserve's number for a unit that replaced one. The design
stored with the result is the list's design with `n` the number of units
used, no `reserve`, and the frame rows found unavailable recorded as
`unavailable`;
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
of it gives those rows probability 0.

## See also

[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md),
[`add_fieldwork()`](https://elkronos.github.io/dRawn/reference/add_fieldwork.md)
to add what the visits found to the list.

## Examples

``` r
pop <- data.frame(id = 1:200, spend = (1:200) %% 37)
s <- draw(pop, design_simple(n = 20, reserve = 5), seed = 1, weights = TRUE)
table(s$.reserve)
#> 
#>  0  1  2  3  4  5 
#> 20  1  1  1  1  1 

# Two units of the main sample could not be reached, nor could the first
# reserve; reserves 2 and 3 replaced them, and 4 and 5 were not visited
gone <- c(s$id[c(3, 11)], s$id[s$.reserve == 1])
found <- ifelse(s$.reserve <= 3, s$id %in% gone, NA)
used <- use_reserves(s, found)
table(used$.reserve)
#> 
#>  0  2  3 
#> 18  1  1 
ht_total(used, "spend")
#> Horvitz-Thompson total  (simple design, n = 20)
#>   estimate 3,319.45
#>   se       514.6565  (analytic)
#>   95% CI  2,242.261 to 4,396.639  (t, 19 df)
#>   deff     1  (about the same as simple random sampling)
#>   unreachable  3 frame rows at probability 0, not covered by this estimate
#>                (3 found unavailable; rows not visited are taken as available)
```
