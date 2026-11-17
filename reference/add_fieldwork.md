# Join fieldwork results to a weighted sample

Fieldwork results usually come back as a table of their own, keyed by an
id and in another order. `add_fieldwork()` adds its columns to `sample`,
row by row through the id, and keeps the sample's rows, their order and
the attributes that hold its design.
[`merge()`](https://rdrr.io/r/base/merge.html) reorders the rows and
drops those attributes, and a join that matches a row twice would repeat
it.

## Usage

``` r
add_fieldwork(sample, field, by)
```

## Arguments

- sample:

  A data frame returned by
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with
  `weights = TRUE`, or by
  [`as_drawn()`](https://elkronos.github.io/dRawn/reference/as_drawn.md).

- field:

  A data frame of results: the `by` columns and the columns to add.

- by:

  The names of one or more columns of `sample` and `field` that identify
  a sampled row.

## Value

`sample` with the other columns of `field` added after its own.

## Details

A sampled row with no row in `field` gets missing values, as a
nonrespondent should:
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
refuses them, and `ht_mean(na_rm = TRUE)` estimates the mean over the
rows with a result (see "Missing values and nonresponse" in
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)).
An id that appears twice in `field`, an id that is not in `sample`, a
missing id, and a column of `field` that `sample` already has are
refused. Ids are compared as in
[`as_drawn()`](https://elkronos.github.io/dRawn/reference/as_drawn.md).

A list from a design with `reserve` takes fieldwork columns too, such as
which units were found unavailable, before
[`use_reserves()`](https://elkronos.github.io/dRawn/reference/use_reserves.md)
turns it into the sample to estimate from.

## See also

[`as_drawn()`](https://elkronos.github.io/dRawn/reference/as_drawn.md),
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md),
[`use_reserves()`](https://elkronos.github.io/dRawn/reference/use_reserves.md)

## Examples

``` r
pop <- data.frame(id = 1:200, site = rep(c("a", "b"), times = c(150, 50)))
s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)

# Results for 36 of the 40 rows, in the order they came in
field <- data.frame(id = rev(s$id)[1:36], spend = seq(20, 370, by = 10))
s <- add_fieldwork(s, field, by = "id")
ht_mean(s, "spend", na_rm = TRUE)
#> Hajek mean  (stratified design, n = 36)
#>   missing  4 sampled rows left out, assumed missing completely at random
#>   estimate 195
#>   se       10.22661  (analytic)
#>   95% CI  174.2973 to 215.7027  (t, 38 df)
```
