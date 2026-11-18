# Rebuild a weighted sample from its rows

A weighted sample carries its design, a copy of its frame and its
inclusion probabilities as attributes. A CSV or spreadsheet file, a
database table and [`merge()`](https://rdrr.io/r/base/merge.html) keep
none of them, so
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md),
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
and
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
refuse the rows that come back. `as_drawn()` puts the attributes back.
It matches each row of `data` to its row of `frame`, works out the
inclusion probabilities again from `design`, and sets `.prob`, `.weight`
and `.row` as
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) sets
them.

## Usage

``` r
as_drawn(data, design, frame, by = NULL)
```

## Arguments

- data:

  The sampled rows as read back: a data frame with the `.row` column
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) added,
  or with the id column or columns named in `by`.

- design:

  The design the sample was drawn with.

- frame:

  The data frame the sample was drawn from, as it was then.

- by:

  `NULL` to match rows through `.row`, or the names of one or more
  columns of `data` and `frame` that identify a frame row.

## Value

A data frame of the class of `data`, ready for
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md),
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
and
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md).

## Matching rows to the frame

With `by = NULL`, each row is matched through its `.row` column, the
frame row [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)
recorded for it. `frame` must then be the frame the sample was drawn
from, with its rows in the same order. A frame in another order, such as
a database table read without `ORDER BY`, holds other rows at those
positions, and its strata and clusters would be wrong for the sample.
Each column `data` shares with `frame` is therefore compared with the
frame rows `.row` names, and a column whose values differ is refused.
Numbers need to agree to about eight significant digits, fewer than a
CSV file keeps, and a number read back from text, such as `"007"`, is
compared as the number 7. Text such as `"T"`, `"F"`, `"t"` or `"f"` that
a file gave back as a logical column is compared as `TRUE` or `FALSE`,
as `readr::read_csv()` reads it. A date-time read back in another time
zone agrees when it is the same instant or shows the same clock time,
since `readr::read_csv()` and
[`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)
read the clock time a CSV file holds as UTC. A missing value on either
side is not compared, so a frame column that fieldwork fills in is
accepted. Nor is a date that a file gave back as text, whose format a
spreadsheet may have changed. A column changed after the draw, such as a
corrected address, is refused too; rename it in `data`, or match on an
id. A `data` that shares no column with `frame` besides those
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) adds
cannot be checked, so keep an id column in the file.

With `by`, rows are matched on the values of one or more id columns that
`data` and `frame` share. Values are compared as
[`match()`](https://rdrr.io/r/base/match.html) compares them, and as
text when the two columns have different classes, so an id stored as
`"007"` does not match the number 7 that a CSV file gives back.

When `data` has both `.row` and the `by` columns, they must agree. A
frame row found through `by` that is not the one `.row` names means
`frame` is not in the order it had when the sample was drawn, which
changes the probabilities of designs that depend on that order, such as
[`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md);
the rebuild is refused.

A row with a missing id, an id that appears twice in `data`, an id that
`frame` does not have or holds more than once, and a row the design
gives a probability of 0 are refused, with the ids named.

## Every row drawn

`data` must hold every row that was drawn, with or without a fieldwork
result. A subset, such as the rows that responded, is not a sample from
the design, and its standard errors would be wrong. For a design that
draws a fixed number of rows (simple, stratified, spread, temporal,
spatial and reservoir designs, systematic designs given `n`, and
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
with `method = "systematic"`, alone or as the `rest` of
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md)),
a `data` with another number of rows is refused. A
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)
or
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
sample, also as the `rest` of
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md),
is refused when its rows are not of `n_clusters` clusters or a cluster
does not hold the rows the design takes from it, which are all of its
rows for whole clusters, and the same number from each cluster when that
number does not depend on which clusters were selected. A design
weighted conditional on its clusters (see "Design weights" in
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)) is
rebuilt from the clusters its rows belong to, which must hold the rows
the design takes from each given that selection. A frame row the design
takes with certainty, with probability 1, is in every sample, and a
`data` without it is refused whatever the design. Keep the rows with no
response, with `y` missing, and estimate a mean with
`ht_mean(na_rm = TRUE)`; see "Missing values and nonresponse" in
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md).

A list from a design with `reserve` is rebuilt whole, the main sample
and every reserve, with the `.reserve` column
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) added,
which records the order the reserves are used in. Pass the result to
[`use_reserves()`](https://elkronos.github.io/dRawn/reference/use_reserves.md).

## What comes back

The rows of `data` in their own order, with `.prob`, `.weight` and
`.row` first and the other columns of `data` after them. `.prob` and
`.weight` come from `design`. A numeric `.prob` in `data` must agree
with them, or the design or the frame is not the one the sample was
drawn with. A `.weight` adjusted after the draw is replaced by
`1 / .prob`, so the adjustment has to be made again.

Strata, clusters and the other design columns are read from the copy of
`frame` stored with the result, as for a sample from
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md), so
`data` needs only the columns to estimate from: `y` and the `by` columns
of the estimates.
[`add_fieldwork()`](https://elkronos.github.io/dRawn/reference/add_fieldwork.md)
joins results that come back keyed by an id.

## See also

[`add_fieldwork()`](https://elkronos.github.io/dRawn/reference/add_fieldwork.md)
to join fieldwork results,
[`compact_sample()`](https://elkronos.github.io/dRawn/reference/compact_sample.md)
to drop frame columns before saving or sharing a sample, and
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

## Examples

``` r
pop <- data.frame(id = 1:200, site = rep(c("a", "b"), times = c(150, 50)))
by_site <- design_stratified("site", n = 40)
s <- draw(pop, by_site, seed = 1, weights = TRUE)

# A file keeps the rows and loses the design
path <- tempfile(fileext = ".csv")
write.csv(s, path, row.names = FALSE)
back <- as_drawn(read.csv(path), by_site, pop)
sample_summary(back)
#> Sample of 40 from 200  (stratified design)
#>   sampling fraction  0.2
#>   design weights     5 to 5   (cv 0)
#>   effective n        40  (Kish, from the weights alone)
#>   rebuilt            by as_drawn() of drawn 0.2.0.9000 on 2026-11-18; the seed is not recorded
#> 
#>   by site:
#>     group  drawn  in frame   rate
#>     a         30      150  0.200
#>     b         10       50  0.200

# Without `.row`, match on an id column
ids <- read.csv(path)[, c("id", "site")]
back <- as_drawn(ids, by_site, pop, by = "id")
unlink(path)
```
