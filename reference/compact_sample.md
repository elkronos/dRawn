# Drop the frame columns a sample's design does not need

A weighted sample holds its own copy of the whole frame, every column of
it, so that
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md),
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
and
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
can work from the sample alone. That copy is kept in memory and written
to every file the sample is saved to, which costs as much as the frame
each time, and hands the frame, with any names, addresses or other
personal data in it, to whoever receives the file. `compact_sample()`
replaces the copy with one that holds only the columns the design reads
(strata, clusters, size measures, `order_by`, time stamps, coordinates
and the like) and the columns named in `keep`.

## Usage

``` r
compact_sample(sample, keep = NULL)
```

## Arguments

- sample:

  A data frame returned by
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with
  `weights = TRUE`, or by
  [`as_drawn()`](https://elkronos.github.io/dRawn/reference/as_drawn.md).

- keep:

  Names of other columns of the frame to keep, such as the `by` columns
  of the estimates.

## Value

`sample`, with a smaller copy of the frame.

## Details

Estimates, standard errors,
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
and
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
are unchanged, with two exceptions that read other frame columns. A `by`
domain's size in the frame, which `ht_mean(estimator = "ht")` with `by`
divides by, needs the `by` columns in the stored frame; without them
that estimator is refused with a message naming `keep`. And a stratum
with a single sampled row leaves a domain's variance unmeasured unless
the stored frame's `by` columns show the stratum holds none of the
domain (see "Domains" in
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)).
Keep the `by` columns of the estimates in `keep`, and the columns
[`frame_totals()`](https://elkronos.github.io/dRawn/reference/frame_totals.md)
is asked for. Whether a column of the frame left out holds only 0 and 1,
which makes
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
treat it as a proportion, is recorded.

The sample's own columns are not changed: drop columns of the sample
with `s$col <- NULL`. The inclusion probabilities of every frame row are
kept, one number per row.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md),
[`as_drawn()`](https://elkronos.github.io/dRawn/reference/as_drawn.md)

## Examples

``` r
pop <- data.frame(id = 1:5000, site = rep(c("a", "b"), each = 2500),
                  note = strrep("x", 40), spend = (1:5000) %% 97)
s <- draw(pop, design_stratified("site", n = 100), seed = 1, weights = TRUE)
small <- compact_sample(s)

# The same estimate from a fraction of the bytes
c(full = length(serialize(s, NULL)),
  compact = length(serialize(small, NULL)))
#>    full compact 
#>  395577   95655 
all.equal(ht_total(s, "spend")$se, ht_total(small, "spend")$se)
#> [1] TRUE
```
