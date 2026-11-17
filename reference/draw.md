# Draw a sample

Applies a sampling design to data.

## Usage

``` r
draw(data, design, seed = NULL, weights = FALSE)
```

## Arguments

- data:

  A data frame.
  [`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md)
  also accepts a list, a connection, or a zero-argument generator
  function.

- design:

  A design built by one of the
  [`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md)
  family.

- seed:

  Optional seed, applied only for this draw. See "Seeding".

- weights:

  Attach `.prob`, `.weight` and `.row` columns. See "Design weights".

## Value

The sampled rows. See "What you get back".

## What you get back

A data frame with the same class and the same columns, in the same
order, as `data`. `weights = TRUE` prepends three columns (see "Design
weights"), and some designs qualify it further:

- [`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md)
  prepends a `.replicate` column identifying which replicate each row
  belongs to, so all replicates come back in one frame. Split them with
  `split(out, out$.replicate)` if you need a list. With
  `output = "indices"` or `"weights"` it returns a matrix of row numbers
  or counts instead, with a column per replicate.

- [`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md)
  returns a list when `data` is a list or a stream rather than a data
  frame, because there is nothing to make a data frame out of. With
  `batches = TRUE` it returns the rows sampled from a function's data
  frames as one data frame.

- [`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md)
  and
  [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
  with `reserve` return a list: the main sample, then the reserves, with
  a leading integer `.reserve` column that is 0 for the main sample and
  numbers the reserves 1, 2, ... in the order they are used. With
  `weights = TRUE` it follows `.row`, and `.prob` is the main sample's
  probability for every row.
  [`use_reserves()`](https://elkronos.github.io/dRawn/reference/use_reserves.md)
  turns the list into the sample to estimate from.

Row names of the frame's own, such as the car names of `mtcars`, stay
with their rows. A row drawn more than once, by a design with
`replace = TRUE` or by
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
gets names made unique by
[`make.unique()`](https://rdrr.io/r/base/make.unique.html):
`"Mazda RX4"`, then `"Mazda RX4.1"`. The automatic row names 1, 2, ...
number the sample's rows afresh, and `weights = TRUE` records each row's
place in `data` in `.row`. A tibble has no row names to keep.

## Design weights

`weights = TRUE` prepends three columns: `.prob`, the probability that
the row was included in the sample; `.weight`, its reciprocal, the
number of population units the row stands for; and `.row`, the row of
`data` it came from. `.prob` and `.weight` are what an unbiased estimate
needs: a Horvitz-Thompson total is `sum(y * .weight)`.

They come from the design applied to the population, not from the drawn
sample, and some designs have no closed form for them.
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
documents which, why, and what to do instead.

Two kinds of design are weighted conditional on the clusters the draw
selected:
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)
with `balanced = TRUE`, and
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
without `size` when the rows taken from a cluster depend on which other
clusters were selected. With `a` of `A` clusters selected and `m_i` of
the `N_i` rows of cluster `i` taken, `.prob` is `(a / A) * m_i / N_i`,
the chance of the cluster being selected times the chance of the row
being drawn given these clusters. The total estimated from it is
unbiased given the clusters selected, and so over every selection.
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
says when a sample is weighted this way.

Keep `.prob` as it is:
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md) and
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
refuse a sample whose `.prob` was changed. Weights adjusted after the
draw, for nonresponse say, go in `.weight`, and
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
is the way to estimate from them;
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
and [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
keep to `.prob` and warn when `.weight` no longer matches it.

The sample also carries its design, the frame rows it came from and a
copy of the frame as attributes, which is what lets
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md),
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
and
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
work from the sample alone.

- The frame is copied when the sample is drawn, every column of it.
  Later edits to `data`, including data.table's in-place `:=`,
  `setorder()` and `setkey()`, do not reach the sample. Each sample
  holds its own copy, in memory and in every file it is saved to with
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html), so a list of
  samples from a large frame holds the frame once per sample, and a
  sample sent to someone sends them the frame, with any personal data in
  it.
  [`compact_sample()`](https://elkronos.github.io/dRawn/reference/compact_sample.md)
  keeps only the columns the design reads.

- Reordering the sample by `[`, `dplyr::arrange()` or
  [`data.table::setorder()`](https://rdrr.io/pkg/data.table/man/setorder.html)
  is safe, because each row finds its frame row through `.row`.

- Estimate from the whole sample and use `by` for groups. Row subsets
  ([`head()`](https://rdrr.io/r/utils/head.html), `s[cond, ]`,
  `dplyr::filter()`), [`rbind()`](https://rdrr.io/r/base/cbind.html) of
  samples and joins that repeat rows keep the attributes but no longer
  match the rows that were drawn, and the functions above refuse them.
  [`merge()`](https://rdrr.io/r/base/merge.html), column selection
  (`s[, cols]`), [`subset()`](https://rdrr.io/r/base/subset.html) and
  [`transform()`](https://rdrr.io/r/base/transform.html) drop the
  attributes altogether. Add fieldwork results with
  [`add_fieldwork()`](https://elkronos.github.io/dRawn/reference/add_fieldwork.md),
  which matches them on an id, or with `s$col <- value`.

- A CSV or spreadsheet file and a database table keep the rows and the
  `.row` column but not the attributes.
  [`as_drawn()`](https://elkronos.github.io/dRawn/reference/as_drawn.md)
  rebuilds the sample from them, the design and the frame.

`weights = TRUE` is not available for a design that samples **with
replacement**, or a
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md)
whose `rest` does. `.prob` there is the probability of being selected
*at least once*, but the sample holds duplicates, so `sum(y * .weight)`
over it counts each duplicate at the distinct-unit weight. The total
comes out high by a factor of `(n / N) / (1 - (1 - 1 / N)^n)`, which is
about 5% at `n = N / 10`, 15% at `n = 0.3 * N`, 58% at `n = N` and more
beyond, applied within each stratum of a stratified design. Use
`replace = FALSE`, or the Hansen-Hurwitz form `N / n * sum(y)`, summed
over strata as `N_h / n_h * sum(y_h)` for a stratified design.

## Seeding

With a `seed`, `draw()` saves the caller's random number state, seeds,
draws, and restores the state on exit, so drawing a sample inside a
simulation does not shift the simulation's own random number stream.
With `seed = NULL` it draws from the caller's stream and advances it, as
any function that draws random numbers does.

A seed is applied with R's default generator,
`set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")`,
whatever generator the session has selected, and the session's generator
is restored on exit with its state. A seed therefore draws the same
sample after `RNGkind("L'Ecuyer-CMRG")`, which parallel code often
selects, or under `sample.kind = "Rounding"` as it does in a session
left at R's defaults. It must be a whole number from
`-.Machine$integer.max` to `.Machine$integer.max`. The same holds for
the `seed` of
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md),
[`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)
and
[plot()](https://elkronos.github.io/dRawn/reference/plot.drawn_design.md).
Most designs select rows by their position in `data`, so a seed is sure
to draw the same rows again only from a frame with the same rows in the
same order.

With `weights = TRUE` the sample records how it was drawn in its
`drawn_meta` attribute: `seed` (`NULL` without one), `rng_kind` (the
generator used, as [`RNGkind()`](https://rdrr.io/r/base/Random.html)
gives it), `version` (of drawn), `time`, `frame_rows` (the rows of
`data`) and `made_by` (`"draw()"`, or `"as_drawn()"` for a sample
rebuilt from a file).
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
prints it.

## See also

[designs](https://elkronos.github.io/dRawn/reference/designs.md) for the
full list of constructors,
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
for the probabilities themselves, and
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
to estimate a population total with a standard error.

## Examples

``` r
df <- data.frame(id = 1:100, site = rep(letters[1:4], each = 25))

draw(df, design_simple(n = 10), seed = 1)
#>    id site
#> 1  68    c
#> 2  39    b
#> 3   1    a
#> 4  34    b
#> 5  87    d
#> 6  43    b
#> 7  14    a
#> 8  82    d
#> 9  59    c
#> 10 51    c

# The same design, reused
by_site <- design_stratified(strata = "site", n = 12)
table(draw(df, by_site, seed = 1)$site)
#> 
#> a b c d 
#> 3 3 3 3 
table(draw(df[1:60, ], by_site, seed = 1)$site)
#> 
#> a b c 
#> 5 5 2 

# With design weights, ready for estimation
df$spend <- seq_len(100)
s <- draw(df, by_site, seed = 1, weights = TRUE)
sum(s$spend * s$.weight)   # estimates sum(df$spend) = 5050
#> [1] 4966.667
```
