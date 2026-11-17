# Bootstrap sampling

Generates bootstrap replicates. `"simple"` resamples rows independently.
`"block"` is the moving-block bootstrap (Künsch 1989), which keeps the
serial dependence in ordered data by concatenating
`ceiling(n / block_length)` runs of `block_length` consecutive rows per
replicate, each starting at a row chosen at random. A run must fit
inside the frame, so fewer runs cover the first and last
`block_length - 1` rows, and the first row and the last turn up only
`1 / block_length` as often as a row in the middle. `"circular"`
(Politis and Romano 1992) joins the last row to the first, so a run can
wrap around the end and every row turns up equally often on average.

## Usage

``` r
design_bootstrap(
  n_replicates = 1000L,
  n = NULL,
  method = c("simple", "block", "circular"),
  block_length = NULL,
  output = c("rows", "indices", "weights")
)
```

## Arguments

- n_replicates:

  Number of replicates to generate.

- n:

  Rows per replicate. `NULL` uses `nrow(data)`, the standard
  nonparametric bootstrap. `n = 0` gives replicates with no rows. The
  frame then has no rows but keeps its `.replicate` column, so
  [`split()`](https://rdrr.io/r/base/split.html) returns no replicates
  at all, while
  `split(out, factor(out$.replicate, seq_len(n_replicates)))` returns
  `n_replicates` empty ones. A plan from
  [`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md)
  gives its `n`.

- method:

  `"simple"`, `"block"` or `"circular"`.

- block_length:

  Run length for `method = "block"` and `"circular"`, and an error for
  `method = "simple"`, where it would have no effect. `NULL` uses
  `floor(nrow(data)^(1/3))`, a common rule of thumb.

- output:

  What [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)
  returns.

  - `"rows"`, the default, returns the resampled rows of every replicate
    in one data frame, with a leading `.replicate` column.

  - `"indices"` returns an integer matrix with `n` rows and a column per
    replicate, holding the rows of `data` that each replicate drew, in
    order. Row `i` of replicate `r` is `data[out[i, r], ]`.

  - `"weights"` returns an integer matrix with a row for each row of
    `data` and a column per replicate, holding the number of times each
    replicate drew that row. With the default `n`, these are replicate
    weights for
    `survey::svrepdesign(type = "bootstrap", combined.weights = TRUE)`.

  For one seed, all three hold the same replicates.

## Value

A design object, for use with
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

## Details

By default
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) returns
all replicates in one data frame with a leading `.replicate` column.
Split them with `split(out, out$.replicate)`. That frame holds
`n_replicates * n` rows, a copy of the data for each replicate, and
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) warns,
with class `drawn_bootstrap_size`, when it would exceed 10 million rows.
`output = "indices"` and `output = "weights"` return the same replicates
as a matrix instead.

## Not a survey design

The bootstrap resamples the data you already have to approximate the
sampling distribution of a statistic (Efron 1979). It is not a
probability sample of a finite population, so it has no inclusion
probabilities and cannot be passed to
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md).
For the design-based variance of a total or mean, use the estimators
there.

Resampling rows treats them as independent draws from one distribution.
For such data and a statistic those estimators do not cover, such as a
median, the bootstrap is the standard tool. A sample drawn with strata
or clusters is not such data, and resampling its rows one by one ignores
the strata and treats the rows of a cluster as independent. Resample the
clusters within strata instead, for example with
`survey::as.svrepdesign(as_svydesign(s), type = "subbootstrap")`. Runs
of rows keep short-range dependence in ordered data (Künsch 1989).

## References

Efron, B. (1979). Bootstrap methods: another look at the jackknife.
*Annals of Statistics*, 7, 1–26.

Künsch, H. R. (1989). The jackknife and the bootstrap for general
stationary observations. *Annals of Statistics*, 17, 1217–1241.

Politis, D. N. and Romano, J. P. (1992). A circular block-resampling
procedure for stationary data. In R. LePage and L. Billard (eds.),
*Exploring the Limits of Bootstrap*, 263–270. Wiley.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)

Other designs:
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md),
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md),
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md),
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
[`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md),
[`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md),
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)

## Examples

``` r
df <- data.frame(id = 1:100, value = (1:100) / 10)

reps <- draw(df, design_bootstrap(n_replicates = 5, n = 20), seed = 1)
table(reps$.replicate)
#> 
#>  1  2  3  4  5 
#> 20 20 20 20 20 

# A statistic per replicate
vapply(split(reps, reps$.replicate), function(r) mean(r$value), numeric(1))
#>     1     2     3     4     5 
#> 5.450 5.285 4.395 4.735 5.570 

# The same replicates as row numbers, without copying the data
idx <- draw(df, design_bootstrap(n_replicates = 5, n = 20,
                                 output = "indices"), seed = 1)
apply(idx, 2, function(i) mean(df$value[i]))
#> [1] 5.450 5.285 4.395 4.735 5.570

# Runs of 10 rows that may wrap from the last row to the first
circ <- draw(df, design_bootstrap(n_replicates = 200, method = "circular",
                                  block_length = 10, output = "weights"),
             seed = 1)
range(rowMeans(circ))
#> [1] 0.910 1.185
```
