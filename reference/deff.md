# Design effect

How much precision the design costs against simple random sampling of
the same size (Kish 1965). `deff = 1` means the design is doing as well
as a simple random sample from the rows it can reach; above 1 it is
doing worse, which is the usual price of clustering; below 1 it is doing
better, which is what stratification and
probability-proportional-to-size buy you.

## Usage

``` r
deff(object, ...)
```

## Arguments

- object:

  A result from
  [`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
  [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md),
  [`ht_ratio()`](https://elkronos.github.io/dRawn/reference/ht_ratio.md)
  or
  [`ht_quantile()`](https://elkronos.github.io/dRawn/reference/ht_quantile.md),
  or an object from `survey` that holds a design effect.

- ...:

  Passed to
  [`survey::deff()`](https://rdrr.io/pkg/survey/man/surveysummary.html)
  for an object from `survey`, such as its `quietly`. Not used for an
  estimate made here.

## Value

A single number, or `NA` when the design has no variance estimate, `y`
is constant in the sample, or the estimate is a mean from
`ht_mean(na_rm = TRUE)` with `y` missing in some sampled rows. For an
object from `survey`, what
[`survey::deff()`](https://rdrr.io/pkg/survey/man/surveysummary.html)
returns.

## Details

Read it as an exchange rate on sample size: at `deff = 2`, a sample of
400 carries about as much information as 200 drawn at random.

The simple random sample it is measured against is drawn from the rows
the design gives a probability above 0,
`sum(inclusion_prob(design, frame) > 0)` of them, because those are the
rows a total or a Hajek mean estimates over. Rows past `max_items`,
outside a region or time window, dropped by `na_rm` or in a stratum
allocated no rows are left out. Counting them set an equal-probability
sample of 20 from the first 100 rows of a 400-row frame at a design
effect of 0.05 for a total, "better than simple random sampling", where
the right figure is 1. `survey` measures against the sum of the weights,
which is that count for an equal-probability design, and gives the same
design effect there. A Horvitz-Thompson mean (`estimator = "ht"`)
divides the total by the frame's row count, so its design effect is the
total's.

The design effect is `NA` when there is no variance, and when `y` takes
one value in every sampled row, where there is no spread for a simple
random sample to be measured by. It is `NA` too for a mean from
`ht_mean(na_rm = TRUE)` with `y` missing in some sampled rows, since the
rows with `y` are a domain of the sample, and a domain's design effect
is not estimated.

A table of domain estimates from `by`, of several variables, of the
levels of a factor or of quantiles carries no design effect, and
`deff()` refuses it with a message saying so. A contrast's design effect
is `NA`, and so is a quantile's from
[`ht_quantile()`](https://elkronos.github.io/dRawn/reference/ht_quantile.md).
A ratio from
[`ht_ratio()`](https://elkronos.github.io/dRawn/reference/ht_ratio.md),
and a total from the ratio estimator, are measured against the same
estimator from a simple random sample, through the variance of its
linearised residuals, as `survey::svyratio(deff = TRUE)` measures a
ratio.

`survey` has a `deff()` of its own, and whichever of the two packages is
attached last masks the other's. Either one answers for both. This
`deff()` is a generic that passes an object it has no method for, such
as `survey::svymean(..., deff = TRUE)`, to
[`survey::deff()`](https://rdrr.io/pkg/survey/man/surveysummary.html),
and `survey`'s generic finds the methods for the estimates made here.

## References

Kish, L. (1965). *Survey Sampling*. Wiley.

## See also

[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md),
[`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md),
which takes a design effect as an input.

## Examples

``` r
set.seed(1)
# Sites differ from each other, and rows within a cluster are alike --
# exactly the structure that makes stratifying pay and clustering cost.
pop <- data.frame(
  id = 1:400,
  site = rep(c("a", "b", "c", "d"), each = 100),
  cl = rep(paste0("c", 1:40), each = 10)
)
pop$y <- rep(c(20, 60, 120, 200), each = 100) + round(stats::rnorm(400, 0, 8))

# Stratifying on something that matters buys precision (deff below 1)
deff(ht_total(draw(pop, design_stratified("site", n = 40), seed = 1,
                   weights = TRUE), "y"))
#> [1] 0.01174153

# Clustering usually costs it (deff above 1)
deff(ht_total(draw(pop, design_cluster("cl", n_clusters = 4), seed = 1,
                   weights = TRUE), "y"))
#> [1] 12.91977
```
