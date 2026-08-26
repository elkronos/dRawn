# Estimate a population total from a sample

Forms the Horvitz-Thompson total `sum(y / pi)` and, where the design
allows it, a design-unbiased variance and confidence interval.

## Usage

``` r
ht_total(
  sample,
  y,
  variance = c("auto", "analytic", "jackknife", "none"),
  level = 0.95
)
```

## Arguments

- sample:

  A data frame returned by
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with
  `weights = TRUE`.

- y:

  The variable to total: a column name, or a numeric vector as long as
  `sample`.

- variance:

  How to compute it. `"auto"` uses the analytic estimator when the
  design has one and falls back to the jackknife when it does not;
  `"analytic"` insists on the analytic form, returning `NA` with the
  reason in `note` rather than falling back; `"jackknife"` always
  resamples; `"none"` skips it. The result reports which was used in
  `method` — and reports `"none"` when neither could produce a figure,
  rather than naming a method that declined.

- level:

  Confidence level for the interval.

## Value

A list with a [`print()`](https://rdrr.io/r/base/print.html) method,
holding:

- `total`:

  The Horvitz-Thompson total, `sum(y / pi)`.

- `variance`, `se`, `ci`, `level`:

  Its estimated variance, standard error and confidence interval. `NA`
  where the design supports none.

- `n`, `design`:

  Rows used, and the design's type.

- `deff`:

  The design effect — see
  [`deff()`](https://elkronos.github.io/dRawn/reference/deff.md).

- `method`:

  `"analytic"`, `"jackknife"` or `"none"`.

- `note`:

  Why a variance is missing, or which fallback was taken. `NULL` when
  the analytic estimator applied cleanly.

## Variance

The estimator is chosen to match how the design actually randomises:

- **Fixed-size designs** use the Sen-Yates-Grundy estimator, which is
  non-negative more often than the general Horvitz-Thompson form and is
  the usual choice.

- **Poisson sampling** has a random size and independent rows, so the
  independent-units form `sum((1 - pi) / pi^2 * y^2)` is used instead.

- **Cluster designs** take whole clusters, so the number of *rows* is
  random whenever the clusters differ in size. The cluster is the
  sampling unit, and the estimator is applied at that level —
  algebraically the same thing as the delete-a-cluster jackknife.

- **Certainty designs** hand the problem to `rest` over the rows below
  the threshold, since the certainty rows are in every possible sample
  and contribute nothing to the variance.

Sen-Yates-Grundy can still return a negative number on an unlucky
sample. That is a failure of the estimator rather than a variance, so it
is reported as one: `variance` is `NA`, `note` says what happened, and
`variance = "auto"` falls through to the jackknife.

A variance needs joint inclusion probabilities, and not every design has
them — see
[`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md).
Where they are unavailable the estimate is still returned, with
`variance` as `NA` and a note saying why. Systematic sampling is the
notable case: most pairs of rows can never co-occur, so no
design-unbiased variance estimator exists at all.

## See also

[`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md),
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)

## Examples

``` r
set.seed(1)
pop <- data.frame(
  id = 1:200,
  site = rep(c("a", "b"), times = c(150, 50)),
  spend = round(stats::runif(200, 10, 500))
)

s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
ht_total(s, "spend")
#> Horvitz-Thompson total  (stratified design, n = 40)
#>   estimate 54,000
#>   se       3,713.698  (analytic)
#>   95% CI  46,721.29 to 61,278.71
#>   deff     1.03  (about the same as simple random sampling)

sum(pop$spend)   # the truth
#> [1] 52732
```
