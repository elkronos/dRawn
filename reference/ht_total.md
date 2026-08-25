# Estimate a population total from a sample

Forms the Horvitz-Thompson total \`sum(y / pi)\` and, where the design
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

  A data frame returned by \[draw()\] with \`weights = TRUE\`.

- y:

  The variable to total: a column name, or a numeric vector as long as
  \`sample\`.

- variance:

  How to compute it. \`"auto"\` uses the analytic estimator when the
  design has one and falls back to the jackknife when it does not;
  \`"analytic"\` insists on the analytic form and errors otherwise;
  \`"jackknife"\` always resamples; \`"none"\` skips it. The result
  reports which was used.

- level:

  Confidence level for the interval.

## Value

A list with \`total\`, \`variance\`, \`se\`, \`ci\`, \`n\`, \`method\`
and \`note\`, with a \`print()\` method.

## Variance

For fixed-size designs the Sen-Yates-Grundy estimator is used, which is
non-negative more often than the general Horvitz-Thompson form and is
the usual choice. Poisson sampling has a random size, so the
independent-units form \`sum((1 - pi) / pi^2 \* y^2)\` is used instead.

A variance needs joint inclusion probabilities, and not every design has
them — see \[joint_prob()\]. Where they are unavailable the estimate is
still returned, with \`variance\` as \`NA\` and a note saying why.
Systematic sampling is the notable case: most pairs of rows can never
co-occur, so no design-unbiased variance estimator exists at all.

## See also

\[joint_prob()\], \[inclusion_prob()\]

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
#>   total    54,000
#>   se       3,713.698  (analytic)
#>   95% CI  46,721.29 to 61,278.71

sum(pop$spend)   # the truth
#> [1] 52732
```
