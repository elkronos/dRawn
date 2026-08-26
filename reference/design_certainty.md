# Take some rows with certainty, sample the rest

Composes a design out of two parts: rows at or above a threshold are all
taken, and any other design is applied to what remains. This is the
standard shape of an audit or financial sample — every invoice over a
hundred thousand gets examined, the long tail below is sampled — and it
is the cheapest way to cut variance when a few units dominate the total.

## Usage

``` r
design_certainty(above, threshold, rest, na_rm = FALSE)
```

## Arguments

- above:

  Column holding the size measure.

- threshold:

  Rows with `above >= threshold` are taken with certainty.

- rest:

  A design applied to the rows below the threshold. Its `n` is the
  number drawn from *those* rows, not the total.

- na_rm:

  Drop rows whose size measure is `NA` instead of raising an error.

## Value

A design object, for use with
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

## Details

Certainty rows have an inclusion probability of exactly 1, so they
contribute their own value to a Horvitz-Thompson total with a weight of
1 and add nothing to its variance. That is the whole point: the
uncertainty in the estimate comes only from the part you sampled.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md),
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)
for probability-proportional-to-size, which handles dominant units by
taking them with certainty automatically.

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md),
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md),
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
[`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md),
[`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md),
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)

## Examples

``` r
set.seed(1)
invoices <- data.frame(
  id = 1:500,
  site = rep(c("a", "b"), times = c(300, 200)),
  value = round(stats::rlnorm(500, 7, 1.4))
)

# Everything over 20,000 examined; 40 drawn from the rest
d <- design_certainty("value", threshold = 20000,
                      rest = design_stratified("site", n = 40))
s <- draw(invoices, d, seed = 1, weights = TRUE)

table(certain = s$.prob == 1)
#> certain
#> FALSE  TRUE 
#>    40    13 
sum(invoices$value >= 20000)
#> [1] 13
```
