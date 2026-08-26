# Joint inclusion probabilities

The probability that a *pair* of rows both land in the sample.
First-order probabilities from
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
give you an unbiased total; second-order probabilities are what let you
put a standard error on it.

## Usage

``` r
joint_prob(data, design, rows = NULL, simulate = FALSE, R = 5000, seed = NULL)
```

## Arguments

- data:

  A data frame.

- design:

  A design object.

- rows:

  Optional row indices. Supply these — usually the rows you drew — to
  get the submatrix for them instead of the full `nrow(data)` square,
  which is what makes this usable on a large population.

- simulate:

  Estimate the probabilities by repeated draws rather than in closed
  form. Works for every probability design, including the ones with no
  closed form, at the cost of Monte Carlo error. Refused for
  [`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
  where every row appears in some replicate and the count converges to 1
  for all of them.

- R:

  Number of simulated draws when `simulate = TRUE`.

- seed:

  Optional seed for the simulation.

## Value

A square matrix with one row and column per element of `rows` (or per
row of `data`). The diagonal holds first-order probabilities.

## Which designs have them

Closed forms exist, and are used, for the designs whose selection is
either independent across groups or a simple random sample within them:

|  |  |
|----|----|
| [`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md) | `n(n-1) / (N(N-1))` for a pair, without replacement |
| [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md) | within a stratum as above; across strata, independent |
| [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md) | same cluster: `a/A`; different clusters: `a(a-1)/(A(A-1))` |
| [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md) | the two stages multiplied, where the per-cluster take is constant |
| [`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md) | `1` between certainty rows; otherwise the other row's own `pi` |
| [`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md) | simple random sampling over the first `max_items` rows |
| [`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md) | within an interval as above; across intervals, independent |
| [`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md) | simple random sampling inside the region |
| [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md) | `"poisson"` only, where rows are independent: `pi_i * pi_j` |
| [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md) | `1/interval` for rows sharing a residue class, otherwise **zero** |

Systematic sampling is the awkward one. Most pairs can never co-occur,
so their joint probability is genuinely 0 and no design-unbiased
variance estimator exists.
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
says so rather than returning a number. Its residue classes follow the
order the design walks, so `order_by` changes which pairs can co-occur.

`design_weighted(method = "systematic")` has joint probabilities, but
they depend on the order units are visited and need a dedicated
algorithm. Use `sampling::UPsystematicpi2()` for those.

## See also

[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md),
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)

## Examples

``` r
df <- data.frame(id = 1:10)
round(joint_prob(df, design_simple(n = 4)), 3)
#>        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
#>  [1,] 0.400 0.133 0.133 0.133 0.133 0.133 0.133 0.133 0.133 0.133
#>  [2,] 0.133 0.400 0.133 0.133 0.133 0.133 0.133 0.133 0.133 0.133
#>  [3,] 0.133 0.133 0.400 0.133 0.133 0.133 0.133 0.133 0.133 0.133
#>  [4,] 0.133 0.133 0.133 0.400 0.133 0.133 0.133 0.133 0.133 0.133
#>  [5,] 0.133 0.133 0.133 0.133 0.400 0.133 0.133 0.133 0.133 0.133
#>  [6,] 0.133 0.133 0.133 0.133 0.133 0.400 0.133 0.133 0.133 0.133
#>  [7,] 0.133 0.133 0.133 0.133 0.133 0.133 0.400 0.133 0.133 0.133
#>  [8,] 0.133 0.133 0.133 0.133 0.133 0.133 0.133 0.400 0.133 0.133
#>  [9,] 0.133 0.133 0.133 0.133 0.133 0.133 0.133 0.133 0.400 0.133
#> [10,] 0.133 0.133 0.133 0.133 0.133 0.133 0.133 0.133 0.133 0.400

# Only for the rows you drew
joint_prob(df, design_simple(n = 4), rows = c(2, 5, 7))
#>           [,1]      [,2]      [,3]
#> [1,] 0.4000000 0.1333333 0.1333333
#> [2,] 0.1333333 0.4000000 0.1333333
#> [3,] 0.1333333 0.1333333 0.4000000
```
