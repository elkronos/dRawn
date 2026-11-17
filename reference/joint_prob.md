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

  Optional row indices, usually the rows you drew, which are `s$.row`
  for a sample `s` from `draw(weights = TRUE)`. The result is then the
  submatrix for these rows instead of the full `nrow(data)` square,
  which is what makes this usable on a large population. See "Size".
  Indices must be whole numbers. A row given more than once is one unit,
  so its entry with itself, off the diagonal too, is its own inclusion
  probability.

- simulate:

  Estimate the probabilities by repeated draws rather than in closed
  form. Works for every probability design, including the ones with no
  closed form, at the cost of Monte Carlo error. Refused for
  [`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
  where every row appears in some replicate and the count converges to 1
  for all of them.

- R:

  Number of simulated draws when `simulate = TRUE`. A simulation costs
  `R` draws; see "What a simulation costs" in
  [`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md).

- seed:

  Optional seed for the simulation, applied as in
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) (see
  "Seeding" there). Checked whether or not `simulate` is `TRUE`.

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
| [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md) | the two stages multiplied, where the per-cluster take is constant and `size` is not given |
| [`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md) | `1` between certainty rows; otherwise the other row's own `pi` |
| [`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md) | simple random sampling over the first `max_items` rows |
| [`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md) | within an interval as above; across intervals, independent |
| [`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md) | simple random sampling inside the region |
| [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md) | `"poisson"` only, where rows are independent: `pi_i * pi_j` |
| [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md) | `1/interval` for rows sharing a residue class, otherwise **zero**; with `n`, the overlap of the two rows' stretches, from 0 to `n / N`; with `n_starts = m`, `1/interval` for rows sharing a residue class modulo `m * interval`, otherwise `m(m-1) / (mk(mk-1))` with `k` the interval |

Systematic sampling is the awkward one. At a step of 2 rows or more,
which is every `interval` above 1 and every `n` up to half the rows
walked, neighbouring rows and many other pairs can never co-occur, so
their joint probability is genuinely 0 and no design-unbiased variance
estimator exists;
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
uses the successive-difference approximation instead, and says so. Its
residue classes follow the order the design walks, so `order_by` changes
which pairs can co-occur. Given `n`, the design walks at a step of
`N / n` rows from a start `u` uniform on (0, 1), each row covering a
stretch of length `n / N`, and a pair is drawn together when `u` falls
where both stretches, moved into (0, 1), overlap. With `n` above `N / 2`
the stretches are longer than 1/2, any two of them overlap, and every
pair has a positive joint probability;
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
uses the successive-difference approximation there too. With a fixed
`start` the design is deterministic, and it is refused as in
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md).
With several starts (`n_starts`), rows in different residue classes can
be drawn together too, which is what gives a replicated sample its
unbiased variance. Within the strata of
`design_stratified(within = design_systematic())`, pairs in one stratum
get the overlap of their stretches, as with `n`, and pairs in two strata
the product of their probabilities.

Three designs have joint probabilities that are positive but have no
closed form: `design_weighted(method = "systematic")`, which shuffles
the rows before walking them so every pair can co-occur, with or without
`strata`;
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
with `size`, which does the same with the clusters; and
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md).
For all three,
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
uses a variance approximation that needs only first-order probabilities,
and `simulate = TRUE` here estimates the joint probabilities themselves.

## Size

The matrix has a row and a column for each of `rows`, by default every
row of `data`, and holds 8 bytes for each pair: 800 MB for 10,000 rows.
Past about 1 GB, from about 11,200 rows, `joint_prob()` refuses and
gives the size, rather than risk exhausting memory, which can end the R
session.
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
and [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
do not build this matrix for any design in the package.

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
