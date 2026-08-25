# Joint inclusion probabilities

The probability that a \*pair\* of rows both land in the sample.
First-order probabilities from \[inclusion_prob()\] give you an unbiased
total; second-order probabilities are what let you put a standard error
on it.

## Usage

``` r
joint_prob(data, design, rows = NULL)
```

## Arguments

- data:

  A data frame.

- design:

  A design object.

- rows:

  Optional row indices. Supply these — usually the rows you drew — to
  get the submatrix for them instead of the full \`nrow(data)\` square,
  which is what makes this usable on a large population.

## Value

A square matrix with one row and column per element of \`rows\` (or per
row of \`data\`). The diagonal holds first-order probabilities.

## Which designs have them

Closed forms exist, and are used, for the designs whose selection is
either independent across groups or a simple random sample within them:

|  |  |
|----|----|
| \[design_simple()\] | \`n(n-1) / (N(N-1))\` for a pair, without replacement |
| \[design_stratified()\] | within a stratum as above; across strata, independent |
| \[design_cluster()\] | same cluster: \`a/A\`; different clusters: \`a(a-1)/(A(A-1))\` |
| \[design_multistage()\] | the two stages multiplied, equal allocation only |
| \[design_reservoir()\] | identical to simple random sampling |
| \[design_temporal()\] | within an interval as above; across intervals, independent |
| \[design_spatial()\] | simple random sampling inside the region |
| \[design_weighted()\] | \`"poisson"\` only, where rows are independent: \`pi_i \* pi_j\` |
| \[design_systematic()\] | \`1/interval\` for rows sharing a residue class, otherwise \*\*zero\*\* |

Systematic sampling is the awkward one. Most pairs can never co-occur,
so their joint probability is genuinely 0 and no design-unbiased
variance estimator exists. \[ht_total()\] says so rather than returning
a number.

\`design_weighted(method = "systematic")\` has joint probabilities, but
they depend on the order units are visited and need a dedicated
algorithm. Use \`sampling::UPsystematicpi2()\` for those.

## See also

\[inclusion_prob()\], \[ht_total()\]

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
