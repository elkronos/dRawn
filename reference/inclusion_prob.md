# Inclusion probabilities and design weights

`inclusion_prob()` returns the first-order inclusion probability of
every row of `data` under `design` — the probability that the row lands
in a sample — without drawing one. `sampling_weight()` returns
`1 / inclusion_prob()`, the number of population units each sampled row
stands for.

## Usage

``` r
inclusion_prob(data, design, simulate = FALSE, R = 5000, seed = NULL)

sampling_weight(data, design, simulate = FALSE, R = 5000, seed = NULL)
```

## Arguments

- data:

  A data frame.

- design:

  A design object.

- simulate:

  Estimate the probabilities by repeated draws rather than in closed
  form. Required for the designs listed above; allowed for any
  probability design, which is a convenient way to check the exact
  formulas.

- R:

  Number of simulated draws when `simulate = TRUE`.

- seed:

  Optional seed for the simulation.

## Value

A numeric vector with one element per row of `data`.

Rows the design can never select — outside the region, outside the time
window, at zero weight — get `0` from `inclusion_prob()` and `NA` from
`sampling_weight()`, since no finite number of population rows is
represented by a row that cannot be drawn.
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
counts them.

## Details

These are what make a sample usable for estimation. A Horvitz-Thompson
total is `sum(y / pi)` over the sampled rows, and its unbiasedness rests
entirely on `pi` being the design's real inclusion probability rather
than a plausible guess.

## Which designs have a closed form

Most do, and those are computed exactly:

|  |  |
|----|----|
| [`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md) | `n / N` |
| [`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md) | `n_h / N_h` within each stratum |
| [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md) | `1 / interval`, for every row, exactly |
| [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md) | `n_clusters / N_clusters` |
| [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md) | `(n_clusters / N_clusters) * (n_h / N_h)`, when the per-cluster take is constant |
| [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md) | `n * p_i` for the `"systematic"` and `"poisson"` methods |
| [`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md) | `1` above the threshold, `rest`'s own below it |
| [`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md) | `n / min(N, max_items)`, and `0` past `max_items` |
| [`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md) | `per_interval / N_bucket` within each interval |
| [`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md) | `n / N_in_region` |

Five cases have no closed form, and the package refuses to invent one:

- `design_cluster(balanced = TRUE)` — the per-cluster take is the
  smallest *selected* cluster's size, which is itself random. Simulation
  on clusters of 2/4/6/8 gives 0.50, 0.41, 0.34, 0.25, against the 0.50
  a naive `n_clusters / N_clusters` would claim for every row.

- `design_multistage(allocation = "proportional")` — the stage-two
  allocation depends on which clusters were selected. Simulation on
  clusters of 3/5/7/9 gives 0.17, 0.20, 0.17, 0.15, against a naive
  0.33, 0.20, 0.14, 0.11.

- [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
  where the per-cluster take is not constant — `n` not divisible by
  `n_clusters`, or a cluster smaller than `n / n_clusters`. The
  allocation runs over the *selected* clusters, so the remainder goes to
  the largest of them and an undersized cluster is capped with its
  shortfall dealt to whichever clusters came with it. A row's
  probability then depends on the company it keeps. On clusters of
  20/10/8/6 with `n = 7` over 2 clusters, the average take is out by up
  to 17%, and a cluster smaller than the take can push the figure above
  1.

- `design_weighted(method = "successive")` — the default. Successive
  sampling has no closed-form inclusion probability; this is the whole
  reason the other two methods exist. Note that its realised
  probabilities are not merely unknown but genuinely different from
  proportional-to-size: for weights `1:10` at `n = 5`, simulation puts
  them up to 35% away.

- [`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md)
  — resampling with replacement from the sample is not a probability
  sample of a finite population, so there is no `pi` to report.

For all but the bootstrap, pass `simulate = TRUE` to estimate them by
Monte Carlo instead. The estimate carries `R`-sized error, which is fine
for checking a design and not fine for publishing a variance. Simulating
a bootstrap is refused rather than answered: every row turns up in some
replicate, so the count converges to 1 for all of them.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md), which
attaches these to a sample when `weights = TRUE`.

## Examples

``` r
df <- data.frame(id = 1:20, site = rep(c("a", "b"), times = c(15, 5)))

# 8 rows allocated proportionally: 6 of 15, then 2 of 5
inclusion_prob(df, design_stratified("site", n = 8))
#>  [1] 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4
#> [20] 0.4

sampling_weight(df, design_stratified("site", n = 8))
#>  [1] 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5
#> [20] 2.5

# The default weighted design has no closed form; ask for a simulation
w <- data.frame(id = 1:5, w = c(1, 1, 1, 1, 16))
inclusion_prob(w, design_weighted("w", n = 2), simulate = TRUE, R = 2000,
               seed = 1)
#> [1] 0.2785 0.2310 0.2700 0.2465 0.9740

# Systematic PPS does have one, and it is exactly proportional to size
inclusion_prob(w, design_weighted("w", n = 2, method = "systematic"))
#> [1] 0.25 0.25 0.25 0.25 1.00
```
