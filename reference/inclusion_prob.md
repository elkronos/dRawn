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

  Number of simulated draws when `simulate = TRUE`. See "What a
  simulation costs".

- seed:

  Optional seed for the simulation, applied as in
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) (see
  "Seeding" there). Checked whether or not `simulate` is `TRUE`.

## Value

A numeric vector with one element per row of `data`.

A row the design can never select gets `0` from `inclusion_prob()` and
`NA` from `sampling_weight()`, since no finite number of population rows
is represented by a row that cannot be drawn. Such a row lies outside
the region or the time window, comes past `max_items`, or was dropped
under `na_rm = TRUE`.
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
counts them. A weight of 0 is not one of these. The designs that select
in proportion to size refuse it, so drop such rows from `data` first.

A design that
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) refuses,
such as an `n` larger than the frame or than the rows inside a region,
is refused here too, with the same message, rather than given
probabilities above 1 or a census.

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
| [`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md) | `1 / interval` with a random `start`, or `n / N` with `n`, for every row walked, exactly |
| [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md) | `n_clusters / N_clusters`; with `balanced = TRUE`, that times the cluster's average take over its size, worked out exactly (see [`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)) |
| [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md) | `(n_clusters / N_clusters) * (n_h / N_h)`, when the per-cluster take `n_h` is constant; with replacement, `(n_clusters / N_clusters) * (1 - (1 - 1 / N_h)^n_h)`; with `size`, the cluster's `pi_h = n_clusters * M_h / sum(M)` in place of `n_clusters / N_clusters` |
| [`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md) | `n * p_i` for the `"systematic"` and `"poisson"` methods, with rows past 1 taken with certainty and the rest rescaled to sum to `n`; with `strata`, `n_h * p_i` within each stratum |
| [`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md) | `1` above the threshold, `rest`'s own below it |
| [`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md) | `n / min(N, max_items)`, and `0` past `max_items` |
| [`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md) | `per_interval / N_bucket` within each interval |
| [`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md) | `n / N_in_region` |
| [`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md) | `n / N`, or `n * p_i` with a `size` column; with `strata`, the same within each stratum |

A design with `reserve` gives the probabilities of its main sample,
which the reserves do not change. The design of a sample from
[`use_reserves()`](https://elkronos.github.io/dRawn/reference/use_reserves.md)
gives `n / (N - u)`, and `0` for the `u` rows found unavailable.

Four cases have no closed form, and the package refuses to invent one:

- `design_multistage(allocation = "proportional")` — the stage-two
  allocation depends on which clusters were selected. Simulation on
  clusters of 3/5/7/9 gives 0.17, 0.20, 0.17, 0.15, against a naive
  0.33, 0.20, 0.14, 0.11.

- [`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
  where the per-cluster take is not constant — `n` not divisible by
  `n_clusters`, or, without replacement, a cluster smaller than
  `n / n_clusters`. The allocation runs over the *selected* clusters, so
  the remainder goes to the largest of them and an undersized cluster is
  capped with its shortfall dealt to whichever clusters came with it. A
  row's probability then depends on the company it keeps. On clusters of
  20/10/8/6 with `n = 7` over 2 clusters, the average take is out by up
  to 17%, and a cluster smaller than the take can push the figure above
  1.

  For these two, the exact probability is the probability given the
  clusters selected, averaged over every possible selection, and that
  average has no formula. `draw(weights = TRUE)` weights the sample
  conditional on the clusters it selected instead, which gives an
  unbiased total without the average; see "Design weights" in
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

- `design_weighted(method = "successive")` — the default. Successive
  sampling has no closed-form inclusion probability; this is the whole
  reason the other two methods exist. Note that its realised
  probabilities are not merely unknown but genuinely different from
  proportional-to-size: for weights `1:10` at `n = 5`, enumerating every
  ordered draw puts them up to about 36% away.

- [`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md)
  — resampling with replacement from the sample is not a probability
  sample of a finite population, so there is no `pi` to report.

A
[`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md)
with a fixed `start` is refused too, for another reason. It is
deterministic: each row has probability 0 or 1, and `simulate = TRUE`
reports exactly that. Without `simulate` those values are refused,
because `draw(weights = TRUE)` would give every sampled row a weight of
1 and
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
would report the sum over the sample as the total of the frame. So are
`sampling_weight()`,
[`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md),
`draw(weights = TRUE)` and the estimators. Leave `start` as `NULL` for a
probability sample.

For all but the bootstrap, pass `simulate = TRUE` to estimate them by
Monte Carlo instead. The estimate carries `R`-sized error, which is fine
for checking a design and not fine for publishing a variance. Simulating
a bootstrap is refused rather than answered: every row turns up in some
replicate, so the count converges to 1 for all of them.

## What a simulation costs

`simulate = TRUE` makes `R` draws, 5,000 by default, so it takes `R`
times as long as one
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md). Most
designs draw from a frame of 10,000 rows in well under a second. The
local pivotal method of
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
does work in proportion to `N^2` on every draw, so at 10,000 rows each
draw takes seconds and the default `R` takes hours. When the draws after
the first two, which carry one-off costs, project a run of more than 30
seconds, a message says how long it will take, and an interactive
session shows a progress bar. Lower `R` for a quicker answer with more
Monte Carlo error: its standard error is `sqrt(p * (1 - p) / R)`, at
most 0.007 for `R = 5000`.

A warning that every draw raises, such as a stratum allocated no rows or
`max_items` cutting a frame short, is given once at the end, with the
number of draws that raised it.

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
