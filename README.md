# drawn

<!-- badges: start -->
[![R-CMD-check](https://github.com/elkronos/dRawn/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/elkronos/dRawn/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Take a sample, then say something defensible about the population it came
from.**

Pulling rows out of a data frame is easy. The hard part is what comes next: your
boss asks what the total was across all 40,000 invoices, and you sampled 200.
Answering that honestly needs to know how likely each sampled row was to be
picked — and that number lives in the *design*, not in the sample.

`drawn` keeps the design around so it can answer:

```r
plan <- design_stratified(strata = "site", n = 200)

s <- draw(invoices, plan, seed = 1, weights = TRUE)
ht_total(s, "value")
#> Horvitz-Thompson total  (stratified design, n = 200)
#>   total    1,284,300
#>   se       47,910
#>   95% CI  1,190,398 to 1,378,202
```

That standard error is the point. Without it you have a number; with it you have
an estimate.

## Why you might want this

**You need an estimate, not just a subset.** `dplyr::slice_sample()` gives you
rows. It cannot tell you the sampling variance of a total computed from them,
because it doesn't retain what the sampling scheme was. `drawn` does, so
`ht_total()` can produce a standard error and a confidence interval.

**You want to check a plan before committing to it.** Because the design is a
value, you can interrogate it against your population without drawing anything:

```r
inclusion_prob(invoices, plan) |> tapply(invoices$site, unique)
#>  east north south  west
#>  0.10  0.10  0.10  0.10
```

If a site was going to be sampled at 40% and you expected 10%, you find out now
rather than after the fieldwork.

**You want the same plan applied consistently.** A design is reusable — across
months, across teams, across data sets — and prints as a description of itself,
which makes it reviewable in a way a call buried in a script is not.

**You want to be told when the answer isn't available.** Some designs have no
closed-form inclusion probability. `drawn` says so rather than returning a
plausible-looking number:

```r
inclusion_prob(invoices, design_weighted("value", n = 200))
#> Error: `design_weighted(method = "successive")` has no closed-form
#> inclusion probability.
#> Use method = "systematic" or "poisson" for a design whose inclusion
#> probabilities really are proportional to the weights.
```

## Installation

```r
# install.packages("remotes")
remotes::install_github("elkronos/dRawn")
```

## How to use it

Three steps, always the same.

**1. Describe the design.**

```r
library(drawn)

plan <- design_stratified(strata = "site", n = 200)
plan
#> <sampling design: stratified>
#>   strata           "site"
#>   n                200
#>   allocation       "proportional"
#>   min_per_stratum  0
#>   replace          FALSE
#>   na_rm            FALSE
```

**2. Draw from it.** Ask for `weights = TRUE` if you intend to estimate.

```r
s <- draw(invoices, plan, seed = 1, weights = TRUE)
```

You get back a data frame with the same class and the same columns, in the same
order, as the one you passed in — plus `.prob` (the chance that row was
included) and `.weight` (its reciprocal: how many population rows it stands for).

**3. Estimate.**

```r
ht_total(s, "value")
```

## The ten designs

| Design | Selects | Size argument |
|---|---|---|
| `design_simple()` | rows uniformly at random | `n` |
| `design_stratified()` | a share of each stratum | `n` (total) |
| `design_systematic()` | every *k*-th row from a random start | `interval` |
| `design_cluster()` | whole clusters | `n_clusters` |
| `design_multistage()` | clusters, then rows within them | `n_clusters` and `n` |
| `design_weighted()` | rows with probability driven by a weight | `n` |
| `design_reservoir()` | a fixed-size sample from a stream, in one pass | `n` |
| `design_bootstrap()` | resampled replicates | `n_replicates`, `n` |
| `design_temporal()` | a share of each time interval | `per_interval` |
| `design_spatial()` | rows inside a region | `n` |

```r
draw(data, design_simple(n = 100), seed = 1)
draw(data, design_cluster("site", n_clusters = 5), seed = 1)
draw(data, design_multistage("site", n_clusters = 5, n = 100), seed = 1)
draw(data, design_temporal("when", from = "2024-01-01", to = "2024-01-15",
                           interval = 6, per_interval = 2, unit = "hours"),
     seed = 1)
```

Streams are not data frames, and reservoir sampling treats them differently: a
data frame takes a direct vectorised path, while a connection or a
zero-argument generator function runs Algorithm L in a single pass.

```r
draw(gen, design_reservoir(n = 1000), seed = 1)   # gen() returns NULL when done
```

Bootstrap replicates come back in one frame with a leading `.replicate` column:

```r
reps <- draw(data, design_bootstrap(n_replicates = 500), seed = 1)
vapply(split(reps, reps$.replicate), function(r) mean(r$value), numeric(1))
```

## One contract across all of them

Arguments mean the same thing everywhere:

- **`n` is always the total drawn**, never a per-group figure. `design_temporal()`
  says `per_interval` precisely because that one *is* per group.
- **`allocation`** always says how a total is split across groups:
  `"proportional"`, `"equal"`, or `"neyman"`, which puts more rows where the
  values vary most and minimises the variance of a total for a fixed `n`.

  ```r
  design_stratified("site", n = 60, allocation = "neyman", allocation_by = "value")
  ```
- **`na_rm`** always decides whether missing keys are dropped or raise an error.
  It is never silently assumed.
- **`seed`** is local to the draw. `.Random.seed` is restored on exit, so
  sampling inside a simulation does not shift the simulation's own stream.
- **The result keeps the input's class and column order.** A tibble in, a tibble
  out.

Allocation is exact. Splitting 60 across strata of 300/180/90/30 returns 60
rows, not 61 — the largest-remainder method is used rather than independent
per-stratum rounding.

## Weighted sampling: choose the method deliberately

This is the one place where the obvious call is probably not the one you want.

```r
design_weighted("size", n = 100)                          # successive (default)
design_weighted("size", n = 100, method = "systematic")   # exact piPS, fixed n
design_weighted("size", n = 100, method = "poisson")      # exact piPS, random n
```

| `method` | Inclusion probabilities | Sample size | Estimable |
|---|---|---|---|
| `"successive"` | not proportional to weight | fixed | no closed form |
| `"systematic"` | exactly `n * p_i` | fixed | yes, jackknife variance |
| `"poisson"` | exactly `n * p_i` | random, mean `n` | yes, analytic variance |

The default is what `base::sample(prob = )` does. It biases selection toward
heavy units, which is often all you want — but the weights govern each
sequential *draw*, not the probability of ending up in the sample, so the result
is not probability-proportional-to-size and has no closed-form inclusion
probability. Use `"systematic"` or `"poisson"` if the sample will be estimated
from.

Units heavy enough that `n * p_i > 1` are taken with certainty and the remainder
rescaled, repeatedly, until every probability is valid.

## Seeing what a design does

```r
par(mfrow = c(2, 1), mar = c(2, 1, 2, 1))
plot(design_simple(n = 60), invoices, seed = 1)
plot(design_systematic(interval = 7), invoices, seed = 1)
```

Every frame row is a dot, in frame order, with the selected ones filled in.
Designs look distinct: simple random sampling scatters, systematic makes a
lattice, cluster sampling takes solid contiguous runs, and size-proportional
selection thickens wherever the weight is large. If your frame is sorted by
something meaningful, an unintended pattern shows up straight away.

The second view plots inclusion probability against frame position — flat means
everyone had the same chance, steps mean strata, a slope means
size-proportional, and anything at zero is a row the design can never reach:

```r
plot(design_stratified("site", n = 60), invoices, type = "probability")
```

Base graphics, so there is no plotting dependency to install.

## What can and cannot be estimated

`inclusion_prob()` gives first-order probabilities, `joint_prob()` gives
second-order ones, and `ht_total()` needs both to produce a standard error.

| Design | Inclusion probability | Variance |
|---|---|---|
| simple, stratified, cluster, reservoir, temporal, spatial | exact | yes |
| multistage, `allocation = "equal"` | exact | yes, when `n` divides by `n_clusters` |
| `design_weighted(method = "poisson")` | exact | yes |
| `design_weighted(method = "systematic")` | exact | jackknife |
| multistage, `n` not divisible by `n_clusters` | exact | jackknife |
| `design_systematic()` | exact | none — most pairs can never co-occur |
| `design_cluster(balanced = TRUE)` | none | no |
| `design_multistage(allocation = "proportional")` | none | no |
| `design_weighted(method = "successive")` | none | no |
| `design_bootstrap()` | none | no |

Where there is no closed form, `simulate = TRUE` estimates it by Monte Carlo:

```r
inclusion_prob(data, design_cluster("site", n_clusters = 4, balanced = TRUE),
               simulate = TRUE, R = 2000, seed = 1)
```

Variances use the Sen-Yates-Grundy estimator for fixed-size designs and the
independent-units form for Poisson sampling. Where no analytic form exists but
inclusion probabilities do, `ht_total()` falls back to a delete-a-group
jackknife and says so; it reports which method it used. Every analytic estimator
was checked against the empirical sampling variance of its own estimator over
4,000 replications, and the jackknife reproduces the analytic value exactly for
simple and cluster designs.

Systematic sampling gets neither: it has a single primary sampling unit — the
random start — so deleting rows misrepresents the design, and `ht_total()`
declines rather than returning a misleading number.

One caveat worth knowing: cluster designs with few clusters have few effective
degrees of freedom, so the normal-approximation interval undercovers. With 8
clusters drawn from 24, observed coverage of a nominal 95% interval is about
89%. Treat those intervals as indicative.

## Spatial sampling and the antimeridian

`design_spatial()` takes `coords = c(x, y)` — **longitude first**, matching
`sf::st_as_sf()`.

Under spherical geometry, consecutive polygon vertices are joined by the
*shortest* great-circle path. An edge from longitude −179 to +179 therefore
spans the 2 degrees across the antimeridian, not the 358 the coordinates
suggest, so a "whole world" rectangle collapses to a narrow pole-to-pole strip
of about 2.8 million km² against the globe's 510 million. `drawn` warns when a
region has an edge spanning more than 180 degrees. Split the region at the
antimeridian, or use `sf::sf_use_s2(FALSE)`.

## See also

[`sampling`](https://cran.r-project.org/package=sampling) is the deeper library
for classical design-based sampling: a dozen unequal-probability algorithms,
joint inclusion probabilities for several of them, calibration, and balanced
sampling via the cube method. Reach for it when you need Brewer, Midzuno,
Sampford, Tillé, pivotal or maximum-entropy sampling.

[`survey`](https://cran.r-project.org/package=survey) analyses complex survey
data once you have it. A `drawn` sample carries the `.weight` column that
`svydesign()` expects.

## Learn more

```r
vignette("choosing-a-design", package = "drawn")
```

## Credit

The sampling routines began as an R port of
[sample_py](https://github.com/elkronos/sample_py).

## License

GPL-3
