# drawn 0.1.0

First release.

## Designs are objects

Sampling is expressed as a design you build and then apply, rather than a
function call with the design buried in its arguments:

```r
plan <- design_stratified(strata = "site", n = 500)
draw(data, plan, seed = 1)
```

A design can be printed, stored, reused across data sets, and — the part that
matters for estimation — asked what it does before you draw anything.

Ten designs ship: `design_simple()`, `design_stratified()`,
`design_systematic()`, `design_cluster()`, `design_multistage()`,
`design_weighted()`, `design_reservoir()`, `design_bootstrap()`,
`design_temporal()` and `design_spatial()`.

## One contract across all of them

* `n` is always the total number of rows drawn, never a per-group figure.
  `design_temporal()` uses `per_interval` precisely because that one *is* per
  group.
* `allocation` always says how a total is split across groups.
* `na_rm` always decides whether missing keys are dropped or raise an error.
* `draw()` restores the caller's random number stream before returning, so
  sampling inside a simulation does not shift the simulation's own draws.
* What comes back has the same class and the same columns, in the same order,
  as what went in.

## Inclusion probabilities and estimation

`inclusion_prob()` gives the first-order inclusion probability of every
population row under a design, without drawing anything — a quick way to check
a plan before committing to it. `draw(weights = TRUE)` attaches those
probabilities and the corresponding design weights to a sample as `.prob` and
`.weight`.

`joint_prob()` gives second-order probabilities, and `ht_total()` uses them to
return a Horvitz-Thompson total with a standard error and confidence interval:

```r
s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
ht_total(s, "spend")
```

Variances use the Sen-Yates-Grundy estimator for fixed-size designs and the
independent-units form for Poisson sampling. Each was checked against the
empirical sampling variance of its own estimator over 4,000 replications.

## Saying when it cannot

Not every design has a closed-form inclusion probability, and the ones that
don't report that rather than returning an approximation:

* `design_cluster(balanced = TRUE)` — the per-cluster take is the smallest
  *selected* cluster's size, which is random.
* `design_multistage(allocation = "proportional")` — the second-stage
  allocation depends on which clusters were selected.
* `design_weighted(method = "successive")` — successive sampling has no closed
  form. This is the default, and the reason the other two methods exist.
* `design_bootstrap()` — resampling a sample is not a probability sample of a
  finite population.

`simulate = TRUE` estimates any of them by Monte Carlo. Systematic sampling has
first-order probabilities but no design-unbiased variance, because most pairs of
rows can never co-occur; `ht_total()` says so instead of returning a number.

## Probability proportional to size

`design_weighted()` offers three methods, and they are not interchangeable:

| `method` | Inclusion probabilities | Sample size |
|---|---|---|
| `"successive"` (default) | Not proportional to weight, no closed form | Fixed |
| `"systematic"` | Exactly `n * p_i` | Fixed |
| `"poisson"` | Exactly `n * p_i` | Random, mean `n` |

The default is what `base::sample(prob = )` does: a good way to bias selection
toward heavy units, but the weights govern each sequential draw rather than the
probability of ending up in the sample. Use `"systematic"` or `"poisson"` when
the sample will be estimated from. Units heavy enough that `n * p_i > 1` are
taken with certainty and the remainder rescaled, iterating until every
probability is valid.

## Notes

* `sf` is a suggested dependency, needed only by `design_spatial()`. The other
  nine designs have no geospatial requirement.
* Under spherical geometry a longitude/latitude polygon whose edge spans more
  than 180 degrees is drawn the short way across the antimeridian, so a "whole
  world" rectangle collapses to a narrow strip. `design_spatial()` warns when a
  region has such an edge.
* Reservoir sampling takes a vectorised path for a data frame, whose length is
  already known, and Algorithm L for a genuine stream — a connection or a
  generator function.
* Calendar units in `design_temporal()` step by calendar months and years, not
  by fixed 30.44-day durations.
