# cran-comments

## Test environments

* local: macOS (aarch64), R 4.6.1
* GitHub Actions: macOS, Windows and Ubuntu, on R-release, R-devel and oldrel-1

## R CMD check results

0 errors | 0 warnings | 0 notes

## What this package is for

`drawn` expresses a sampling design as a reusable object and applies it with a
single verb:

```r
plan <- design_stratified(strata = "site", n = 500)
draw(data, plan, seed = 1, weights = TRUE)
```

Eleven designs share one contract: `n` always means the total drawn,
`allocation` always means how a total is split, `na_rm` always means the same
thing on both the draw and the probability path, `replace` always rules out
design weights, the caller's RNG stream is restored on exit, and the result is a
data frame with the input's class and column order. A design can also report its
own first- and second-order inclusion probabilities, so `ht_total()` and
`ht_mean()` return an estimate with a standard error, a confidence interval and
a design effect. `plan_size()` solves for the sample size a target margin of
error needs, and `as_svydesign()` hands a sample to `survey` for the analysis
this package does not do.

## Relationship to existing packages

This overlaps a well-established area, and I want to be straightforward about
where it sits.

**`sampling`** (Tillé and Matei) is the reference library for design-based
sampling in R and is considerably deeper than this package on the classical
core. It implements a dozen unequal-probability algorithms where `drawn` has
three, provides joint inclusion probabilities for several of them, and adds
calibration, balanced sampling via the cube method, and ratio and regression
estimators. Anyone who needs Brewer, Midzuno, Sampford, Tillé, pivotal or
maximum-entropy sampling should use `sampling`, and the documentation for
`design_weighted()` and `joint_prob()` says so by name.

`drawn` differs in three ways rather than in statistical novelty:

1. **Interface.** `sampling` is a library of functions over vectors that return
   index tables, which you then join back to your data. `drawn` takes a data
   frame and returns a data frame, and the design is a value you can hold,
   print and reuse.

2. **Scope beyond classical survey designs.** Reservoir sampling from a stream,
   moving-block bootstrap, sampling within time intervals, sampling within a
   spatial region, and a composite certainty-plus-sample design are in this
   package and not in `sampling`.

3. **Declining to approximate.** Five situations have no closed-form inclusion
   probability, and `inclusion_prob()` raises an informative error naming an
   alternative rather than returning a plausible number. A Monte Carlo estimate
   is available on request via `simulate = TRUE`, except for the bootstrap,
   where it would converge to 1 for every row and is refused. Similarly,
   `ht_total()` returns `NA` for the variance under systematic sampling and
   explains that no design-unbiased estimator exists, rather than silently
   substituting the simple random sampling formula.

**`survey`** (Lumley) analyses complex survey data given a design specification;
it does not draw samples. The two are complementary, and `as_svydesign()`
expresses a `drawn` design in `survey`'s own terms — strata, cluster ids,
`survey::poisson_sampling()`, or a taken-whole certainty stratum — so that the
two packages' totals and standard errors agree to floating point wherever
`survey` models the same design. The tests assert that agreement design by
design.

**`drawsample`** selects rows so a subsample matches target distributional
characteristics. That is a data-shaping tool rather than probability sampling,
so the overlap is in the name only.

**`sampler`** computes sample sizes and draws simple and stratified samples,
then reports margins of error. `plan_size()` covers similar ground for sizing,
with finite population, design effect and non-response corrections.

## Naming

The package is named `drawn`, which does not differ only in case from any
existing CRAN or Bioconductor package. `drawr` was avoided because of `DRaWR`,
and `sdraw` because of `SDraw`.

## Correctness checks

The test suite (875 tests) is built around empirical verification rather than
fixed expected values, because the failure mode this package most needs to avoid
is a plausible-looking wrong number:

* every inclusion probability and joint inclusion probability is checked against
  the observed frequency over thousands of simulated draws;
* every variance estimator is checked against the empirical sampling variance of
  its own estimator, and confidence interval coverage against its nominal level;
* `plan_size()` is checked by drawing the size it recommends and confirming the
  margin it promised is achieved;
* every design `survey` can express is checked against `survey::svytotal()` for
  both the total and its standard error.

Two guards run in CI alongside `R CMD check`: `tools/check-docs.R` fails the
build on stale or unresolvable documentation, and `tools/check-duplicates.R`
fails it on a function defined in more than one file.

## Known differences from `survey`

Documented in `?as_svydesign` and asserted in the tests:

* multistage designs — `survey` uses the ultimate-cluster approximation, which
  ignores second-stage sampling; `ht_total()` uses the exact two-stage form and
  is around 5–10% larger;
* systematic sampling — `survey`, told only `ids = ~1`, returns the
  simple-random variance; `ht_total()` declines, because no design-unbiased
  estimator exists.

## First submission

This is a first submission, so there are no reverse dependencies.
