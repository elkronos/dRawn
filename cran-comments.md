# cran-comments

## Test environments

* local: Ubuntu 24.04, R 4.3.3
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

Twelve designs share one contract: `n` always means the total drawn,
`allocation` always means how a total is split, `na_rm` always means the same
thing on both the draw and the probability path, `replace` always rules out
design weights, the caller's RNG stream is restored on exit, and the result is a
data frame with the input's class and column order. A design can also report its
own first- and second-order inclusion probabilities, so `ht_total()` and
`ht_mean()` return an estimate — overall or by domain — with a standard error
from the variance estimator the literature recommends for that design, a
confidence interval on the design's degrees of freedom, and a design effect. `plan_size()` solves for the sample size a target margin of
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
   moving-block bootstrap, sampling within time intervals, and a composite
   certainty-plus-sample design over any other design are in this package and
   not in `sampling`. `design_spread()` implements the local pivotal method
   (Grafström, Lundström and Schelin 2012), which `BalancedSampling` also
   provides as a vector function; here it is integrated with exact inclusion
   probabilities and a local-mean variance estimator in the same
   draw-and-estimate workflow as every other design.

3. **A variance estimator matched to each design, and named.**
   Sen-Yates-Grundy where joint probabilities are known, Deville's
   approximation for systematic PPS, successive differences for systematic
   samples, local means for spatially balanced ones, and a stratified (JKn)
   jackknife otherwise; t intervals on the design's degrees of freedom.

4. **Declining to mislead.** Five situations have no closed-form inclusion
   probability, and `inclusion_prob()` raises an informative error naming an
   alternative rather than returning a plausible number. A Monte Carlo estimate
   is available on request via `simulate = TRUE`, except for the bootstrap,
   where it would converge to 1 for every row and is refused. `ht_total()`
   declines a variance where a stratum holds a single sampled row, rather than
   report the understated figure the Sen-Yates-Grundy sum gives, and `draw()`
   warns when a stratum is allocated no rows and so could never be sampled.

**`survey`** (Lumley) analyses complex survey data given a design specification;
it does not draw samples. The two are complementary, and `as_svydesign()`
expresses a `drawn` design in `survey`'s own terms — strata, cluster ids, two
stages, `survey::poisson_sampling()`, Brewer's PPS approximation, or a
taken-whole certainty stratum — so that the two packages' totals, standard
errors, domain estimates, degrees of freedom and proportion intervals agree to
floating point wherever `survey` models the same design. The tests assert that
agreement design by design.

**`spsurvey`** implements GRTS spatially balanced designs and an environmental
monitoring workflow; it is the deeper tool for that field.

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

The test suite is built around empirical verification rather than fixed
expected values, because the failure mode this package most needs to avoid is a
plausible-looking wrong number:

* every inclusion probability and joint inclusion probability is checked against
  the observed frequency over thousands of simulated draws;
* every variance estimator is checked against the empirical sampling variance of
  its own estimator, and confidence interval coverage against its nominal level;
* `plan_size()` is checked by drawing the size it recommends and confirming the
  margin it promised is achieved;
* every design `survey` can express is checked against `survey::svytotal()`,
  `survey::svyby()`, `survey::degf()` and `survey::svyciprop()`, and Deville's
  approximation against `sampling::varest()`.

Two guards run in CI alongside `R CMD check`: `tools/check-docs.R` fails the
build on stale or unresolvable documentation, and `tools/check-duplicates.R`
fails it on a function defined in more than one file.

## Known differences from `survey`

Documented in `?as_svydesign` and asserted in the tests:

* systematic PPS — neither package has its joint probabilities; `ht_total()`
  uses Deville's approximation and `survey` Brewer's, and they agree within
  about 0.2%;
* `design_spread()` — `survey` has no estimator for a spatially balanced sample
  and is given the first-order design, so its standard error is larger;
* systematic sampling — `survey`, told only `ids = ~1`, returns the
  simple-random variance; `ht_total()` uses the successive-difference
  approximation. Neither is design-unbiased, because no such estimator exists.

## First submission

This is a first submission, so there are no reverse dependencies.
