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

Ten designs share one contract: `n` always means the total drawn, `allocation`
always means how a total is split, `na_rm` always means the same thing, the
caller's RNG stream is restored on exit, and the result is a data frame with the
input's class and column order. A design can also report its own first- and
second-order inclusion probabilities, so `ht_total()` returns a
Horvitz-Thompson total with a standard error.

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
   moving-block bootstrap, sampling within time intervals, and sampling within
   a spatial region are in this package and not in `sampling`.

3. **Declining to approximate.** Four designs have no closed-form inclusion
   probability, and `inclusion_prob()` raises an informative error naming an
   alternative rather than returning a plausible number. A Monte Carlo estimate
   is available on request via `simulate = TRUE`. Similarly, `ht_total()`
   returns `NA` for the variance under systematic sampling and explains that no
   design-unbiased estimator exists, rather than silently substituting the
   simple random sampling formula.

**`survey`** (Lumley) analyses complex survey data given a design specification;
it does not draw samples. The two are complementary, and a `drawn` sample
carries the `.weight` column that `survey::svydesign()` expects.

**`drawsample`** selects rows so a subsample matches target distributional
characteristics. That is a data-shaping tool rather than probability sampling,
so the overlap is in the name only.

**`sampler`** computes sample sizes and draws simple and stratified samples,
then reports margins of error.

## Naming

The package is named `drawn`, which does not differ only in case from any
existing CRAN or Bioconductor package. `drawr` was avoided because of `DRaWR`,
and `sdraw` because of `SDraw`.

## Correctness checks

Inclusion probabilities are verified against 20,000-draw simulations, and each
variance estimator against the empirical sampling variance of its own estimator
over 4,000 replications. These run as part of the test suite (416 tests).

## First submission

This is a first submission, so there are no reverse dependencies.
