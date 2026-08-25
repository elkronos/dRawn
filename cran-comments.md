# cran-comments

## Test environments

* local: Ubuntu 24.04, R 4.3.3
* (fill in win-builder / R-hub / GitHub Actions results before submitting)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes for the reviewer

### Relationship to existing sampling packages

CRAN already hosts several packages with "sample" or "draw" in the name. The
closest by name is **drawsample**, and the closest by domain are **sampler** and
**sampling**. `drawn` differs from each in what it is for:

* **drawsample** selects rows so that the resulting subsample matches target
  distributional characteristics (skewness, kurtosis, and so on). It is a
  data-shaping tool. `drawn` implements design-based probability sampling: the
  selection mechanism is specified in advance and every unit has a known
  inclusion probability.
* **sampler** computes sample sizes and draws simple and stratified samples,
  then reports margins of error. `drawn` covers ten designs rather than two, and
  its distinguishing feature is that a design object can report its own
  first-order inclusion probabilities via `inclusion_prob()`, which `draw()`
  attaches to a sample as `.prob` and `.weight` so that Horvitz-Thompson
  estimation is possible.
* **sampling** provides a large library of selection algorithms as free
  functions. `drawn` wraps a smaller set behind a single reusable design object
  and one verb, and is deliberately explicit about which designs have no
  closed-form inclusion probability.

Where a design has no closed form — balanced cluster sampling, proportional
multi-stage allocation, successive weighted sampling, and the bootstrap —
`inclusion_prob()` raises an informative error naming an alternative rather than
returning an approximation. A Monte Carlo estimate is available on request via
`simulate = TRUE`.

### First submission

This is a first submission, so there are no reverse dependencies.

The package name deliberately avoids `sampleR`, which differs from the existing
CRAN package `sampler` only in the case of its final letter, and `drawr`, which
would collide with `DRaWR` the same way.
