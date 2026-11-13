# drawn (development version)

## Wrong answers, fixed

* **Reordering, subsetting or binding a weighted sample gave wrong standard
  errors.** The sample matched its rows to the frame by position, so a sort
  moved each row onto another row's stratum, cluster or systematic position.
  Sorting one stratified sample by `y` cut its standard error from 3,015 to
  1,332, and `head()`, `s[cond, ]` and `rbind()` were estimated as if drawn.
  Each row now finds its frame row through the new `.row` column, so
  reordering by `[`, `dplyr::arrange()` or `data.table::setorder()` changes
  nothing, and `ht_total()`, `ht_mean()`, `as_svydesign()` and
  `sample_summary()` refuse a subset, a sample bound to another, or a join
  that repeats rows.
* **Editing the frame in place after drawing changed the sample's
  estimates.** The sample kept a reference to the frame, so a data.table
  `setkey()`, `setorder()` or `:=` on it afterwards changed the standard
  error, interval, degrees of freedom and design effect of a sample already
  drawn. In one example `setkey()` halved a cluster sample's standard error;
  in another `setorder()` took a certainty sample's design effect from 1.04
  to 171. `draw()` now stores a copy of the frame.
* **A certainty design over a with-replacement `rest` was weighted anyway.**
  `draw(weights = TRUE)` checked `replace` only on the outer design, so
  `design_certainty()` over `design_simple(replace = TRUE)` or
  `design_stratified(replace = TRUE)` returned duplicates weighted as distinct
  units, and in one example a total a third too high. It is now refused like
  any other with-replacement design. Drawing it without weights is still
  allowed.
* **Renaming, recoding or dropping a design column in a sample changed its
  standard errors.** The degrees of freedom, the jackknife, `as_svydesign()`
  and `sample_summary()` read strata and cluster labels from the sample's
  columns, while the analytic variance read them from the frame. Dropping a
  cluster sample's cluster column set the degrees of freedom to 0, which
  removed the interval, and made the jackknife delete rows instead of
  clusters, 61% low in one example. Recoding stratum labels broke the
  jackknife and `as_svydesign()`, or gave them the wrong strata. Every path
  now reads design columns from the copy of the frame stored with the sample,
  so only `y` and `by` come from the sample's own columns, and
  `as_svydesign()` passes the cluster labels to `survey` in a `.cluster`
  column.
* **A spatial sample analysed under another `sf::sf_use_s2()` setting got
  the wrong standard error.** Which rows lie inside the region was worked out
  again in the current session. A sample drawn under `sf::sf_use_s2(FALSE)`,
  the documented workaround for a region that crosses the antimeridian, and
  analysed after switching back or after `readRDS()` in a new session got
  standard errors for a different region: 5.7 times too large in one example,
  and 0 with a zero-width interval in another. The point estimate was not
  affected. `ht_total()`, `ht_mean()`, `as_svydesign()` and
  `sample_summary()` now check that the design still gives the frame the
  inclusion probabilities it gave at draw time, and refuse the sample when it
  does not, naming `sf::sf_use_s2()` for a spatial design. `?design_spatial`
  says the setting must be the same at draw time and after.
* **The jackknife and the local-mean variance treated rows taken with
  certainty as random.** A row that is in every possible sample was deleted,
  or compared with its neighbours, as if it could have been left out. With
  the units `design_weighted(method = "systematic")` takes with certainty,
  the jackknife averaged 181 times the true variance on one frame, and
  `variance = "auto"` fell back to it when one sampled row was below
  certainty. `design_spread(size =)` with dominant sizes averaged 1,270 times
  the true variance. Both estimators now hold certain units fixed, and the
  jackknife's finite population correction counts only the other units. With
  one unit below certainty the variance is reported as missing instead.
* **A census got a nonzero jackknife or local-mean standard error.**
  `variance = "jackknife"` on a sample of every row or every cluster applied
  no finite population correction, giving 114 for a six-cluster census in one
  example, and `design_spread()` with `n` equal to the frame reported a
  standard error of about 5% of the exact total. Both now give 0, as the
  analytic estimators did.
* **`design_weighted()` could give a certain unit a probability a rounding
  error below 1.** Depending on the scale of the sizes, a unit whose share came to
  exactly 1 got `1 - 1e-16`, so `.prob == 1` missed it and a census reported
  a small standard error under a "Deville approximation" note. A share within
  `1e-9` of 1 is now taken with certainty at exactly 1, in
  `design_weighted()` and `design_spread(size =)`.
* **A certainty design whose `rest` drew no rows reported an exact total.**
  When a `rest` of random size, such as Poisson sampling with a small `n`,
  drew none of the rows below the threshold, `ht_total()` and `ht_mean()`
  gave a standard error of 0 and a zero-width interval although those rows
  went unobserved. The total is exact only when no row lies below the
  threshold. Otherwise the variance is now missing, with a note naming the
  rows `rest` did not reach.

## Refusals and clearer errors

* `as_svydesign()` and `sample_summary()` give the same reason as
  `ht_total()` when a sample has lost its design, and all three name the
  steps that drop it (`merge()`, column selection, `subset()`, `transform()`,
  `cbind()`) and how to add fieldwork columns instead.
* The refusal of `weights = TRUE` for a with-replacement design, and
  `?draw`, gave the bias as "about 15%". It is a factor of
  `(n / N) / (1 - (1 - 1 / N)^n)`, about 5% at `n = N / 10`, 58% at `n = N`
  and more beyond, applied within each stratum, and both now say so and give
  the per-stratum Hansen-Hurwitz form.
* `ht_total()`, `ht_mean()`, `as_svydesign()` and `sample_summary()` refuse a
  sample whose `.prob` column was removed, is not numeric, holds a value
  outside (0, 1], or differs from the probabilities recorded at draw time. A
  removed `.prob` gave a total of 0, and a missing value gave `NA` with a
  note saying the jackknife had run. The message says that adjusted weights
  belong in `.weight`.
* `ht_total()` and `ht_mean()` refuse an infinite or `NaN` `y`, with the
  count and the rows, instead of returning an infinite total or standard
  error. A variance that overflows is reported as missing, and the note no
  longer says the jackknife ran when it gave no figure.
* `ht_total()` and `ht_mean()` warn, with class `drawn_weight_edited`, when
  `.weight` is no longer `1 / .prob`. They use the design probabilities in
  `.prob`, while `as_svydesign()` hands `.weight` to `survey`, so after a
  nonresponse adjustment the two gave different totals without a word.
  `?ht_total` and `?as_svydesign` say which column each uses.
* `ht_total()` and `ht_mean()` warn, with class `drawn_grouped_df`, when the
  sample is grouped with `dplyr::group_by()`, which they ignore, and give the
  `by` that estimates per group.
* A `by` column that the sample does not have is named in the error.
* `draw(weights = TRUE)` refuses a `design_certainty()` whose `rest` gives
  every row below the threshold a probability of 0, such as one with
  `n = 0`, because the weighted total would leave those rows out. Drawing
  without weights still gives the cut-off sample of the rows above the
  threshold.
* An error raised by a certainty design's `rest` in `draw()` or
  `inclusion_prob()` now says which rows `rest` was given, as in "In `rest`,
  applied to the 230 rows with `w` below 100". `rest` counts only those rows,
  so its own message spoke of 230 rows on a frame of 240.

## Behaviour changes

* `draw(weights = TRUE)` adds a third column, `.row`, after `.prob` and
  `.weight`. It holds the row of `data` each sampled row came from. A frame
  that already has a `.row` column is refused, as one with `.prob` or
  `.weight` is.
* A sample drawn from a data.table takes `:=` without the "Invalid
  .internal.selfref" warning. data.table joins `Suggests`.

## Handing off to survey

* `as_svydesign()` hands `survey` a copy of the sample, so an in-place sort
  of the sample afterwards no longer changes the survey design's standard
  errors.

# drawn 0.2.0

An adversarial review of 0.1.0 found estimates that were wrong without saying
so. Every fix below was reproduced by simulation first, and each now has a
test that fails on the old behaviour.

## Wrong answers, fixed

* **One sampled row per stratum reported a variance of 0.** Sen-Yates-Grundy
  has no within-stratum pair to sum over, so a stratified sample with one row
  in every stratum printed a standard error of exactly zero and a zero-width
  interval. `ht_total()` and `ht_mean()` now decline with the reason, as
  `survey` does by default, for strata and for `design_temporal()` intervals.
* **A stratum allocated no rows silently biased the total.** Proportional
  allocation gives nothing to a stratum smaller than about `N / n`, and Neyman
  allocation gives nothing to one whose auxiliary variable is constant. Its
  rows then have probability 0 and drop out of every estimate; one example was
  84% low. `draw()` now warns, naming the strata and the fix. The
  documentation previously claimed `min_per_stratum = 0` "leaves allocation
  unbiased" and that a floor biases a total; both were wrong, and are
  corrected.
* **The jackknife ignored strata.** It deleted rows across the whole sample,
  counting the differences *between* strata as sampling variance: sixty times
  too large on a well-stratified frame. It is now the stratified (JKn)
  jackknife, and equals the analytic figure exactly for `design_stratified()`
  and `design_temporal()`.
* **One row per cluster understated a multistage variance by ~30%.** The
  exact two-stage estimator cannot see within-cluster variation when only one
  row is taken per cluster. It now declines, and `variance = "auto"` falls back
  to the delete-a-cluster (ultimate-cluster) jackknife.
* **Neyman allocation abandoned its own rule when a stratum was capped.** The
  surplus from a stratum too small for its share was dealt out round-robin;
  it is now re-split in proportion to `size * sd`, as in Cochran (1977, 5.9).
* `ht_mean(estimator = "ht")` printed itself as a Hajek mean.
* `design_systematic(order_by = )` returned rows in sorted order, breaking the
  documented promise that rows come back in frame order.
* The documentation for `design_weighted(method = "systematic")` said some
  pairs of rows can never be drawn together and pointed to
  `sampling::UPsystematicpi2()` for joint probabilities. The rows are shuffled
  before the walk, so every pair can co-occur and that function answers a
  different design.

## Variance estimators from the literature

Every estimate now reports which estimator produced its variance in `method`,
and names any approximation when it prints.

* **Confidence intervals use t on the design's degrees of freedom** — primary
  sampling units minus strata (Korn and Graubard 1999), matching
  `survey::degf()`. At four clusters the normal interval covered about 87%;
  the t interval restores nominal coverage. `df = Inf` gives the normal
  interval.
* **Systematic PPS uses Deville's (1999) approximation**, which needs only
  first-order probabilities and scales each row by its own `1 - pi`. It
  equals `sampling::varest()` exactly and tracks the empirical variance at
  0.98, where the previous jackknife ran at 0.89.
* **Systematic samples get the successive-difference approximation** (Wolter
  2007) rather than no variance at all. On a frame sorted by something related
  to `y` it stays near the truth where the simple random formula overstates
  two- to thirtyfold. Its limits — a cycle matching the interval, or a trend
  so smooth the random start is the only variation — are documented.

## New

* `design_spread()` draws a spatially balanced sample by the local pivotal
  method (Grafström, Lundström and Schelin 2012): exact first-order
  inclusion probabilities, equal or proportional to `size`, with neighbours
  competing so the sample covers the space evenly. On a smooth surface its
  variance was 14–25% of a simple random sample's. It spreads across any
  numeric columns, so it also balances a sample on covariates without
  choosing strata. `ht_total()` uses the local mean variance estimator
  (Stevens and Olsen 2003) for it.
* `ht_total()` and `ht_mean()` take `by =` for domain estimates computed from
  the whole sample rather than a subset of it. Standard errors match
  `survey::svyby()` exactly.
* `ht_mean()` accepts a logical `y`, so a proportion is one call.
* `plot(design, data, type = "map")` places every row by two columns and
  fills in the selected ones — the view that shows whether a sample covers the
  ground.
* Key functions carry `@references` to the methods they implement.

## Handing off to survey

* Multistage designs are expressed as two stages, so `survey`'s standard error
  now matches `ht_total()` exactly; it previously differed by up to 5–10%.
* Systematic PPS, alone or as a certainty design's `rest`, maps to
  `pps = "brewer"` and agrees with `ht_total()` within 0.2%.
* A reservoir capped by `max_items` uses the rows the stream reached as its
  population, not the whole frame.

## Documentation and site

* The README is generated from `README.Rmd`, so every output shown is real.
* New articles: a getting-started walkthrough, how the estimates work (every
  estimator, its reference and its simulation check), an audit-sampling
  walkthrough, and a spatially balanced field survey.
* `graphics` and `utils` are declared in `Imports`.

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

Eleven designs ship: `design_simple()`, `design_stratified()`,
`design_systematic()`, `design_cluster()`, `design_multistage()`,
`design_weighted()`, `design_certainty()`, `design_reservoir()`,
`design_bootstrap()`, `design_temporal()` and `design_spatial()`.

`design_certainty()` composes: rows at or above a threshold are all taken, and
any other design is applied to what remains. That is the standard shape of an
audit sample, and the cheapest way to cut variance when a few units dominate a
total. Certainty rows have inclusion probability exactly 1, so they contribute
their own value with a weight of 1 and add nothing to the variance; inclusion
and joint probabilities compose exactly across the two parts.

## One contract across all of them

* `n` is always the total number of rows drawn, never a per-group figure.
  `design_temporal()` uses `per_interval` precisely because that one *is* per
  group.
* `allocation` always says how a total is split across groups.
* `na_rm` always decides whether missing keys are dropped or raise an error,
  and means the same thing whether you are drawing or asking for probabilities.
* `replace` always rules out `weights = TRUE`: an inclusion probability
  describes distinct units, and a sample holding duplicates cannot be weighted
  by one.
* `draw()` restores the caller's random number stream before returning, so
  sampling inside a simulation does not shift the simulation's own draws.
* What comes back has the same class and the same columns, in the same order,
  as what went in — tibble in, tibble out, with or without weights. Rows come
  back in frame order except where draw order is meaningful: `design_simple()`
  and `design_weighted()` return draw order, and `design_bootstrap()` returns
  replicates behind a leading `.replicate` column.

## Inclusion probabilities and estimation

`inclusion_prob()` gives the first-order inclusion probability of every
population row under a design, without drawing anything — a quick way to check
a plan before committing to it. `draw(weights = TRUE)` attaches those
probabilities and the corresponding design weights to a sample as `.prob` and
`.weight`.

`joint_prob()` gives second-order probabilities, and `ht_total()` and
`ht_mean()` use them to return an estimate with a standard error and confidence
interval:

```r
s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
ht_total(s, "spend")
ht_mean(s, "spend")
```

`ht_mean()` defaults to the Hajek estimator, `sum(y/pi) / sum(1/pi)`, which
divides by the estimated population size rather than the known one;
`estimator = "ht"` divides by the true `N`. The two coincide exactly whenever
the weights of the drawn rows sum to `N` — every fixed-size equal-probability
design, and every stratified design without replacement — and part company when
the sample size is random or the weights vary within a fixed-size sample. Hajek
is usually the steadier of the two, which is why it is the default; `?ht_mean`
sets out when it is not, and why `"ht"` is unbiased only where every row is
reachable.

`deff()` reports the design effect on either result — the design's variance
against simple random sampling of the same size, which reads as an exchange rate
on sample size. Stratifying on a variable that matters can put it far below 1;
clustering on the same variable puts it well above.

The variance estimator matches how each design randomises: Sen-Yates-Grundy for
fixed-size designs, the independent-units form for Poisson sampling, the
cluster-level form for cluster designs — whose row count is random whenever
clusters differ in size — and, for a certainty design, whatever `rest` uses,
since certainty rows are in every possible sample and contribute nothing.

Where a design has inclusion probabilities but no closed-form joint ones,
`ht_total()` falls back to a delete-a-group jackknife and reports which method
it used, or `"none"` where neither could produce a figure. The jackknife is
declined for systematic sampling, which has a single primary sampling unit, and
for Poisson sampling, which has an exact variance already.

Every analytic estimator is checked in the test suite against the empirical
sampling variance of its own estimator, every inclusion and joint probability
against the observed frequency over thousands of draws, and every design that
`survey` can express against `survey::svytotal()`.

## Planning a sample, and checking one

`plan_size()` solves for `n` given the precision wanted, rather than asking for
a guess. It takes a margin of error and a measure of spread, and applies three
corrections: a finite population correction from `N`, an inflation for the
design from `deff`, and an inflation for non-response from `response`. It works
for a mean, a proportion or a total, and says so plainly when the margin asked
for is unreachable by sampling.

```r
plan_size(margin = 5, sd = 40, N = 20000, deff = 2.5, response = 0.7)
```

`sample_summary()` reports what was drawn against what was there: the sampling
fraction, the range and coefficient of variation of the design weights, the
per-stratum or per-cluster counts against the frame, and the number of rows the
design could never have reached.

## Handing off to the survey package

`as_svydesign()` builds a `survey::svydesign()` object from a drawn sample, so
the analysis this package does not do — subpopulation estimates, regression,
calibration, quantiles with proper standard errors — can be done by the package
that does.

The two packages compute variance from different starting points: this one from
the design's joint inclusion probabilities, `survey` from the design's shape.
Each design is therefore expressed in `survey`'s own terms — strata for a
stratified or temporal design, cluster ids for a cluster design,
`survey::poisson_sampling()` for Poisson, and a taken-whole stratum for the
certainty rows. Totals and standard errors then agree to floating point for
every design but two: `survey` uses the ultimate-cluster approximation for
multistage designs, and returns the conservative simple-random figure for
systematic sampling where `ht_total()` declines to return one. `?as_svydesign`
has the table. Compositions with no single `survey` equivalent — a certainty
design over a cluster, multistage or Poisson `rest` — are refused rather than
approximated.

## Seeing a design

`plot(design, data)` shows which rows a design selects — every frame row as a
dot, selected ones filled in — and `type = "probability"` plots inclusion
probability against frame position. Base graphics, no plotting dependency.

## Allocation

`design_stratified()` accepts `allocation = "neyman"` with `allocation_by`,
giving strata shares proportional to `size * sd`. That minimises the variance of
a total for a fixed `n` by putting more rows where the values vary most.

## Saying when it cannot

Not every design has a closed-form inclusion probability, and the ones that
don't report that rather than returning an approximation:

* `design_cluster(balanced = TRUE)` — the per-cluster take is the smallest
  *selected* cluster's size, which is random.
* `design_multistage(allocation = "proportional")` — the second-stage
  allocation depends on which clusters were selected.
* `design_multistage()` where the per-cluster take is not constant — `n` not
  divisible by `n_clusters`, or a cluster smaller than `n / n_clusters`. The
  allocation runs over the *selected* clusters, so a row's probability depends
  on which others were drawn.
* `design_weighted(method = "successive")` — successive sampling has no closed
  form. This is the default, and the reason the other two methods exist.
* `design_bootstrap()` — resampling a sample is not a probability sample of a
  finite population.

`simulate = TRUE` estimates any of them by Monte Carlo, for first-order
probabilities and, in `joint_prob()`, for second-order ones — which is the
general answer where no formula exists, including
`design_weighted(method = "systematic")`, whose joint probabilities depend on
the order units are visited. `design_bootstrap()` is the one exception, and
simulating it is refused: every row appears in some replicate, so the count
converges to 1 for all of them. Systematic sampling has first-order
probabilities but no design-unbiased variance at all, because most pairs of rows
can never co-occur; `ht_total()` says so instead of returning a number.

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

* `sf` is a suggested dependency, needed only by `design_spatial()`, and
  `survey` only by `as_svydesign()`. Nothing else requires either.
* Under spherical geometry a longitude/latitude polygon whose edge spans more
  than 180 degrees is drawn the short way across the antimeridian, so a "whole
  world" rectangle collapses to a narrow strip. `design_spatial()` warns when a
  region has such an edge.
* Reservoir sampling takes a vectorised path for a data frame, whose length is
  already known, and Algorithm L for a genuine stream — a connection or a
  generator function.
* Calendar units in `design_temporal()` step by calendar months and years, not
  by fixed 30.44-day durations.
