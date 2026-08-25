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
plan <- design_stratified(strata = "site", n = 60)

s <- draw(invoices, plan, seed = 1, weights = TRUE)
ht_total(s, "value")
#> Horvitz-Thompson total  (stratified design, n = 60)
#>   estimate 380,046.5
#>   se       49,499.7  (analytic)
#>   95% CI  283,028.9 to 477,064.1
#>   deff     0.994  (about the same as simple random sampling)
```

That standard error is the point. Without it you have a number; with it you have
an estimate — and the true total, 428,704.5, sits inside that interval.

Every example below uses the same 600-row frame, so you can follow along:

```r
set.seed(42)
invoices <- data.frame(
  id    = 1:600,
  site  = rep(c("north", "south", "east", "west"), times = c(300, 180, 90, 30)),
  team  = rep(paste0("t", 1:20), each = 30),
  value = round(rlnorm(600, meanlog = 6, sdlog = 1.1), 2)
)
```

## Why you might want this

**You need an estimate, not just a subset.** `dplyr::slice_sample()` gives you
rows. It cannot tell you the sampling variance of a total computed from them,
because it doesn't retain what the sampling scheme was. `drawn` does, so
`ht_total()` can produce a standard error and a confidence interval.

**You don't know how big the sample should be.** You rarely do. What you know is
the margin of error you can live with, and `plan_size()` solves for `n` from
that:

```r
plan_size(margin = 0.03, N = 20000, target = "proportion")
#> Sample size for a proportion
#>   draw           1,014
#>   margin         +/- 0.03 at 95% confidence
#>   assuming       sd 0.5, N 20,000
```

**You want to check a plan before committing to it.** Because the design is a
value, you can interrogate it against your population without drawing anything:

```r
tapply(inclusion_prob(invoices, plan), invoices$site, unique)
#>  east north south  west
#>   0.1   0.1   0.1   0.1
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
inclusion_prob(invoices, design_weighted("value", n = 60))
#> Error: `design_weighted(method = "successive")` has no closed-form
#> inclusion probability.
#> Successive sampling has no closed-form inclusion probability.
#> Use method = "systematic" or "poisson" for a design whose inclusion
#> probabilities really are proportional to the weights.
#> Or pass simulate = TRUE to estimate it by Monte Carlo.
```

## Installation

```r
# install.packages("remotes")
remotes::install_github("elkronos/dRawn")
```

## How many rows do you need?

Start here, before choosing a design. `plan_size()` inverts the usual question:
you supply the precision you want and it returns the sample size that buys it.

```r
# A proportion, no prior guess about its value, 20,000 in the frame
plan_size(margin = 0.03, N = 20000, target = "proportion")

# A mean, when a pilot put the spread near 40
plan_size(margin = 5, sd = 40, N = 20000)
#> Sample size for a mean
#>   draw           243
#>   margin         +/- 5 at 95% confidence
#>   assuming       sd 40, N 20,000
```

Three corrections matter and each maps to a real-world fact:

- **`N`** applies a finite population correction. Drawing 400 from 500 is a very
  different proposition from 400 from 500,000, and past a point a bigger frame
  stops mattering at all.
- **`deff`** inflates for the design. A clustered sample of 400 may carry the
  information of 150; pass the `deff()` measured on a comparable past sample
  (see *Estimating* below) and the arithmetic accounts for it.
- **`response`** inflates for non-response, so you draw enough to *end up* with
  what you need.

```r
plan_size(margin = 5, sd = 40, N = 20000, deff = 2.5, response = 0.7)
#> Sample size for a mean
#>   draw           853
#>   to analyse     597  (after 70% response)
#>   margin         +/- 5 at 95% confidence
#>   assuming       sd 40, deff 2.5, N 20,000
```

If the margin you asked for is unreachable by sampling, it says so rather than
returning a number larger than your frame.

## How to use it

Three steps, always the same.

**1. Describe the design.**

```r
library(drawn)

plan <- design_stratified(strata = "site", n = 60)
plan
#> <sampling design: stratified>
#>   strata           "site"
#>   n                60
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

## The eleven designs

| Design | Selects | Size argument |
|---|---|---|
| `design_simple()` | rows uniformly at random | `n` |
| `design_stratified()` | a share of each stratum | `n` (total) |
| `design_systematic()` | every *k*-th row from a random start | `interval` |
| `design_cluster()` | whole clusters | `n_clusters` |
| `design_multistage()` | clusters, then rows within them | `n_clusters` and `n` |
| `design_weighted()` | rows with probability driven by a weight | `n` |
| `design_certainty()` | everything above a threshold, plus a sample of the rest | `rest`'s own |
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
- **`replace` always rules out `weights = TRUE`.** An inclusion probability
  describes distinct units, and a sample holding duplicates cannot be weighted
  by one — `sum(y * .weight)` over it would come out about 15% high. `draw()`
  refuses rather than returning weights that look fine.
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
  out — with or without `weights = TRUE`. Rows come back in frame order, except
  for `design_simple()` and `design_weighted()`, which return them in draw
  order, and `design_bootstrap()`, which returns replicates in order behind a
  leading `.replicate` column.

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

**PPS is a bet, not a free improvement.** It buys precision when the size
measure predicts what you are measuring, and costs precision when it doesn't —
it spends the sample on large units, which only pays off if large units carry
more of the quantity you are totalling. On the same frame, estimating a cost
that scales with invoice size, PPS beat proportional stratification by a factor
of six; estimating a count that does *not* scale with size, it was nearly four
times worse. `deff()` tells you which way the bet went.

## Certainty strata: the units too big to leave to chance

When a handful of units dominate the total — a few enormous invoices, a few
huge stores — leaving them to chance is what makes an estimate wobble. Take them
all, and sample the rest:

```r
d <- design_certainty("value", threshold = 2000,
                      rest = design_stratified("site", n = 60))

s <- draw(invoices, d, seed = 1, weights = TRUE)
table(certain = s$.prob == 1)
#> certain
#> FALSE  TRUE
#>    60    36
```

The 36 rows above the threshold have inclusion probability exactly 1, so they
carry a weight of 1 and add **nothing** to the variance. All the uncertainty in
the estimate comes from the part you actually sampled — which is the whole
point, and why this is the standard shape of an audit or financial sample.

```r
ht_total(s, "value")
#> Horvitz-Thompson total  (certainty design, n = 96)
#>   estimate 458,373.7
#>   se       34,691.5  (analytic)
#>   95% CI  390,379.6 to 526,367.8
#>   deff     0.314  (better than simple random sampling)
```

`rest` takes any design, and its `n` is the number drawn from the rows *below*
the threshold, not the total. Inclusion and joint probabilities compose
correctly across the two parts, so the variance is exact rather than an
approximation.

## Estimating: totals, means, and what the design cost

`ht_total()` estimates a population total; `ht_mean()` estimates a mean. Both
return a standard error, a confidence interval, and a design effect. To show
what a design choice is worth, here is a frame where the grouping variable
genuinely matters — four regions on very different value levels, and delivery
routes that each sit inside one region:

```r
set.seed(11)
pop <- data.frame(
  id     = 1:400,
  region = rep(c("north", "south", "east", "west"), each = 100),
  route  = rep(paste0("r", 1:40), each = 10)
)
pop$value <- rep(c(120, 260, 480, 900), each = 100) + round(rnorm(400, 0, 40))

by_region <- draw(pop, design_stratified("region", n = 40), seed = 1, weights = TRUE)
by_route  <- draw(pop, design_cluster("route", n_clusters = 4), seed = 1, weights = TRUE)

ht_mean(by_region, "value")
#> Hajek mean  (stratified design, n = 40)
#>   estimate 436.775
#>   se       6.478315  (analytic)
#>   95% CI  424.0777 to 449.4723
#>   deff     0.0221  (better than simple random sampling)

ht_mean(by_route, "value")
#> Hajek mean  (cluster design, n = 40)
#>   estimate 503.8
#>   se       216.1592  (analytic)
#>   95% CI  80.13575 to 927.4643
#>   deff     12.9  (worse than simple random sampling)
```

Same population, same 40 rows, standard errors a factor of thirty-three apart.

**`deff()` is the exchange rate on sample size.** It compares the design's
variance against simple random sampling of the same size: `deff = 2` means a
sample of 400 carries about as much information as 200 drawn at random. Above 1
is the usual price of clustering; below 1 is what stratification and
size-proportional selection buy you.

```r
deff(ht_mean(by_region, "value"))   #> 0.0221
deff(ht_mean(by_route,  "value"))   #> 12.9
```

That number is exactly what `plan_size(deff = )` wants for the next study, which
closes the loop: measure the design effect once, size the next sample honestly.

**Which mean?** `ht_mean()` defaults to the **Hájek** estimator,
`sum(y/pi) / sum(1/pi)`, which divides by the *estimated* population size rather
than the known one; `estimator = "ht"` divides by the true `N`.

The two coincide *exactly* whenever the weights of the rows you drew sum to `N`
— which is every fixed-size equal-probability design, and every stratified
design without replacement. There is nothing to choose between them there. They
part company when the sample size is random (Poisson, clusters of unequal size)
or the weights vary within a fixed-size sample (PPS). Hájek is usually the
steadier of the two, which is why it is the default; the exception is when `y`
is close to proportional to the size measure that drove selection, where `y/pi`
is nearly constant and dividing by the known `N` wins. One caveat on `"ht"`: it
is unbiased for the frame mean only if every row could have been selected —
rows at probability 0 sit in the `N` it divides by but can never enter the
numerator. `sample_summary()` counts them.

## Checking what you got

`sample_summary()` reports what was drawn against what was there. Worth a look
before analysing:

```r
sample_summary(draw(invoices, design_stratified("site", n = 60), seed = 1,
                    weights = TRUE))
#> Sample of 60 from 600  (stratified design)
#>   sampling fraction  0.1
#>   design weights     10 to 10   (cv 0)
#>
#>   by site:
#>     group  drawn  in frame   rate
#>     east       9       90  0.100
#>     north     30      300  0.100
#>     south     18      180  0.100
#>     west       3       30  0.100
```

It surfaces the three things that quietly ruin an estimate: strata that came up
short, weights that vary far more than you expected, and rows the design could
never have reached at all.

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
second-order ones, and `ht_total()` and `ht_mean()` need both to produce a
standard error.

| Design | Inclusion probability | Variance |
|---|---|---|
| simple, stratified, cluster, reservoir, temporal, spatial | exact | yes |
| `design_certainty()` | exact | whatever `rest` supports |
| multistage, equal allocation, constant per-cluster take | exact | yes |
| `design_weighted(method = "poisson")` | exact | yes |
| `design_weighted(method = "systematic")` | exact | jackknife |
| `design_systematic()` | exact | none — most pairs can never co-occur |
| `design_cluster(balanced = TRUE)` | none | no |
| `design_multistage(allocation = "proportional")` | none | no |
| multistage where the per-cluster take varies | none | no |
| `design_weighted(method = "successive")` | none | no |
| `design_bootstrap()` | none | no |

"Constant per-cluster take" means `n` divides by `n_clusters` and no cluster is
smaller than `n / n_clusters`. Otherwise the allocation runs over the *selected*
clusters — the remainder goes to the largest of them, an undersized one is
capped and its shortfall dealt to whichever clusters came with it — so a row's
probability depends on which other clusters were drawn. There is no closed form
for that, and the package says so instead of averaging.

Where there is no closed form, `simulate = TRUE` estimates it by Monte Carlo —
for first-order probabilities and for joint ones alike:

```r
inclusion_prob(data, design_cluster("site", n_clusters = 4, balanced = TRUE),
               simulate = TRUE, R = 2000, seed = 1)

# Second-order too, for the rows you drew. This is the general answer where no
# formula exists -- slower and noisier, but available for every probability
# design.
joint_prob(data, design_weighted("value", n = 60, method = "systematic"),
           rows = drawn_rows, simulate = TRUE, R = 5000, seed = 1)
```

`design_bootstrap()` is the exception, and simulation is refused for it rather
than answered. Every row turns up in some replicate, so the count converges to 1
for all of them — a confident-looking number that means nothing.

The variance estimator matches how the design actually randomises. Fixed-size
designs get Sen-Yates-Grundy. Poisson sampling, whose size is random and whose
rows are independent, gets the independent-units form. **Cluster designs are
estimated at the cluster level**, because taking whole clusters makes the row
count random whenever clusters differ in size — a row-level formula understates
the variance by more than half on a frame whose clusters run from 2 rows to 10.
Certainty designs hand the problem to `rest`, since rows taken with certainty
are in every possible sample and contribute nothing.

Where no analytic form exists but inclusion probabilities do, `ht_total()` falls
back to a delete-a-group jackknife and says so, reporting which method it used —
and reporting `"none"` when neither could produce a figure, rather than naming a
method that declined.

Systematic sampling gets neither: it has a single primary sampling unit — the
random start — so deleting rows misrepresents the design, and `ht_total()`
declines rather than returning a misleading number. Poisson sampling is refused
the jackknife for the opposite reason: it has an exact variance already, and
deleting rows from a sample whose size is itself random understates it
threefold.

Every analytic estimator is checked in the test suite against the empirical
sampling variance of its own estimator over thousands of replications, and every
inclusion and joint probability against the observed frequency over thousands of
draws. Where `survey` models the same design, the two packages' standard errors
are compared and must agree to floating point.

One caveat worth knowing: the interval is a normal approximation, so it
undercovers when the number of *sampling units* is small. That bites hardest on
cluster designs, where the unit is the cluster rather than the row. Drawing 8
clusters from 24, a nominal 95% interval covers about 92%; at 3 or 4 clusters it
drops into the high 70s and low 80s. The variance itself is right — measured
against the empirical sampling variance it sits within a few percent of 1 — but
treat the interval as indicative and prefer more clusters over more rows within
them.

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
data once you have it — subpopulation estimates, regression, calibration,
quantiles with proper standard errors. `as_svydesign()` hands a sample straight
over.

The two packages compute variance from different starting points: this one from
the design's joint inclusion probabilities, `survey` from the design's *shape*.
So `as_svydesign()` expresses each design in `survey`'s own terms — strata for a
stratified or temporal design, cluster ids for a cluster design,
`survey::poisson_sampling()` for Poisson, and a taken-whole stratum for the
certainty rows — rather than handing over a weight column and hoping.

```r
library(survey)

des <- as_svydesign(s)
svytotal(~value, des)
#>        total    SE
#> value 380046 49500

# The same standard error this package reports, to floating point
c(drawn = ht_total(s, "value")$se, survey = as.numeric(SE(svytotal(~value, des))))
#>   drawn  survey
#> 49499.7 49499.7

# Now the analysis this package does not do
svyby(~value, ~site, des, svymean)
svyglm(value ~ site, des)
```

Standard errors match exactly for simple, stratified, temporal, cluster,
reservoir, spatial, both PPS methods and certainty designs. Two cases differ,
and the help page says why: `survey` uses the ultimate-cluster approximation for
multistage designs, giving a figure around 5–10% smaller; and for systematic
sampling it returns the conservative simple-random figure where `ht_total()`
declines to return anything at all. Compositions with no single `survey`
equivalent are refused rather than approximated.

One thing to watch: `survey` exports its own `deff()`, so `library(survey)`
masks this package's. Call `drawn::deff()` afterwards.

Draw here, analyse there. This package's job ends where `survey`'s begins.

## Learn more

```r
vignette("choosing-a-design", package = "drawn")
```

## Credit

The sampling routines began as an R port of
[sample_py](https://github.com/elkronos/sample_py).

## License

GPL-3
