# Choosing a sampling design

``` r

library(drawn)
```

Every function in this package answers one question — *which rows do I
look at?* — and the interesting part is what you can say afterwards
about the rows you didn’t.

## One population, eleven questions

We’ll use a synthetic audit population throughout: invoices from a
handful of sites, over a fortnight, with a value attached.

``` r

set.seed(42)
n <- 600
invoices <- data.frame(
  id    = seq_len(n),
  site  = rep(c("north", "south", "east", "west"), times = c(300, 180, 90, 30)),
  team  = rep(paste0("t", 1:20), each = 30),
  when  = rep(seq(as.POSIXct("2024-03-01", tz = "UTC"), by = "day",
                  length.out = 15), each = 40),
  value = round(stats::rlnorm(n, meanlog = 6, sdlog = 1.1), 2)
)

# Two things we could only learn by opening the invoice: what the review cost
# (scales with size) and how many checks it needed (doesn't)
invoices$cost <- round(invoices$value * stats::runif(n, 0.02, 0.06), 2)
invoices$checks <- stats::rpois(n, 3) + 1

str(invoices)
#> 'data.frame':    600 obs. of  7 variables:
#>  $ id    : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ site  : chr  "north" "north" "north" "north" ...
#>  $ team  : chr  "t1" "t1" "t1" "t1" ...
#>  $ when  : POSIXct, format: "2024-03-01" "2024-03-01" ...
#>  $ value : num  1823 217 602 809 629 ...
#>  $ cost  : num  65.75 8.08 28 22.09 33.69 ...
#>  $ checks: num  6 3 3 8 4 4 3 2 8 5 ...
sum(invoices$value)
#> [1] 428704.5
```

The population total is what we’ll try to recover from a sample of 60.

## How many rows? Ask before you choose a design

Sixty is a stand-in for a number you should really derive.
[`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md)
inverts the usual question: rather than guessing `n` and hoping, say
what margin of error you can live with and let it solve for the size.

``` r

plan_size(margin = 100, sd = stats::sd(invoices$value), N = nrow(invoices))
#> Sample size for a mean
#>   draw           260
#>   margin         +/- 100 at 95% confidence
#>   assuming       sd 1091, N 600
```

The frame size matters — `N` applies a finite population correction, and
600 invoices is small enough that it bites hard. So do two other facts
about the real world:

``` r

# A clustered design carries less information per row, and people don't respond
plan_size(margin = 100, sd = stats::sd(invoices$value), N = nrow(invoices),
          deff = 2.5, response = 0.7)
#> Sample size for a mean
#>   draw           563
#>   to analyse     394  (after 70% response)
#>   margin         +/- 100 at 95% confidence
#>   assuming       sd 1091, deff 2.5, N 600
```

`deff` is not a guess you have to invent — it is measurable from a
comparable past sample with
[`deff()`](https://elkronos.github.io/dRawn/reference/deff.md), which we
come back to below. If the margin you ask for cannot be reached by
sampling at all,
[`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md)
says so rather than returning a number bigger than your frame:

``` r

plan_size(margin = 2, sd = stats::sd(invoices$value), N = nrow(invoices))
#> Sample size for a mean
#>   draw           600
#>   margin         +/- 2 at 95% confidence
#>   assuming       sd 1091, N 600
#> 
#> 
#>   That is a census: the frame is not large enough to reach this margin
#>   by sampling, so every row is needed. Widen `margin`, or accept the
#>   precision a full count gives.
```

We’ll stay with 60 for the rest of this vignette, because it makes the
comparisons between designs legible.

## A design is a value, not a function call

Designs are built, then applied. That separation is the point: the same
design can be printed, stored, reused on next month’s data, and — the
part that matters — asked what it does.

``` r

plan <- design_stratified(strata = "site", n = 60)
plan
#> <sampling design: stratified>
#>   strata           "site"
#>   n                60
#>   allocation       "proportional"
#>   min_per_stratum  0
#>   replace          FALSE
#>   na_rm            FALSE

sample_1 <- draw(invoices, plan, seed = 1)
table(sample_1$site)
#> 
#>  east north south  west 
#>     9    30    18     3
```

Proportional allocation splits 60 across sites in proportion to their
size, and the total is exact — no rounding drift.

``` r

nrow(sample_1)
#> [1] 60
```

## Estimating from the sample

A sample is only useful for estimation if you know how likely each row
was to be in it. `weights = TRUE` attaches that.

``` r

s <- draw(invoices, plan, seed = 1, weights = TRUE)
head(s[, c(".prob", ".weight", "site", "value")], 4)
#>   .prob .weight  site   value
#> 1   0.1      10 north  296.88
#> 2   0.1      10 north 1723.52
#> 3   0.1      10 north   56.86
#> 4   0.1      10 north 3244.49
```

`.weight` is `1 / .prob`: the number of population rows each sampled row
stands for. Multiply and add up, and you have a Horvitz–Thompson total.

``` r

ht <- sum(s$value * s$.weight)
c(estimate = round(ht), truth = sum(invoices$value))
#> estimate    truth 
#> 380046.0 428704.5
```

One draw won’t land on the truth.
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
gives you the total *and* a standard error, so you know how far off it
might be:

``` r

ht_total(s, "value")
#> Horvitz-Thompson total  (stratified design, n = 60)
#>   estimate 380,046.5
#>   se       49,499.7  (analytic)
#>   95% CI  283,028.9 to 477,064.1
#>   deff     0.994  (about the same as simple random sampling)
```

That interval comes from second-order inclusion probabilities — the
chance a *pair* of rows both land in the sample — which
[`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)
computes:

``` r

round(joint_prob(invoices, plan, rows = 1:4), 4)
#>        [,1]   [,2]   [,3]   [,4]
#> [1,] 0.1000 0.0097 0.0097 0.0097
#> [2,] 0.0097 0.1000 0.0097 0.0097
#> [3,] 0.0097 0.0097 0.1000 0.0097
#> [4,] 0.0097 0.0097 0.0097 0.1000
```

The estimator is unbiased — the average over many draws converges.

``` r

ests <- vapply(1:400, function(i) {
  si <- draw(invoices, plan, seed = i, weights = TRUE)
  sum(si$value * si$.weight)
}, numeric(1))

c(mean_estimate = round(mean(ests)), truth = sum(invoices$value))
#> mean_estimate         truth 
#>      435056.0      428704.5
```

And the reported standard error tracks the estimator’s real spread:

``` r

ses <- vapply(1:400, function(i) {
  ht_total(draw(invoices, plan, seed = i, weights = TRUE), "value")$se
}, numeric(1))

c(mean_reported_se = round(mean(ses)), actual_sd = round(stats::sd(ests)))
#> mean_reported_se        actual_sd 
#>            74637            80154
```

You can also ask the design about the population *before* drawing
anything, which is a quick way to check a plan is doing what you meant:

``` r

p <- inclusion_prob(invoices, plan)
tapply(p, invoices$site, unique)
#>  east north south  west 
#>   0.1   0.1   0.1   0.1
```

Every site is sampled at very nearly 10%, which is what proportional
allocation means.

## Means, and what the design cost you

[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md) is
the same machinery pointed at an average instead of a total:

``` r

ht_mean(s, "value")
#> Hajek mean  (stratified design, n = 60)
#>   estimate 633.4108
#>   se       82.49949  (analytic)
#>   95% CI  471.7148 to 795.1069
#>   deff     0.994  (about the same as simple random sampling)
mean(invoices$value)
#> [1] 714.5076
```

It defaults to the Hájek estimator, `sum(y/pi) / sum(1/pi)`, which
divides by the *estimated* population size rather than the known one.
`estimator = "ht"` divides by the true `N`. The two agree exactly
whenever the weights of the rows you drew sum to `N` — which includes
the stratified design above, so here they are the same number:

``` r

c(hajek = ht_mean(s, "value", variance = "none")$mean,
  ht    = ht_mean(s, "value", "ht", variance = "none")$mean)
#>    hajek       ht 
#> 633.4108 633.4108
```

They part company when the sample size is random. Poisson sampling is
the clean case — and which of the two wins depends on what you are
averaging:

``` r

pois <- design_weighted("value", n = 60, method = "poisson")

spread_both <- function(column) {
  b <- vapply(1:300, function(i) {
    si <- draw(invoices, pois, seed = i, weights = TRUE)
    c(hajek = ht_mean(si, column, variance = "none")$mean,
      ht    = ht_mean(si, column, "ht", variance = "none")$mean)
  }, numeric(2))
  apply(b, 1, stats::sd)
}

round(sapply(c(checks = "checks", cost = "cost", value = "value"),
             spread_both), 2)
#>       checks cost  value
#> hajek   0.41 5.14 129.38
#> ht      1.01 2.99  72.22
```

Selection here is driven by `value`. For `checks`, which has nothing to
do with invoice size, Hájek is more than twice as steady — a draw that
happens to catch too many low-probability rows inflates numerator and
denominator together, and they partly cancel. For `cost`, which scales
with size, and for `value` itself, that cancellation is exactly what you
*don’t* want: `y / pi` is already nearly constant, so dividing by the
known `N` wins.

Hájek is the default because the first case is the common one — you
usually select on size and measure something else. But it is a tendency,
not a law, and this is how you check which side of it you are on.

Every estimate also reports a **design effect**: the variance under this
design against simple random sampling of the same size. Read it as an
exchange rate on sample size — at `deff = 2`, a sample of 400 carries
the information of 200.

To see it move, we need a frame where the grouping variable genuinely
matters. Regions on very different value levels, and delivery routes
that each sit inside one region:

``` r

set.seed(11)
regions <- data.frame(
  id     = 1:400,
  region = rep(c("north", "south", "east", "west"), each = 100),
  route  = rep(paste0("r", 1:40), each = 10)
)
regions$value <- rep(c(120, 260, 480, 900), each = 100) +
  round(stats::rnorm(400, 0, 40))

by_region <- draw(regions, design_stratified("region", n = 40), seed = 1,
                  weights = TRUE)
by_route <- draw(regions, design_cluster("route", n_clusters = 4), seed = 1,
                 weights = TRUE)

c(stratified = deff(ht_mean(by_region, "value")),
  clustered  = deff(ht_mean(by_route,  "value")))
#>  stratified   clustered 
#>  0.02207873 12.86168571
```

Same population, same forty rows. Stratifying on what varies removes
nearly all the variance; clustering on it keeps nearly all of it. That
is the range a design choice moves you across — and that number is
exactly what `plan_size(deff = )` wants for the next study, which closes
the loop.

## Checking what you actually got

Before analysing, look at what came back against what was there:

``` r

sample_summary(s)
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

This one is healthy: every stratum came back at its intended 10%, and
the weights are constant. That is what you want to confirm before
analysing, because the three things that quietly ruin an estimate all
show up here — strata that came up short, design weights that vary far
more than expected, and rows the design could never have reached at all.
Temporal designs are the usual source of that last one, since anything
outside the sampling window has probability zero:

``` r

window <- design_temporal("when", from = "2024-03-01", to = "2024-03-08",
                          interval = 1, per_interval = 4, unit = "days")
sample_summary(draw(invoices, window, seed = 1, weights = TRUE))
#> Sample of 28 from 600  (temporal design)
#>   sampling fraction  0.0467
#>   design weights     10 to 10   (cv 0)
#>   unreachable rows   320  (probability 0 under this design)
```

## When you want the rare stratum represented

`west` has only 30 invoices. Proportional allocation gives it three. If
you need coverage rather than efficiency, ask for a floor — and note
that this is a deliberate bias, not a free lunch:

``` r

covered <- draw(invoices, design_stratified("site", n = 60, min_per_stratum = 8),
                seed = 1)
table(covered$site)
#> 
#>  east north south  west 
#>     8    28    16     8
```

The default of `min_per_stratum = 0` leaves allocation unbiased. Setting
it over-represents small strata, which is fine when the goal is “look at
every site” and wrong when the goal is “estimate the total”.

## Sampling proportional to size

For an audit you usually want large invoices sampled more often. The
obvious call is not quite the right one:

``` r

naive <- design_weighted("value", n = 60)
naive
#> <sampling design: weighted>
#>   weights  "value"
#>   n        60
#>   replace  FALSE
#>   method   "successive"
#>   na_rm    FALSE
```

This is what `base::sample(prob = )` does. It’s a perfectly good way to
bias selection toward large invoices, but the weights govern each
sequential *draw*, not the probability of ending up in the sample — so
it has no closed-form inclusion probability, and the package refuses to
invent one:

``` r

inclusion_prob(invoices, naive)
#> Error:
#> ! `design_weighted(method = "successive")` has no closed-form inclusion probability.
#> Successive sampling has no closed-form inclusion probability.
#> Use method = "systematic" or "poisson" for a design whose inclusion
#> probabilities really are proportional to the weights.
#> Or pass simulate = TRUE to estimate it by Monte Carlo.
```

Ask for a design whose inclusion probabilities really are proportional
to size:

``` r

pps <- design_weighted("value", n = 60, method = "systematic")
pi <- inclusion_prob(invoices, pps)
round(head(pi / invoices$value, 5), 8)
#> [1] 0.00014386 0.00014386 0.00014386 0.00014386 0.00014386
```

A constant ratio to `value` is what “proportional to size” means. And
now the estimator works:

``` r

ests_pps <- vapply(1:400, function(i) {
  si <- draw(invoices, pps, seed = i, weights = TRUE)
  sum(si$value * si$.weight)
}, numeric(1))

c(mean_estimate = round(mean(ests_pps)), truth = sum(invoices$value))
#> mean_estimate         truth 
#>      428705.0      428704.5
```

Compare its spread against the stratified estimator’s, and something
surprising happens:

``` r

c(stratified_sd = round(stats::sd(ests)),
  pps_sd        = round(stats::sd(ests_pps)))
#> stratified_sd        pps_sd 
#>         80154             0
```

Zero. Not “very small” — exactly zero, in every sample. That is not a
fluke, it is the limiting case, and it is worth understanding because it
says exactly what PPS is for. Under systematic PPS the inclusion
probability is `pi_i = n * value_i / sum(value)`, so each sampled row
contributes `value_i / pi_i = sum(value) / n` — the same amount,
whichever rows you happen to draw. Totalling the very variable that
drove selection gives the answer back exactly.

Real work is never that kind. You select on a size measure you already
hold for the whole frame, and you total something you can only observe
once you open the file. Whether PPS helps then depends entirely on
whether the two are related:

``` r

spread_of <- function(design, column) {
  stats::sd(vapply(1:400, function(i) {
    si <- draw(invoices, design, seed = i, weights = TRUE)
    sum(si[[column]] * si$.weight)
  }, numeric(1)))
}

round(rbind(
  cost   = c(stratified = spread_of(plan, "cost"),
             pps        = spread_of(pps,  "cost")),
  checks = c(stratified = spread_of(plan, "checks"),
             pps        = spread_of(pps,  "checks"))
), 1)
#>        stratified   pps
#> cost       3360.9 517.3
#> checks      133.1 495.4
```

Review cost scales with invoice size, and PPS estimates its total
several times more precisely than stratification does. The number of
checks an invoice needs does not scale with size — and there PPS is
markedly *worse*, because it spends most of the sample on large invoices
that carry no more information about checks than small ones do.

That is the trade PPS actually offers. It is a bet that your size
measure predicts what you are measuring, and it is a bad bet when it
doesn’t. [`deff()`](https://elkronos.github.io/dRawn/reference/deff.md),
below, puts a number on which way the bet went.

## The units too big to leave to chance

PPS is one answer to dominant units. The blunter one is to stop sampling
them at all: take every invoice above a threshold, and sample what is
left.

``` r

audit <- design_certainty("value", threshold = 2000,
                          rest = design_stratified("site", n = 60))
audit
#> <sampling design: certainty>
#>   above      "value"
#>   threshold  2000
#>   rest       <drawn_design_stratified[7]>
#>   na_rm      FALSE

s_audit <- draw(invoices, audit, seed = 1, weights = TRUE)
table(certain = s_audit$.prob == 1)
#> certain
#> FALSE  TRUE 
#>    60    36
```

`rest` takes any design, and its `n` counts rows drawn from *below* the
threshold, not the total. The rows above it have inclusion probability
exactly 1 — so they carry a weight of 1, contribute their own value to
the total, and add **nothing** to its variance:

``` r

ht_total(s_audit, "value")
#> Horvitz-Thompson total  (certainty design, n = 96)
#>   estimate 458,373.7
#>   se       34,691.5  (analytic)
#>   95% CI  390,379.6 to 526,367.8
#>   deff     0.314  (better than simple random sampling)
```

All the uncertainty comes from the part that was actually sampled. That
is the whole point, and it is why this is the standard shape of an audit
or financial sample. Compare it against stratifying alone on the same
frame:

``` r

c(stratified = ht_total(s, "value")$se,
  certainty  = ht_total(s_audit, "value")$se)
#> stratified  certainty 
#>    49499.7    34691.5
```

Inclusion and joint probabilities compose exactly across the two parts —
a certainty row co-occurs with another row exactly when that row is
drawn, and two certainty rows always co-occur — so that standard error
is analytic, not an approximation.

## The designs at a glance

``` r

designs <- list(
  simple      = design_simple(n = 60),
  stratified  = design_stratified("site", n = 60),
  systematic  = design_systematic(interval = 10),
  cluster     = design_cluster("team", n_clusters = 4),
  multistage  = design_multistage("team", n_clusters = 6, n = 60),
  weighted    = design_weighted("value", n = 60, method = "systematic"),
  certainty   = design_certainty("value", 2000, design_simple(n = 60)),
  temporal    = design_temporal("when", from = "2024-03-01", to = "2024-03-15",
                                interval = 1, per_interval = 4, unit = "days")
)

vapply(designs, function(d) nrow(draw(invoices, d, seed = 1)), numeric(1))
#>     simple stratified systematic    cluster multistage   weighted  certainty 
#>         60         60         60        120         60         60         96 
#>   temporal 
#>         56
```

[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)
has no `n`: you choose how many clusters, and the row count follows from
which ones you get. Use
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
when you need to control both.

## Which designs can be estimated from

Most designs have exact inclusion probabilities. Five situations do not,
and the package says so rather than approximating:

``` r

no_form <- list(
  `balanced clusters`      = design_cluster("team", n_clusters = 4, balanced = TRUE),
  `proportional multistage`= design_multistage("team", n_clusters = 6, n = 60,
                                               allocation = "proportional"),
  `uneven multistage`      = design_multistage("team", n_clusters = 6, n = 61),
  `successive weighting`   = design_weighted("value", n = 60),
  `bootstrap`              = design_bootstrap(n_replicates = 10)
)

for (nm in names(no_form)) {
  msg <- tryCatch(inclusion_prob(invoices, no_form[[nm]]),
                  error = function(e) conditionMessage(e))
  cat(nm, "->", strsplit(msg, "\n")[[1]][1], "\n")
}
#> balanced clusters -> `design_cluster(balanced = TRUE)` has no closed-form inclusion probability. 
#> proportional multistage -> `design_multistage(allocation = "proportional")` has no closed-form inclusion probability. 
#> uneven multistage -> `design_multistage()` with n = 61 over 6 clusters has no closed-form inclusion probability. 
#> successive weighting -> `design_weighted(method = "successive")` has no closed-form inclusion probability. 
#> bootstrap -> A bootstrap resamples the sample; it is not a probability sample of a finite population, so it has no inclusion probability.
```

In each case the probability depends on something random — which
clusters were picked, or the order draws happened in. The multistage
pair is worth dwelling on, because the arithmetic looks innocent. The
second-stage allocation runs over the clusters that were *selected*, so
a leftover row goes to the largest of them, and a cluster smaller than
its target is capped with the shortfall dealt to whichever clusters came
with it. A row’s chance therefore depends on the company it keeps.
Averaging over the possibilities is not a closed form, and the average
is not the answer.

Where you need a number anyway, ask for a simulation and accept the
Monte Carlo error:

``` r

sim <- inclusion_prob(invoices, no_form$`balanced clusters`,
                      simulate = TRUE, R = 200, seed = 1)
round(tapply(sim, invoices$team, unique)[1:5], 3)
#>    t1   t10   t11   t12   t13 
#> 0.190 0.200 0.225 0.155 0.205
```

The one thing simulation cannot rescue is the bootstrap, and it is
refused rather than answered — every row turns up in some replicate, so
the count would converge to 1 for all of them.

``` r

inclusion_prob(invoices, no_form$bootstrap, simulate = TRUE, R = 100)
#> Error:
#> ! A bootstrap is not a probability sample of a finite population, so simulating it does not estimate an inclusion probability: every row appears in some replicate, and the count converges to 1 for all of them.
#> If you want the expected number of times each row appears, that is n_replicates * n / nrow(data).
```

Second-order probabilities can be simulated too, which is the general
answer where no formula exists. Systematic PPS is the case that needs
it: it has exact first-order probabilities, but its joint ones depend on
the order units are visited.

``` r

joint_prob(invoices, pps, rows = 1:4)
#> Error:
#> ! `design_weighted(method = "systematic")` has no closed-form joint inclusion probability.
#> Its joint probabilities depend on the order units are visited and need a
#> dedicated algorithm. `sampling::UPsystematicpi2()` computes them.
```

``` r

round(joint_prob(invoices, pps, rows = 1:4, simulate = TRUE, R = 2000,
                 seed = 1), 4)
#>        [,1]   [,2]   [,3]   [,4]
#> [1,] 0.2520 0.0115 0.0225 0.0275
#> [2,] 0.0115 0.0390 0.0035 0.0050
#> [3,] 0.0225 0.0035 0.0795 0.0110
#> [4,] 0.0275 0.0050 0.0110 0.1190
```

The diagonal holds the first-order probabilities, and simulating them is
a useful check on the closed form:

``` r

round(cbind(
  simulated = diag(joint_prob(invoices, pps, rows = 1:4, simulate = TRUE,
                              R = 2000, seed = 1)),
  exact     = inclusion_prob(invoices, pps)[1:4]
), 4)
#>      simulated  exact
#> [1,]    0.2520 0.2622
#> [2,]    0.0390 0.0312
#> [3,]    0.0795 0.0865
#> [4,]    0.1190 0.1164
```

## Clusters, when fieldwork costs more than precision

Visiting five teams is cheaper than visiting twenty.
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)
takes whole clusters, so the row count follows from which ones you get:

``` r

by_team <- design_cluster("team", n_clusters = 4)
res <- draw(invoices, by_team, seed = 1)
c(teams = length(unique(res$team)), rows = nrow(res))
#> teams  rows 
#>     4   120
```

[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
adds a second stage when you need to control the total:

``` r

ms <- design_multistage("team", n_clusters = 6, n = 60)
res_ms <- draw(invoices, ms, seed = 1)
table(res_ms$team)
#> 
#>  t1 t13 t19  t2  t4  t7 
#>  10  10  10  10  10  10
```

Both are estimable, but cluster designs pay for their convenience in
precision — and in degrees of freedom. With only a handful of clusters
the normal-approximation interval undercovers, so treat it as
indicative:

``` r

ht_total(draw(invoices, by_team, seed = 1, weights = TRUE), "value")
#> Horvitz-Thompson total  (cluster design, n = 120)
#>   estimate 446,199.5
#>   se       46,834.3  (analytic)
#>   95% CI  354,406 to 537,993.1
#>   deff     0.783  (better than simple random sampling)
```

## Streams, when the data does not fit

A data frame is not a stream: its length is already known, so
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md)
takes a direct vectorised path for one. The streaming path is for data
you cannot hold — pass a connection, or a function that returns the next
item and `NULL` when exhausted.

``` r

i <- 0
next_invoice <- function() {
  i <<- i + 1
  if (i > 1e5) NULL else i
}

drawn_ids <- unlist(draw(next_invoice, design_reservoir(n = 5), seed = 1))
drawn_ids
#> [1] 93918 86723 83503 32995 42017
```

One pass, five items held, a hundred thousand seen. Algorithm L skips
ahead geometrically rather than drawing a random number per item.

## Bootstrap, for the uncertainty of a statistic

The other designs estimate a population total. The bootstrap estimates
the sampling distribution of whatever you like — here the median, which
has no convenient closed form:

``` r

reps <- draw(invoices, design_bootstrap(n_replicates = 400), seed = 1)
meds <- vapply(split(reps, reps$.replicate), function(r) median(r$value),
               numeric(1))

c(observed = median(invoices$value),
  boot_se  = round(stats::sd(meds), 1),
  lower    = round(unname(stats::quantile(meds, 0.025)), 1),
  upper    = round(unname(stats::quantile(meds, 0.975)), 1))
#> observed  boot_se    lower    upper 
#>   387.41    16.60   351.60   420.70
```

Replicates come back in one frame with a leading `.replicate` column.
For ordered data where neighbouring rows are correlated,
`method = "block"` concatenates randomly chosen runs instead of
independent rows:

``` r

design_bootstrap(n_replicates = 100, method = "block", block_length = 20)
#> <sampling design: bootstrap>
#>   n_replicates  100
#>   method        "block"
#>   block_length  20
```

## Reproducibility without side effects

`seed` is local to the draw. Sampling inside a simulation does not shift
the simulation’s own random stream:

``` r

set.seed(99)
before <- runif(3)

set.seed(99)
invisible(draw(invoices, plan, seed = 7))
after <- runif(3)

identical(before, after)
#> [1] TRUE
```

That is deliberate: a sampling call is not supposed to be visible to the
code around it.

## Where this package stops

`drawn` draws samples and estimates totals and means from them. It does
not do subpopulation estimates, regression, calibration, or quantiles
with proper standard errors —
[`survey`](https://cran.r-project.org/package=survey) does all of that,
and
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
hands a sample straight over.

``` r

library(survey)

des <- as_svydesign(s)
svytotal(~value, des)
#>        total    SE
#> value 380046 49500
```

The mapping is direct: design weights become `weights`, the
stratification or clustering column becomes `strata` or `ids`, and the
population counts become `fpc` — so `survey` applies the same finite
population correction this package does, and the two agree.

``` r

c(drawn  = ht_total(s, "value")$se,
  survey = as.numeric(SE(svytotal(~value, des))))
#>   drawn  survey 
#> 49499.7 49499.7
```

From there it is `survey`’s vocabulary:

``` r

svyby(~value, ~site, des, svymean)
#>        site    value       se
#> east   east 853.5622 246.9130
#> north north 577.2890 116.2225
#> south south 693.3406 150.9955
#> west   west 174.5967  47.0911
```

Draw here, analyse there.
