# Choosing a sampling design

``` r

library(drawn)
```

Every function in this package answers one question — *which rows do I
look at?* — and the interesting part is what you can say afterwards
about the rows you didn’t.

## One population, ten questions

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

str(invoices)
#> 'data.frame':    600 obs. of  5 variables:
#>  $ id   : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ site : chr  "north" "north" "north" "north" ...
#>  $ team : chr  "t1" "t1" "t1" "t1" ...
#>  $ when : POSIXct, format: "2024-03-01" "2024-03-01" ...
#>  $ value: num  1823 217 602 809 629 ...
sum(invoices$value)
#> [1] 428704.5
```

The population total is what we’ll try to recover from a sample of 60.

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
#>   total    380,046.5
#>   se       49,499.7  (analytic)
#>   95% CI  283,028.9 to 477,064.1
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

Because selection probability tracks the thing being totalled, the PPS
estimator is also far more precise than the stratified one here:

``` r

c(stratified_sd = round(stats::sd(ests)),
  pps_sd        = round(stats::sd(ests_pps)))
#> stratified_sd        pps_sd 
#>         80154             0
```

## The designs at a glance

``` r

designs <- list(
  simple      = design_simple(n = 60),
  stratified  = design_stratified("site", n = 60),
  systematic  = design_systematic(interval = 10),
  cluster     = design_cluster("team", n_clusters = 4),
  multistage  = design_multistage("team", n_clusters = 6, n = 60),
  weighted    = design_weighted("value", n = 60, method = "systematic"),
  temporal    = design_temporal("when", from = "2024-03-01", to = "2024-03-15",
                                interval = 1, per_interval = 4, unit = "days")
)

vapply(designs, function(d) nrow(draw(invoices, d, seed = 1)), numeric(1))
#>     simple stratified systematic    cluster multistage   weighted   temporal 
#>         60         60         60        120         60         60         56
```

[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md)
has no `n`: you choose how many clusters, and the row count follows from
which ones you get. Use
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md)
when you need to control both.

## Which designs can be estimated from

Nine of the ten have exact inclusion probabilities. Four situations do
not, and the package says so rather than approximating:

``` r

no_form <- list(
  `balanced clusters`      = design_cluster("team", n_clusters = 4, balanced = TRUE),
  `proportional multistage`= design_multistage("team", n_clusters = 6, n = 60,
                                               allocation = "proportional"),
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
#> successive weighting -> `design_weighted(method = "successive")` has no closed-form inclusion probability. 
#> bootstrap -> A bootstrap resamples the sample; it is not a probability sample of a finite population, so it has no inclusion probability.
```

In each case the reason is that the probability depends on something
random — which clusters were picked, or the order draws happened in.
Where you need a number anyway, ask for a simulation and accept the
Monte Carlo error:

``` r

sim <- inclusion_prob(invoices, no_form$`balanced clusters`,
                      simulate = TRUE, R = 200, seed = 1)
round(tapply(sim, invoices$team, unique)[1:5], 3)
#>    t1   t10   t11   t12   t13 
#> 0.190 0.200 0.225 0.155 0.205
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
#>   total    446,199.5
#>   se       46,834.3  (analytic)
#>   95% CI  354,406 to 537,993.1
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
