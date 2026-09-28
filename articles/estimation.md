# How the estimates work

``` r

library(drawn)
```

This article sets out what
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md)
and [`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md)
compute, which variance estimator each design gets and why, and how each
one was checked. Everything here is design-based: the population values
are fixed, and the only randomness is which rows the design selected
(Särndal, Swensson and Wretman 1992).

## The estimators

A design gives every frame row \\k\\ an inclusion probability \\\pi_k\\
— the chance it lands in the sample — and every pair of rows a joint
probability \\\pi\_{kl}\\.
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
and
[`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)
compute them from the design and the frame, before anything is drawn.

**The total** is the Horvitz-Thompson estimator (Horvitz and Thompson
1952):

\\\hat Y = \sum\_{k \in s} \frac{y_k}{\pi_k}.\\

It is unbiased for the frame total *provided every row has \\\pi_k \>
0\\*. That proviso is why
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) warns
when a stratum is allocated no rows, and why
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
counts unreachable rows.

**The mean** defaults to the Hájek estimator (Hájek 1971), a ratio of
two Horvitz-Thompson totals:

\\\hat{\bar Y} = \frac{\sum\_{s} y_k / \pi_k}{\sum\_{s} 1 / \pi_k}.\\

Its variance is the variance of the Horvitz-Thompson total of the
linearised residuals \\z_k = (y_k - \hat{\bar Y}) / \hat N\\ (Deville
1999), so every variance estimator below serves means as well as totals.
A proportion is the mean of a 0/1 variable.

**Domains** — sites, regions, months — use the whole sample with \\y\\
set to zero outside the domain: \\z_k = y_k \cdot 1\[k \in d\]\\ for a
total, and the corresponding residual for a mean (Särndal et al. 1992,
section 10.3). Cutting the sample into pieces first would treat the
number of rows that fell in each piece as fixed, which it is not.

## Which variance estimator, for which design

| Design | Estimator | `method` | Status |
|----|----|----|----|
| simple, stratified, temporal, spatial, reservoir, multistage | Sen-Yates-Grundy with exact \\\pi\_{kl}\\ | `"analytic"` | exact, unbiased |
| cluster | between-cluster form over cluster totals | `"analytic"` | exact, unbiased |
| weighted, `"poisson"` | \\\sum (1-\pi_k) y_k^2 / \pi_k^2\\ | `"analytic"` | exact, unbiased |
| weighted, `"systematic"` | Deville (1999) | `"deville"` | approximation |
| systematic | successive differences (Wolter 2007) | `"successive difference"` | approximation |
| spread | local means (Stevens and Olsen 2003) | `"local mean"` | approximation, conservative |
| certainty | whatever `rest` uses | as `rest` | as `rest` |
| any, on request or as fallback | jackknife, stratified where there are strata | `"jackknife"` | approximation |

### Sen-Yates-Grundy

For a fixed-size design with known joint probabilities,

\\\hat V\_{\mathrm{SYG}} = \frac12 \sum\_{k \ne l \in s} \frac{\pi_k
\pi_l - \pi\_{kl}}{\pi\_{kl}} \left(\frac{y_k}{\pi_k} -
\frac{y_l}{\pi_l}\right)^2\\

is unbiased whenever every \\\pi\_{kl} \> 0\\ (Sen 1953; Yates and
Grundy 1953). For stratified simple random sampling it reduces to the
textbook \\\sum_h N_h^2 (1 - f_h) s_h^2 / n_h\\.

The condition matters. **A stratum with one sampled row** has no pair
inside it, so it silently contributes nothing: with one row in every
stratum the sum is exactly zero. `drawn` declines in that case rather
than report it:

``` r

set.seed(1)
pop <- data.frame(g = rep(1:10, each = 20), y = rnorm(200, rep(1:10, each = 20), 3))
one_each <- draw(pop, design_stratified("g", n = 10), seed = 1, weights = TRUE)
ht_total(one_each, "y")
#> Horvitz-Thompson total  (stratified design, n = 10)
#>   estimate 957.7274
#>   se       NA
#> 
#>   No variance is available for this sample. 10 strata `1`, `10`, `2`,
#>   `3`, `4`, ... each have a single sampled row, so the variation within
#>   them cannot be measured and leaving it out would understate the
#>   standard error. Draw at least two rows per stratum (min_per_stratum =
#>   2). The stratified jackknife needs two units per stratum as well.
```

### Cluster sampling

Taking whole clusters makes the number of *rows* random whenever
clusters differ in size, which violates the fixed-size assumption behind
Sen-Yates-Grundy. The cluster is the sampling unit, so the estimator is
applied to cluster totals \\u_c = \sum\_{k \in c} y_k / \pi_k\\: \\\hat
V = (1 - a/A)\\ a\\ s_u^2\\, the simple-random form over \\a\\ of \\A\\
clusters.

### Poisson sampling

Rows are selected independently, so \\\pi\_{kl} = \pi_k \pi_l\\ and the
Horvitz-Thompson variance collapses to \\\sum_s (1 - \pi_k) y_k^2 /
\pi_k^2\\.

### Systematic PPS: Deville’s approximation

`design_weighted(method = "systematic")` shuffles the frame and walks
the cumulative probabilities (Madow 1949; Hartley and Rao 1962). Its
joint probabilities have no closed form. Deville’s estimator needs only
the first-order ones:

\\\hat V = \frac{1}{1 - \sum_s a_k^2} \sum\_{k \in s} c_k
\left(\frac{y_k}{\pi_k} - \hat A\right)^2, \quad c_k = 1 - \pi_k,\\ a_k
= \frac{c_k}{\sum_s c_l},\\ \hat A = \sum_s a_k \frac{y_k}{\pi_k}.\\

Each row is weighted by its own \\1 - \pi_k\\, so rows taken with
certainty drop out, as they should. Matei and Tillé (2005) found it the
most reliable of the standard approximations for high-entropy designs.
It equals
[`sampling::varest()`](https://rdrr.io/pkg/sampling/man/varest.html).

### Systematic sampling: successive differences

A systematic sample has one random start. Most pairs of rows can never
be drawn together (\\\pi\_{kl} = 0\\), so no design-unbiased variance
estimator exists at all. The successive-difference estimator compares
each sampled row with its neighbour in the order the design walked
(Wolter 2007):

\\\hat V = (1 - f)\\ \frac{n}{2(n-1)} \sum\_{j=2}^{n}
\left(\frac{y\_{(j)}}{\pi} - \frac{y\_{(j-1)}}{\pi}\right)^2.\\

Sorting the frame on something related to \\y\\ is what makes systematic
sampling efficient, and this estimator can see that where the
simple-random formula cannot:

``` r

set.seed(12)
frame <- data.frame(x = runif(400, 0, 100))
frame$y <- frame$x + rnorm(400, 0, 20)
sorted <- design_systematic(interval = 10, order_by = "x")

# The truth, by enumerating all ten starts
starts <- vapply(1:10, function(st) {
  10 * sum(draw(frame, design_systematic(10, start = st, order_by = "x"))$y)
}, numeric(1))
true_var <- mean((starts - mean(starts))^2)

est <- vapply(1:200, function(i) {
  s <- draw(frame, sorted, seed = i, weights = TRUE)
  c(successive = ht_total(s, "y")$variance,
    as_if_srs  = 400^2 * 0.9 * var(s$y) / nrow(s))
}, numeric(2))
round(rowMeans(est) / true_var, 2)
#> successive  as_if_srs 
#>       1.59       5.48
```

Both are ratios to the true variance; the successive-difference figure
is far closer. Its limits are shared by every estimator from a single
systematic sample: a frame that cycles with a period matching the
interval, and a trend so smooth that the random start is almost the only
source of variation, both make it run low.

### Spatially balanced samples: local means

[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md)
selects by the local pivotal method (Grafström, Lundström and Schelin
2012), which makes neighbours unlikely to be sampled together. That is
exactly the gain, and a simple-random variance would throw it away. The
local mean estimator compares each sampled row with the mean over its
neighbourhood \\N_k\\ — itself and its three nearest sampled neighbours:

\\\hat V = \sum\_{k \in s} \frac{\|N_k\|}{\|N_k\| - 1}
\left(\frac{y_k}{\pi_k} - \frac{1}{\|N_k\|}\sum\_{l \in N_k}
\frac{y_l}{\pi_l}\right)^2\\

(Stevens and Olsen 2003; Grafström and Schelin 2014). It errs
conservative, more so as the sampling fraction grows.

### The jackknife

The delete-one-unit jackknife works from the sample alone. Units are
clusters where the design has them and rows otherwise, and where the
design has strata they are deleted *within* their stratum, with only
that stratum reweighted — the stratified or JKn jackknife (Wolter 2007,
chapter 4). Deleting across strata would count the differences between
strata as sampling variance, which is what stratifying removes. Each
stratum carries its own finite population correction. For a stratified
design the result equals the analytic figure exactly:

``` r

strat <- draw(pop, design_stratified("g", n = 40), seed = 2, weights = TRUE)
c(analytic  = ht_total(strat, "y", variance = "analytic")$se,
  jackknife = ht_total(strat, "y", variance = "jackknife")$se)
#>  analytic jackknife 
#>  78.99364  78.99364
```

`variance = "auto"` uses the design’s own estimator, and the jackknife
only when that declines — for instance a multistage sample with one row
per cluster, where the within-cluster variation cannot be measured and
the delete-a-cluster jackknife gives the standard ultimate-cluster
approximation. It is refused for systematic samples (one primary
sampling unit) and Poisson samples (whose exact variance it
understates).

## Confidence intervals

Intervals use Student’s \\t\\ on the design’s degrees of freedom —
primary sampling units minus strata — rather than the normal
distribution (Korn and Graubard 1999). With hundreds of units the two
agree; with a handful of clusters the normal interval is badly short:

``` r

set.seed(5)
clus <- data.frame(cl = rep(1:24, each = 10))
clus$y <- rep(rnorm(24, 50, 20), each = 10) + rnorm(240, 0, 5)
truth <- sum(clus$y)

hits <- vapply(1:500, function(i) {
  r <- ht_total(draw(clus, design_cluster("cl", n_clusters = 4), seed = i,
                     weights = TRUE), "y")
  c(t = r$ci[1] <= truth && truth <= r$ci[2],
    normal = abs(r$total - truth) <= qnorm(0.975) * r$se)
}, logical(2))
rowMeans(hits)
#>      t normal 
#>  0.960  0.846
```

Pass `df = Inf` for the normal interval.

## How it was checked

Every estimator is checked in the package’s tests against the empirical
sampling variance of its own estimator over many repeated draws, and
every inclusion and joint probability against the observed frequency. A
condensed version, run when this page was built — the ratio of the
average estimated variance to the actual variance of the estimates (1 is
perfect), and the coverage of the nominal 95% interval:

``` r

set.seed(7)
N <- 300
frame <- data.frame(
  site = rep(c("a", "b", "c"), times = c(150, 100, 50)),
  cl   = rep(paste0("c", 1:30), each = 10),
  x    = runif(N), yy = runif(N)
)
frame$size <- rgamma(N, 3, 0.2)
frame$y <- 2 * frame$size + 30 * frame$x + rep(rnorm(30, 0, 8), each = 10) +
  rnorm(N, 0, 6)
truth <- sum(frame$y)

designs <- list(
  simple        = design_simple(n = 40),
  stratified    = design_stratified("site", n = 40),
  cluster       = design_cluster("cl", n_clusters = 6),
  multistage    = design_multistage("cl", n_clusters = 8, n = 32),
  pps_poisson   = design_weighted("size", n = 40, method = "poisson"),
  pps_systematic = design_weighted("size", n = 40, method = "systematic"),
  spread        = design_spread(c("x", "yy"), n = 40, scale = FALSE)
)

check <- function(d, R = 400) {
  out <- vapply(seq_len(R), function(i) {
    r <- ht_total(draw(frame, d, seed = i, weights = TRUE), "y")
    c(r$total, r$variance, r$ci[1] <= truth && truth <= r$ci[2])
  }, numeric(3))
  c(variance_ratio = mean(out[2, ]) / var(out[1, ]),
    coverage = mean(out[3, ]))
}
round(t(vapply(designs, check, numeric(2))), 2)
#>                variance_ratio coverage
#> simple                   1.00     0.93
#> stratified               1.05     0.95
#> cluster                  0.85     0.94
#> multistage               0.97     0.96
#> pps_poisson              0.87     0.93
#> pps_systematic           0.96     0.93
#> spread                   1.22     0.96
```

With a few hundred replicates the ratios carry Monte Carlo error of
several percent; the test suite uses more.

## Agreement with `survey`

[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
expresses each design in `survey`’s own terms (Lumley 2004), and the
tests assert that the two packages agree to floating point on totals,
standard errors, domain estimates and degrees of freedom for every
design `survey` can represent exactly — and within a fraction of a
percent for systematic PPS, where `survey` uses Brewer’s approximation
and `drawn` uses Deville’s.
[`?as_svydesign`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
has the full table.

## References

Deville, J.-C. (1999). Variance estimation for complex statistics and
estimators: linearization and residual techniques. *Survey Methodology*,
25, 193–203.

Grafström, A., Lundström, N. L. P. and Schelin, L. (2012). Spatially
balanced sampling through the pivotal method. *Biometrics*, 68, 514–520.

Grafström, A. and Schelin, L. (2014). How to select representative
samples. *Scandinavian Journal of Statistics*, 41, 277–290.

Hájek, J. (1971). Comment on “An essay on the logical foundations of
survey sampling, part one” by D. Basu. In *Foundations of Statistical
Inference*, Godambe and Sprott (eds.), p. 236. Holt, Rinehart and
Winston.

Hartley, H. O. and Rao, J. N. K. (1962). Sampling with unequal
probabilities and without replacement. *Annals of Mathematical
Statistics*, 33, 350–374.

Horvitz, D. G. and Thompson, D. J. (1952). A generalization of sampling
without replacement from a finite universe. *Journal of the American
Statistical Association*, 47, 663–685.

Korn, E. L. and Graubard, B. I. (1999). *Analysis of Health Surveys*.
Wiley.

Lumley, T. (2004). Analysis of complex survey samples. *Journal of
Statistical Software*, 9, 1–19.

Madow, W. G. (1949). On the theory of systematic sampling, II. *Annals
of Mathematical Statistics*, 20, 333–354.

Matei, A. and Tillé, Y. (2005). Evaluation of variance approximations
and estimators in maximum entropy sampling with unequal probability and
fixed sample size. *Journal of Official Statistics*, 21, 543–570.

Särndal, C.-E., Swensson, B. and Wretman, J. (1992). *Model Assisted
Survey Sampling*. Springer.

Sen, A. R. (1953). On the estimate of the variance in sampling with
varying probabilities. *Journal of the Indian Society of Agricultural
Statistics*, 5, 119–127.

Stevens, D. L. and Olsen, A. R. (2003). Variance estimation for
spatially balanced samples of environmental resources. *Environmetrics*,
14, 593–610.

Wolter, K. M. (2007). *Introduction to Variance Estimation*, 2nd
ed. Springer.

Yates, F. and Grundy, P. M. (1953). Selection without replacement from
within strata with probability proportional to size. *Journal of the
Royal Statistical Society B*, 15, 253–261.
