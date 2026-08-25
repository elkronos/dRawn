# drawn

Design-based sampling from data frames.

You describe a design, then draw from it:

```r
draw(data, design_stratified(strata = "site", n = 500), seed = 1)
```

Ten designs share one contract. Arguments mean the same thing everywhere, inputs
are validated before anything is sampled, the caller's random number stream is
left where it was found, and what comes back has the same class and the same
columns in the same order as what went in.

Because the design is an object rather than a pile of arguments, it can also say
how likely each row was to be selected — which is what turns a sample into
something you can estimate from:

```r
s <- draw(data, design_stratified("site", n = 500), seed = 1, weights = TRUE)
sum(s$spend * s$.weight)   # an unbiased Horvitz-Thompson total
```

Formerly `sampleR`, then briefly `sampleframe`. The original name had to change
because [`sampler`](https://cran.r-project.org/package=sampler) has been on CRAN
since 2018 in the same problem domain, and *Writing R Extensions* asks that
package names not differ from an existing one only by case. (`drawr` would have
hit the same rule against [`DRaWR`](https://cran.r-project.org/package=DRaWR).)

## Installation

```r
# install.packages("remotes")
remotes::install_github("elkronos/sampleR")
```

`sf` is only needed for `design_spatial()` and lives in `Suggests`, so a plain
install does not pull in GDAL, GEOS and PROJ.

## Usage

```r
library(drawn)

df <- data.frame(
  id    = 1:1000,
  value = rnorm(1000),
  site  = rep(paste0("s", 1:20), each = 50),
  ts    = seq(as.POSIXct("2024-01-01", tz = "UTC"), by = "hour", length.out = 1000)
)

# 100 rows, uniformly
draw(df, design_simple(n = 100), seed = 1)

# 100 rows split across sites in proportion to their size
draw(df, design_stratified("site", n = 100), seed = 1)

# 5 whole sites
draw(df, design_cluster("site", n_clusters = 5), seed = 1)

# 2 rows from each 6-hour window
draw(df, design_temporal("ts", from = "2024-01-01", to = "2024-01-08",
                         interval = 6, per_interval = 2, unit = "hours"), seed = 1)
```

A design is a value. Build it once, print it, pass it around, reuse it:

```r
monthly_audit <- design_stratified("site", n = 200, min_per_stratum = 1)

monthly_audit
#> <sampling design: stratified>
#>   strata           "site"
#>   n                200
#>   allocation       "proportional"
#>   min_per_stratum  1
#>   replace          FALSE
#>   na_rm            FALSE

jan <- draw(january_data, monthly_audit, seed = 1)
feb <- draw(february_data, monthly_audit, seed = 2)
```

## The designs

| Constructor | Draws | Key arguments |
|---|---|---|
| `design_simple()` | Rows uniformly at random | `n`, `replace` |
| `design_stratified()` | A share of each stratum | `strata`, `n`, `allocation`, `min_per_stratum` |
| `design_systematic()` | Every *k*-th row | `interval`, `start`, `order_by` |
| `design_cluster()` | Whole clusters | `clusters`, `n_clusters`, `balanced` |
| `design_multistage()` | Clusters, then rows within them | `clusters`, `n_clusters`, `n`, `allocation` |
| `design_weighted()` | Rows with probability by weight | `weights`, `n`, `replace` |
| `design_reservoir()` | A fixed-size sample from a stream | `n`, `max_items` |
| `design_bootstrap()` | Resampled replicates | `n_replicates`, `n`, `method`, `block_length` |
| `design_temporal()` | A share of each time interval | `time`, `from`, `to`, `interval`, `per_interval` |
| `design_spatial()` | Rows inside a region | `coords`, `region`, `n`, `crs` |

## What the shared contract means in practice

**`n` is always a total.** Never a per-group figure. `design_stratified("site",
n = 100)` returns 100 rows whether there are 4 sites or 40, and
`design_multistage(..., n = 100)` returns 100 across the selected clusters. The
one deliberate exception is named for what it is: `design_temporal()` takes
`per_interval`.

**`allocation` always splits a total across groups.** `"proportional"` by
stratum size, `"equal"` evenly. Same argument, same meaning, in
`design_stratified()` and `design_multistage()`.

**`na_rm` always decides the same question.** Drop rows whose key is missing, or
raise an error naming the column. Missing keys are never quietly treated as a
group of their own.

**Seeding is local.** `draw(..., seed = 1)` saves `.Random.seed`, seeds, samples,
and restores what was there. Sampling inside a simulation will not shift the
simulation's own stream.

**Two documented exceptions to "same shape out".** `design_bootstrap()` prepends
a `.replicate` column so every replicate comes back in one frame — split it with
`split(out, out$.replicate)`. `design_reservoir()` returns a list when handed a
real stream rather than a data frame.

## Design weights

A sample is only usable for estimation if you know how likely each row was to be
in it. Every design can report that, and `draw()` will attach it:

```r
s <- draw(df, design_stratified("site", n = 100), seed = 1, weights = TRUE)
names(s)[1:2]
#> ".prob" ".weight"

sum(s$spend * s$.weight)    # an unbiased Horvitz-Thompson total
```

`inclusion_prob(data, design)` gives the same probabilities for the whole
population without drawing anything, which is a quick way to sanity-check a
design before you commit to it.

Four designs have **no closed form**, and the package says so instead of guessing:
`design_cluster(balanced = TRUE)`, `design_multistage(allocation =
"proportional")`, `design_weighted(method = "successive")`, and
`design_bootstrap()`. Pass `simulate = TRUE` to estimate them by Monte Carlo.

## Weighted sampling: pick the method deliberately

`design_weighted()` offers three, and they are not interchangeable:

| `method` | Inclusion probabilities | Sample size |
|---|---|---|
| `"successive"` (default) | Not proportional to weight, and no closed form | Fixed |
| `"systematic"` | Exactly `n * p_i` | Fixed |
| `"poisson"` | Exactly `n * p_i` | Random, mean `n` |

The default is what `base::sample(prob=)` does. It is a perfectly good way to
pick rows with a bias toward heavy ones, but the weights govern each sequential
draw rather than the probability of ending up in the sample. How far that lands
from proportional depends on the weights — with one dominant unit the two nearly
coincide, but weights `1:10` at `n = 5` sit up to 35% apart, and treating them as
proportional biases a Horvitz-Thompson total detectably. Use `"systematic"` or
`"poisson"` when the sample is going to be estimated from.

## Credit

Ported from [elkronos/sample_py](https://github.com/elkronos/sample_py).

## License

GPL-3.
