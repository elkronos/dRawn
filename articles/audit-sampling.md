# Walkthrough: an audit sample

``` r

library(drawn)
```

An auditor has 2,000 invoices and time to examine about 150. The
questions are the usual ones: how much is the ledger misstated in total,
with an upper bound the auditor can sign off against; what share of
invoices contain an error; and is any site worse than the others? This
walkthrough builds the sample in the shape auditors use — examine every
large item, sample the rest in proportion to size — and shows what each
choice buys.

## The ledger

Book values are known for every invoice. Misstatements are only found by
examining one. Here the ledger holds them too, so the answers can be
checked, but nothing below uses them except through the sample.

``` r

set.seed(31)
N <- 2000
ledger <- data.frame(
  invoice = sprintf("INV%05d", 1:N),
  site    = sample(c("north", "south", "east", "west"), N, replace = TRUE,
                   prob = c(0.35, 0.3, 0.2, 0.15)),
  book    = round(rlnorm(N, meanlog = 7, sdlog = 1.2), 2)
)
# Errors are more common at the west site; each misstates a share of the
# invoice's book value
error_rate <- c(north = 0.05, south = 0.05, east = 0.04, west = 0.15)
has_error <- runif(N) < error_rate[ledger$site]
ledger$misstatement <- ifelse(has_error,
                              round(ledger$book * runif(N, 0.05, 0.6), 2), 0)

c(book_total = sum(ledger$book), misstatement = sum(ledger$misstatement),
  error_rate = mean(ledger$misstatement > 0))
#>   book_total misstatement   error_rate 
#>  4393229.220    83884.920        0.066
```

A handful of invoices carry a large share of the book value — the usual
shape of a ledger, and the reason simple random sampling is a poor fit:

``` r

top <- sort(ledger$book, decreasing = TRUE)
round(c(top_1pct = sum(top[1:20]), top_5pct = sum(top[1:100])) / sum(top), 2)
#> top_1pct top_5pct 
#>     0.11     0.31
```

## Four ways to spend 150 examinations

- **Simple random**: every invoice equally likely.
- **Stratified by value band**, with Neyman allocation on book value:
  more examinations where values vary most.
- **Certainty plus PPS**: every invoice over 10,000 examined, and the
  rest drawn with probability proportional to book value — the design
  behind monetary-unit sampling.
- **Certainty plus spread**: the same certainty stratum, and the rest
  spread evenly along book value with equal probabilities — a fine
  stratification on size without choosing the band boundaries.

``` r

ledger$band <- cut(ledger$book, c(0, 500, 2000, 10000, Inf),
                   labels = c("<500", "500-2k", "2k-10k", "10k+"))
n_certain <- sum(ledger$book >= 10000)
n_certain
#> [1] 65

designs <- list(
  simple    = design_simple(n = 150),
  banded    = design_stratified("band", n = 150, allocation = "neyman",
                                allocation_by = "book", min_per_stratum = 2),
  cert_pps  = design_certainty("book", 10000,
                               design_weighted("book", n = 150 - n_certain,
                                               method = "systematic")),
  cert_spread = design_certainty("book", 10000,
                                 design_spread("book", n = 150 - n_certain))
)
```

All four examine 150 invoices. Which is most precise for the
misstatement total depends on how misstatements relate to size — here
they scale with book value — and simulation settles it:

``` r

spread_of <- function(design, R = 150) {
  sd(vapply(seq_len(R), function(i) {
    ht_total(draw(ledger, design, seed = i, weights = TRUE), "misstatement",
             variance = "none")$total
  }, numeric(1)))
}
sds <- vapply(designs, spread_of, numeric(1))
round(sds)
#>      simple      banded    cert_pps cert_spread 
#>       41694       32330       32327       40931
```

The designs that let book value decide how likely an invoice is to be
examined — Neyman bands, and certainty plus PPS — beat simple random
sampling by 22% in standard error. That is a real gain but not a
dramatic one, and the reason is worth seeing: only 7% of invoices are
misstated at all, so *which* invoices happen to contain an error
dominates the uncertainty whatever the design. Spreading the sample
evenly along book value with equal probabilities does not help here,
because it does not put more examinations where misstatements are
largest.

The certainty-plus-PPS design is the one auditors know, and the one used
below.

## Draw and check

``` r

audit_design <- designs$cert_pps
audit_design
#> <sampling design: certainty>
#>   above      "book"
#>   threshold  10000
#>   rest       design_weighted(weights = "book", n = 85, method = "systematic")
#>   na_rm      FALSE

s <- draw(ledger, audit_design, seed = 11, weights = TRUE)
sample_summary(s)
#> Sample of 150 from 2,000  (certainty design)
#>   sampling fraction  0.075
#>   design weights     1 to 124.4   (cv 1.44)
table(certain = s$.prob == 1)
#> certain
#> FALSE  TRUE 
#>    85    65
```

## Estimate

Total misstatement, with a one-sided 95% upper bound — which is the
upper end of a two-sided 90% interval:

``` r

m <- ht_total(s, "misstatement", level = 0.90)
m
#> Horvitz-Thompson total  (certainty design, n = 150)
#>   estimate 71,954.44
#>   se       33,290.69  (Deville approximation)
#>   90% CI  16,585.47 to 127,323.4  (t, 84 df)
#>   deff     0.68  (better than simple random sampling)
#> 
#>   Deville's approximation for unequal-probability sampling: systematic
#>   PPS has no closed-form joint inclusion probabilities.
c(upper_bound_95 = m$ci[2], truth = sum(ledger$misstatement))
#> upper_bound_95          truth 
#>      127323.40       83884.92
```

The variance is Deville’s approximation for the PPS part; the certainty
invoices were all examined, so they add nothing to it.

The rate of invoices in error is the mean of a yes/no column. Its
interval is formed on the logit scale, so it stays inside 0 to 1 and is
lopsided the way a small proportion’s uncertainty is:

``` r

ht_mean(s, s$misstatement > 0)
#> Hajek mean  (certainty design, n = 150)
#>   estimate 0.0523489
#>   se       0.02929237  (Deville approximation)
#>   95% CI  0.01678616 to 0.1516344  (logit, t, 84 df)
#>   deff     2.79  (worse than simple random sampling)
#> 
#>   Deville's approximation for unequal-probability sampling: systematic
#>   PPS has no closed-form joint inclusion probabilities.
```

The design effect is well above 1: a design built to estimate a total in
money is a poor one for counting errors, which do not scale with invoice
size. If the error rate mattered as much as the total, stratifying by
value band would serve both better.

And by site, from the same sample:

``` r

ht_mean(s, s$misstatement > 0, by = "site")
#> Hajek mean by domain  (certainty design, 95% CI, logit, t, 84 df)
#>   site  n     mean        se  ci_lower ci_upper  method
#>   east 34 0.031268 0.0282906 0.0050127 0.171355 deville
#>  north 58 0.001608 0.0004325 0.0009417 0.002745 deville
#>  south 37 0.063771 0.0503184 0.0125856 0.266867 deville
#>   west 21 0.166276 0.1444400 0.0245010 0.612949 deville
#> 
#>   Deville's approximation for unequal-probability sampling: systematic
#>   PPS has no closed-form joint inclusion probabilities.
```

West stands out. Because site was not a stratum, the number of west
invoices examined was itself random, and the standard errors account for
that.

## How far to trust the bound

Misstatements are rare and skewed: most examined invoices contribute
zero, a few contribute a lot. The normal and t approximations behind the
interval assume the estimate’s sampling distribution is roughly
symmetric, and that is least true exactly here. So check it — draw the
design many times and count how often the upper bound actually sits
above the true total:

``` r

truth <- sum(ledger$misstatement)
covered <- vapply(1:300, function(i) {
  r <- ht_total(draw(ledger, audit_design, seed = 1000 + i, weights = TRUE),
                "misstatement", level = 0.90)
  r$ci[2] >= truth
}, logical(1))
mean(covered)
#> [1] 0.8666667
```

87% — short of the 95% the bound claims. Most samples contain only a few
misstated invoices, so the estimate is usually a bit low and
occasionally far too high, and a symmetric interval around it is placed
too low exactly when it matters. That is not a defect of the variance
estimator; it is the normal approximation meeting a rare, skewed
quantity.

Audit practice has bounds built for this case — the Stringer bound and
Bayesian methods, which the
[`jfa`](https://cran.r-project.org/package=jfa) package implements. Two
things help within a design-based analysis: a larger sample, which makes
the estimate’s distribution more symmetric, and a lower one-sided level.
`drawn`’s part is the design and a design-based estimate; checking a
bound this way before signing off against it is cheap, and this is the
case that shows why it is worth doing.

## Hand off

For anything further — a ratio of misstatement to book value by site, a
regression of taint on invoice characteristics — the sample goes to
`survey` with its design intact:

``` r

des <- as_svydesign(s)
survey::svyratio(~misstatement, ~book, des)
#> Ratio estimator: svyratio.survey.design2(~misstatement, ~book, des)
#> Ratios=
#>                    book
#> misstatement 0.01637848
#> SEs=
#>                     book
#> misstatement 0.007577524
```
