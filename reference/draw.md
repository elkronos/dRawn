# Draw a sample

Applies a sampling design to data.

## Usage

``` r
draw(data, design, seed = NULL, weights = FALSE)
```

## Arguments

- data:

  A data frame. \[design_reservoir()\] also accepts a list, a
  connection, or a zero-argument generator function.

- design:

  A design built by one of the \[design_simple()\] family.

- seed:

  Optional seed, applied only for this draw.

- weights:

  Attach \`.prob\` and \`.weight\` columns. See "Design weights".

## Value

The sampled rows. See "What you get back".

## What you get back

A data frame with the same class and the same columns, in the same
order, as \`data\`. Two designs qualify that:

\* \[design_bootstrap()\] prepends a \`.replicate\` column identifying
which replicate each row belongs to, so all replicates come back in one
frame. Split them with \`split(out, out\$.replicate)\` if you need a
list. \* \[design_reservoir()\] returns a list when \`data\` is a stream
rather than a data frame, because there is nothing to make a data frame
out of.

## Design weights

\`weights = TRUE\` prepends two columns: \`.prob\`, the probability that
the row was included in the sample, and \`.weight\`, its reciprocal —
the number of population units the row stands for. Those are what an
unbiased estimate needs: a Horvitz-Thompson total is \`sum(y \*
.weight)\`.

They come from the design applied to the population, not from the drawn
sample, and four designs have no closed form for them.
\[inclusion_prob()\] documents which, why, and what to do instead.

## Seeding

\`seed\` is saved, applied, and unwound: \`.Random.seed\` is restored on
exit, so drawing a sample inside a simulation does not shift the
simulation's own random number stream. Leave it \`NULL\` to draw from
the current stream.

## See also

\[designs\] for the full list of constructors, \[inclusion_prob()\] for
the probabilities themselves, and \[ht_total()\] to estimate a
population total with a standard error.

## Examples

``` r
df <- data.frame(id = 1:100, site = rep(letters[1:4], each = 25))

draw(df, design_simple(n = 10), seed = 1)
#>    id site
#> 1  68    c
#> 2  39    b
#> 3   1    a
#> 4  34    b
#> 5  87    d
#> 6  43    b
#> 7  14    a
#> 8  82    d
#> 9  59    c
#> 10 51    c

# The same design, reused
by_site <- design_stratified(strata = "site", n = 12)
table(draw(df, by_site, seed = 1)$site)
#> 
#> a b c d 
#> 3 3 3 3 
table(draw(df[1:60, ], by_site, seed = 1)$site)
#> 
#> a b c 
#> 5 5 2 

# With design weights, ready for estimation
df$spend <- seq_len(100)
s <- draw(df, by_site, seed = 1, weights = TRUE)
sum(s$spend * s$.weight)   # estimates sum(df$spend) = 5050
#> [1] 4966.667
```
