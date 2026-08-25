# Reservoir sampling

Draws a uniform sample of fixed size from a stream of unknown length in
a single pass, using Algorithm L.

## Usage

``` r
design_reservoir(n, max_items = Inf)
```

## Arguments

- n:

  Reservoir size. Fewer items come back if the stream is shorter; the
  result is never padded.

- max_items:

  Stop after reading this many items. \`Inf\` reads the whole stream.
  Warns when the cap actually truncates.

## Value

A design object, for use with \[draw()\].

## Details

A data frame is not a stream: its length is already known and it already
fits in memory, so \[draw()\] takes a direct vectorised path for one.
Reach for the streaming path when the data genuinely does not fit — pass
a connection, or a function that returns the next item and \`NULL\` at
end of stream. For those, \[draw()\] returns a list rather than a data
frame.

## See also

\[draw()\]

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md),
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
[`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md),
[`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md),
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)

## Examples

``` r
nrow(draw(data.frame(id = 1:1000), design_reservoir(n = 10), seed = 1))
#> [1] 10

# A real stream: a generator that yields 1..50 then stops
i <- 0
gen <- function() {
  i <<- i + 1
  if (i > 50) NULL else i
}
unlist(draw(gen, design_reservoir(n = 5), seed = 1))
#> [1] 46 14 11  4 32
```
