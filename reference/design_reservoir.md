# Reservoir sampling

Draws a uniform sample of fixed size from a stream of unknown length in
a single pass, using Algorithm L (Li 1994), which skips ahead
geometrically instead of drawing a random number for every item (Vitter
1985).

## Usage

``` r
design_reservoir(n, max_items = Inf, batches = FALSE)
```

## Arguments

- n:

  Reservoir size, or a plan from
  [`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md).
  Fewer items come back if the stream is shorter; the result is never
  padded.

- max_items:

  Stop after this many items, or rows with `batches = TRUE`. Items past
  the cap are never sampled, so
  [`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
  gives rows past it `0` and gives the rows before it `n / max_items`
  rather than `n / N`. A stream is read one item past the cap; see "What
  is read". `Inf` reads the whole stream. Warns when the cap truncates.

- batches:

  If `TRUE`, `data` is a function that returns a data frame of rows on
  each call and `NULL` at the end, and the reservoir samples rows rather
  than whole data frames.
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) then
  returns a data frame. A data frame given as `data` is sampled by rows
  either way; a list or a connection is refused.

## Value

A design object, for use with
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

## Details

A data frame or a list is not a stream: its length is already known and
it already fits in memory, so
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) picks `n`
positions in it at random and returns them in their original order.
Every element of a list is an item, `NULL` included, and the elements
keep their names. The streaming path is for data that does not fit in
memory. Pass it a connection, or a function that returns the next item
and `NULL` at end of stream. It holds at most `n` items at a time, plus
the block it is reading, and returns them in the order they arrived. For
a list, a connection or a function,
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) returns a
list rather than a data frame, except with `batches = TRUE`.

The two paths use different algorithms, so for one seed a stream gives a
different sample from a data frame or list holding the same items. Each
path on its own gives the same sample for the same seed.

## Items in a stream

Each line of a connection is one item, a character string as
[`readLines()`](https://rdrr.io/r/base/readLines.html) returns it, and a
header line is an item like any other. To leave a CSV header out of the
sample, open the connection, read the header with
`readLines(con, n = 1)` and then draw from the connection, which is read
from where the header left off.
`read.csv(text = c(header, unlist(lines)))` turns the sampled lines into
a data frame (see the examples).
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) reads a
connection in blocks of 10,000 lines, which gives the same sample for a
seed as reading it line by line.

A function is called once for each item, and whatever it returns is that
item, a data frame included, so a function that returns blocks of rows
gives a sample of whole blocks. With `batches = TRUE` the rows are the
items instead: each call returns a data frame of rows, every batch with
the same columns in the same order, and
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) returns
`n` of those rows as one data frame, bound with
[`rbind()`](https://rdrr.io/r/base/cbind.html) in the order they
arrived. [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)
compares only the column names. A column whose type differs between
batches, such as one that
[`read.csv()`](https://rdrr.io/r/utils/read.table.html) reads as logical
in a batch where it is empty, takes the type
[`rbind()`](https://rdrr.io/r/base/cbind.html) gives the sampled rows.
Cutting the same rows into batches of other sizes gives the same sample
for a seed. Between batches it keeps at most `max(2 * n, 1024)` rows,
some of them already dropped from the sample, besides the batch being
read.

## What is read

A connection that is not open is opened, read to its end (or to one line
past a finite `max_items`) and closed, as
[`read.csv()`](https://rdrr.io/r/utils/read.table.html) does, so it
cannot be drawn from twice. A connection that is already open is read
from where it stands and left open.

With a finite `max_items`,
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) reads one
item past the cap, because only that read tells a stream of exactly
`max_items` items from a longer one. On an open connection in text mode
the line goes back with
[`pushBack()`](https://rdrr.io/r/base/pushBack.html), so the next read
starts with it. A connection opened in binary mode (`"rb"`) cannot take
a line back, and loses it. A function is called once past the cap and
what it returns is dropped; with `batches = TRUE`, so are the rows of a
batch that run past the cap.

## References

Vitter, J. S. (1985). Random sampling with a reservoir. *ACM
Transactions on Mathematical Software*, 11, 37–57.

Li, K.-H. (1994). Reservoir-sampling algorithms of time complexity
O(n(1 + log(N/n))). *ACM Transactions on Mathematical Software*, 20,
481–493.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md),
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md),
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
#> [1]  4 11 14 32 46

# A CSV file read line by line, with its header left out of the sample
path <- tempfile(fileext = ".csv")
write.csv(data.frame(id = 1:200, value = (1:200) %% 7), path,
          row.names = FALSE)
con <- file(path, "r")
header <- readLines(con, n = 1)
lines <- draw(con, design_reservoir(n = 5), seed = 1)
close(con)
read.csv(text = c(header, unlist(lines)))
#>    id value
#> 1  11     4
#> 2  32     4
#> 3  99     1
#> 4 105     0
#> 5 160     6

# The same file parsed 50 rows at a time, sampling rows: the same five
# rows as above for the same seed
con <- file(path, "r")
header <- readLines(con, n = 1)
next_rows <- function() {
  lines <- readLines(con, n = 50)
  if (length(lines) == 0) NULL else read.csv(text = c(header, lines))
}
draw(next_rows, design_reservoir(n = 5, batches = TRUE), seed = 1)
#>    id value
#> 1  11     4
#> 2  32     4
#> 3  99     1
#> 4 105     0
#> 5 160     6
close(con)
unlink(path)
```
