# Draw method for a design class

The generic
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) calls to
make a sample from a design. Each design in the package has a method,
and a design of your own can have one too.

## Usage

``` r
draw_design(design, data)
```

## Arguments

- design:

  A design object.

- data:

  The data to sample from.

## Value

The sampled rows. The cluster and multistage methods attach the codes of
the clusters they selected as the attribute `drawn_selection`, which
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) removes
or, for a sample weighted conditional on its clusters, keeps.

## Writing a design of your own

A design is a list with a class vector ending in `"drawn_design"`, so
that
[`is_design()`](https://elkronos.github.io/dRawn/reference/is_design.md)
accepts it, such as
`structure(list(n = 5), class = c("first_rows", "drawn_design"))`. Its
method takes `design` and `data` and returns rows of `data`, with every
column of `data`, as `data[idx, , drop = FALSE]` gives them. It draws
its random numbers from the current stream;
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) applies
the seed around it. A package registers the method with
`@exportS3Method drawn::draw_design`, and a script with
[`.S3method()`](https://rdrr.io/r/base/S3method.html).

That is enough for
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) without
weights, and for
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
and
[`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)
with `simulate = TRUE`, which count how often each row is drawn and find
the rows through a column they add to `data`. The exact probabilities,
`draw(weights = TRUE)` and the estimators need the design's inclusion
probabilities, which the package computes for its own designs only, so
they refuse a design of your own.

## Examples

``` r
first_rows <- structure(list(n = 3), class = c("first_rows", "drawn_design"))
.S3method("draw_design", "first_rows", function(design, data) {
  data[seq_len(design$n), , drop = FALSE]
})

df <- data.frame(id = 1:10)
draw(df, first_rows)
#>   id
#> 1  1
#> 2  2
#> 3  3
inclusion_prob(df, first_rows, simulate = TRUE, R = 20, seed = 1)
#>  [1] 1 1 1 0 0 0 0 0 0 0
```
