# Is this a sampling design?

Is this a sampling design?

## Usage

``` r
is_design(x)
```

## Arguments

- x:

  An object.

## Value

`TRUE` if `x` was built by one of the `design_*()` constructors.

## Examples

``` r
is_design(design_simple(n = 10))
#> [1] TRUE
is_design(mtcars)
#> [1] FALSE
```
