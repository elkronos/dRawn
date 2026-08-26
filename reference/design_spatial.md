# Spatial sampling

Samples rows whose coordinates fall inside a region.

## Usage

``` r
design_spatial(coords, region, n, crs = 4326, na_rm = FALSE)
```

## Arguments

- coords:

  Two column names, `c(x, y)` — longitude then latitude.

- region:

  An `sf` or `sfc` geometry, or a list of them, which is unioned.
  Reprojected to `crs` when it differs.

- n:

  Number of rows to draw from inside the region.

- crs:

  Coordinate reference system of the coordinate columns. Defaults to
  EPSG:4326.

- na_rm:

  Drop rows with missing coordinates instead of raising an error.

## Value

A design object, for use with
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

## Coordinate order

`coords` is `c(x, y)` — longitude first, then latitude — matching
[`sf::st_as_sf()`](https://r-spatial.github.io/sf/reference/st_as_sf.html).
Getting this backwards puts your points somewhere else entirely, so it
is worth checking against a known landmark once.

## Regions that cross the antimeridian

With spherical geometry enabled
([`sf::sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html),
the default), consecutive polygon vertices are joined by the *shortest*
great-circle path. An edge from longitude -179 to +179 therefore spans
the 2 degrees across the antimeridian, not the 358 degrees the
coordinates suggest. A "whole world" rectangle collapses into a
2-degree-wide pole-to-pole strip of 2.8 million km2 – against the
roughly 510 million km2 of the globe – and contains almost nothing.

Ring orientation is not the cause and reversing it does not help: `sf`
normalises winding, so both directions give the same strip. Adding
intermediate vertices does not reliably help either, because the
intermediate edges still bow along geodesics.

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) warns
whenever any edge of `region` spans more than 180 degrees of longitude.
If you hit it, either split the region at the antimeridian into two
polygons, or switch to planar interpretation with
`sf::sf_use_s2(FALSE)`.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md),
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
[`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md),
[`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md),
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)

## Examples

``` r
box <- sf::st_sfc(
  sf::st_polygon(list(cbind(c(0, 10, 10, 0, 0), c(0, 0, 10, 10, 0)))),
  crs = 4326
)
df <- data.frame(id = 1:20,
                 lon = seq(1, 9, length.out = 20),
                 lat = seq(1, 9, length.out = 20))
draw(df, design_spatial(c("lon", "lat"), region = box, n = 5), seed = 1)
#>   id      lon      lat
#> 1  1 1.000000 1.000000
#> 2  2 1.421053 1.421053
#> 3  4 2.263158 2.263158
#> 4  7 3.526316 3.526316
#> 5 13 6.052632 6.052632
```
