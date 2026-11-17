# Sampling within a region

Draws a simple random sample of the rows whose coordinates fall inside a
region. It decides *where* sampling happens, not how evenly it covers
the ground: for a sample spread evenly across space — usually far more
precise when what you measure varies smoothly — use
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md),
on the frame restricted to the region if need be.

## Usage

``` r
design_spatial(coords, region, n, crs = 4326, na_rm = FALSE)
```

## Arguments

- coords:

  Two column names, `c(x, y)` with longitude first, or the name of one
  POINT geometry column. See "Coordinate order".

- region:

  An `sf` or `sfc` geometry of POLYGON or MULTIPOLYGON features, or a
  list of them. A row inside any of them is inside the region. Each is
  reprojected to the coordinates' CRS when its own differs.

- n:

  Number of rows to draw from inside the region, or a plan from
  [`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md).

- crs:

  Coordinate reference system of the coordinate columns. Defaults to
  EPSG:4326. A geometry column named in `coords` keeps its own CRS when
  it has one.

- na_rm:

  Drop rows with missing coordinates instead of raising an error.

## Value

A design object, for use with
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

## Coordinate order

`coords` is `c(x, y)` — longitude first, then latitude — matching
[`sf::st_as_sf()`](https://r-spatial.github.io/sf/reference/st_as_sf.html).
Getting this backwards puts your points somewhere else entirely, so it
is worth checking against a known landmark once. When the coordinates
are longitude and latitude, a latitude outside -90 to 90 or a longitude
outside -180 to 360 is refused, since it means the columns are the wrong
way round or hold projected coordinates.

A POINT geometry column, such as the geometry of an `sf` frame, can
stand in for the two columns: name it alone, as in
`coords = "geometry"`, and each row's position is read with
[`sf::st_coordinates()`](https://r-spatial.github.io/sf/reference/st_coordinates.html).
The column's own CRS is used when it has one, and `region` is
reprojected to it.

## Edges of a longitude/latitude region

With spherical geometry
([`sf::sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html),
the default), each edge of a longitude/latitude `region` follows the
great circle between its two corners, the shortest path on the globe. On
a flat longitude/latitude map that path bows towards the nearer pole
wherever the edge runs east-west away from the equator; edges along a
meridian or the equator do not bow. A wide box is therefore a different
shape from the rectangle it looks like on such a map. The top edge of a
box from -50 to 50 degrees of longitude and -30 to 30 of latitude
reaches latitude 41.9 at longitude 0, so a row at (0, 41) is inside it,
and its bottom edge reaches -41.9. A box from -5 to 5 degrees of
longitude and 40 to 50 of latitude takes in a row at (0, 50.05) and
leaves out one at (0, 40.05).

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md),
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md),
[`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)
and
[plot()](https://elkronos.github.io/dRawn/reference/plot.drawn_design.md)
warn, with class `drawn_geodesic_edges`, when an edge of `region` passes
more than 0.1 degrees (about 11 km) from the straight line between its
corners on a longitude/latitude map. To make the edges follow those
straight lines, add vertices along them in degree space before you build
the design:

    flat <- sf::st_set_crs(region, NA)    # coordinates as plain numbers
    flat <- sf::st_segmentize(flat, 0.1)  # a vertex at least every 0.1 degrees
    region <- sf::st_set_crs(flat, 4326)

An edge 0.1 degrees long bows by about a metre at most. The CRS has to
come off for the middle step: on a region that keeps it,
[`sf::st_segmentize()`](https://r-spatial.github.io/sf/reference/geos_unary.html)
needs the lwgeom package and places the new vertices along the great
circles, under either s2 setting, and the bow stays. A region larger
than a hemisphere needs splitting first; see the next section.

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
vertices along the edges in degree space, as in the previous section,
does not rescue a region this large either. s2 takes the smaller of the
two areas a ring divides the globe into as its inside, so a ring around
more than half the globe encloses the rest of the globe instead.

The short way is the right reading for a region meant to cross the
antimeridian, such as a box from longitude 179 to -179.
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md),
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md),
[`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)
and
[plot()](https://elkronos.github.io/dRawn/reference/plot.drawn_design.md)
warn, with class `drawn_antimeridian`, whenever an edge of `region`
joins two longitudes more than 180 degrees apart, and the warning gives
both readings. An edge from -180 to 180 is one of them: both ends lie on
the same meridian, so under spherical geometry the edge spans no
longitude, and a box from -180 to 180 encloses nothing. An edge along a
pole, where every longitude meets, does not warn, so a cap around a pole
closed along it does not either. If you meant the short way, `region` is
what you drew. If you meant the long way, split the region into parts
smaller than a hemisphere, add vertices to each as in the previous
section, and pass the parts as a list; or switch to planar
interpretation with `sf::sf_use_s2(FALSE)`. Splitting works under the
default setting. Switching has to happen before
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) and stay
in force for every estimate from the sample; see the next section.

## The s2 setting at draw time and after

Which rows fall inside `region` depends on
[`sf::sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html).
With spherical geometry the edges of a longitude/latitude region are
great-circle arcs, and without it they are straight lines in longitude
and latitude, so rows near an edge can be inside under one setting and
outside under the other. Use the same setting when you draw a sample and
whenever you analyse it.

[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md),
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
and
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
find the region's rows again each time they run, because the variance,
the design effect and `survey`'s population count all depend on them.
They compare the result with the inclusion probabilities stored at draw
time, or with the `.prob` column of a sample saved by drawn 0.2.0, and
refuse the sample when the two differ, naming
[`sf::sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html).
The comparison catches any change in which rows are inside, including
one that comes with another version of sf. A region that holds fewer
than `n` rows under the current setting leaves nothing to compare, and
the sample is refused the same way. The point estimate would not have
moved, since it uses the stored `.prob`, but the standard error would.

`sf_use_s2()` is a session setting, and a new R session starts with it
`TRUE` unless the `_SF_USE_S2` environment variable says otherwise. A
sample drawn under `sf::sf_use_s2(FALSE)` and saved with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) needs
`sf::sf_use_s2(FALSE)` again after
[`readRDS()`](https://rdrr.io/r/base/readRDS.html).
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
and
[`joint_prob()`](https://elkronos.github.io/dRawn/reference/joint_prob.md)
on the frame follow the current setting too.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md),
[`design_simple()`](https://elkronos.github.io/dRawn/reference/design_simple.md),
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md),
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
