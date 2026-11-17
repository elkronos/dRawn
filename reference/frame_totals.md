# Population totals from the frame behind a sample

Totals the columns of the frame a sample was drawn from, as stored with
the sample, in the form
[`survey::calibrate()`](https://rdrr.io/pkg/survey/man/calibrate.html)
takes as `population`. A design handed over by
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
can then be calibrated, or post-stratified, to the frame without
building the totals by hand.

## Usage

``` r
frame_totals(sample, formula, table = FALSE)
```

## Arguments

- sample:

  A data frame returned by
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) with
  `weights = TRUE`.

- formula:

  A one-sided formula of frame columns, such as `~ w + site`, as
  [`survey::calibrate()`](https://rdrr.io/pkg/survey/man/calibrate.html)
  takes it. A name in it that is not a column of the frame, such as
  `cutoff` in `~ I(w > cutoff)`, is read from the environment of
  `formula`, as
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html)
  reads it. A term with no frame column in it, as in `~ cutoff`, has no
  value for each row and is refused, and so is one that gives a single
  value for the frame, such as `I(sum(w))`.

- table:

  `FALSE`, the default, for the column totals of the model matrix of
  `formula`, which
  [`survey::calibrate()`](https://rdrr.io/pkg/survey/man/calibrate.html)
  takes. `TRUE` for the count of rows in each cell of the variables in
  `formula`, as a table, which
  [`survey::postStratify()`](https://rdrr.io/pkg/survey/man/postStratify.html)
  and [`survey::rake()`](https://rdrr.io/pkg/survey/man/rake.html) take.

## Value

With `table = FALSE`, a named numeric vector: the totals of the columns
of `model.matrix(formula)` over the rows the design could reach,
`(Intercept)` being the number of those rows. With `table = TRUE`, a
table of counts of those rows from
[`stats::xtabs()`](https://rdrr.io/r/stats/xtabs.html).

## Details

The totals run over the rows the design could have drawn, those with an
inclusion probability above 0 when the sample was drawn: the rows inside
the region of
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md)
or the window of
[`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md),
those
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md)
read before `max_items`, and not the rows `na_rm` dropped. Those are the
rows the sample's weights stand for. For a sample weighted conditional
on the clusters selected (see "Design weights" in
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md)), the
rows of the clusters left out count too, since another selection could
have drawn them. A missing value of a variable in `formula` in any of
those rows is refused, since its total is unknown.

A value of a column in `formula` that no sampled row has cannot be
calibrated to.
[`survey::calibrate()`](https://rdrr.io/pkg/survey/man/calibrate.html)
stops on such a factor level, and leaves out the total of such a
character value with a warning that it reordered the totals. Calibration
itself is left to `survey`, which gives the calibrated estimates their
standard errors.

## See also

[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)

## Examples

``` r
pop <- data.frame(
  id = 1:400,
  site = rep(c("a", "b", "c", "d"), times = c(200, 100, 60, 40)),
  size = 20 + (1:400 * 53) %% 300
)
pop$spend <- 3 * pop$size + (1:400 * 37) %% 90
s <- draw(pop, design_simple(n = 60), seed = 1, weights = TRUE)
des <- as_svydesign(s)

# Calibrate to the frame's row count and total size
totals <- frame_totals(s, ~size)
totals
#> (Intercept)        size 
#>         400       67300 
cal <- survey::calibrate(des, ~size, population = totals)
survey::svytotal(~spend, cal)
#>        total     SE
#> spend 218918 1238.9
ht_total(s, "spend")
#> Horvitz-Thompson total  (simple design, n = 60)
#>   estimate 219,173.3
#>   se       12,239.06  (analytic)
#>   95% CI  194,683 to 243,663.6  (t, 59 df)
#>   deff     1  (about the same as simple random sampling)

# Post-stratify on site
counts <- frame_totals(s, ~site, table = TRUE)
post <- survey::postStratify(des, ~site, counts)
survey::svytotal(~spend, post)
#>        total    SE
#> spend 219786 12124
```
