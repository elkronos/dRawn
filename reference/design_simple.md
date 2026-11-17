# Simple random sampling

Draws `n` rows uniformly at random.

## Usage

``` r
design_simple(n, replace = FALSE, reserve = 0)
```

## Arguments

- n:

  Number of rows to draw. A single non-negative whole number, or a plan
  from
  [`plan_size()`](https://elkronos.github.io/dRawn/reference/plan_size.md).

- replace:

  Sample with replacement? Note that `draw(weights = TRUE)` declines a
  with-replacement design: `.prob` would be the chance of appearing at
  least once, which does not weight a sample holding duplicates. See
  [`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

- reserve:

  Number of reserve units to list after the main sample, to replace
  units found unavailable in the field. A single non-negative whole
  number; `n + reserve` cannot exceed the rows of the frame. Not
  available with `replace = TRUE` or `n = 0`. See "Reserves".

## Value

A design object, for use with
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md).

## Reserves

Fieldwork often finds that some sampled units cannot be reached: a site
with no access, a business that has closed, a plot under water.
`reserve` draws that many more units with the main sample, as a list to
replace them from.
[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md) returns
the `n` units of the main sample followed by the `reserve` units, with
an integer `.reserve` column: 0 for the main sample, and 1, 2, ... for
the reserves, in the order they are used. The list is one random draw of
`n + reserve` units in random order, so the main sample is a simple
random sample of `n`, the one the design without `reserve` draws with
the same seed, and the first `k` units of the list are a simple random
sample of `k` for any `k`.

Visit every unit of the main sample, and for each one found unavailable,
visit the next reserve.
[`use_reserves()`](https://elkronos.github.io/dRawn/reference/use_reserves.md)
takes the list and what the visits found, and returns the sample to
estimate from: the first `n` units of the list found available. They are
a simple random sample of `n` from the units of the frame that are
available, and
[`use_reserves()`](https://elkronos.github.io/dRawn/reference/use_reserves.md)
weights them as one, with probability `n / (N - u)` when `u` listed
units were found unavailable. Replacing units changes what is estimated:
the totals and means are those of the units that can be reached, not of
the whole frame. "What the estimates describe" in
[`use_reserves()`](https://elkronos.github.io/dRawn/reference/use_reserves.md)
has the details.

Until then the list is not a sample to estimate from.
[`inclusion_prob()`](https://elkronos.github.io/dRawn/reference/inclusion_prob.md)
and the `.prob` column give every row the main sample's probability,
`n / N`, and
[`ht_total()`](https://elkronos.github.io/dRawn/reference/ht_total.md),
[`ht_mean()`](https://elkronos.github.io/dRawn/reference/ht_mean.md),
[`as_svydesign()`](https://elkronos.github.io/dRawn/reference/as_svydesign.md)
and
[`sample_summary()`](https://elkronos.github.io/dRawn/reference/sample_summary.md)
refuse the list.

## See also

[`draw()`](https://elkronos.github.io/dRawn/reference/draw.md),
[`use_reserves()`](https://elkronos.github.io/dRawn/reference/use_reserves.md)

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_certainty()`](https://elkronos.github.io/dRawn/reference/design_certainty.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md),
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md),
[`design_spread()`](https://elkronos.github.io/dRawn/reference/design_spread.md),
[`design_stratified()`](https://elkronos.github.io/dRawn/reference/design_stratified.md),
[`design_systematic()`](https://elkronos.github.io/dRawn/reference/design_systematic.md),
[`design_temporal()`](https://elkronos.github.io/dRawn/reference/design_temporal.md),
[`design_weighted()`](https://elkronos.github.io/dRawn/reference/design_weighted.md)

## Examples

``` r
df <- data.frame(id = 1:100)
draw(df, design_simple(n = 10), seed = 1)
#>    id
#> 1  68
#> 2  39
#> 3   1
#> 4  34
#> 5  87
#> 6  43
#> 7  14
#> 8  82
#> 9  59
#> 10 51
draw(df, design_simple(n = 150, replace = TRUE), seed = 1)
#>      id
#> 1    68
#> 2    39
#> 3     1
#> 4    34
#> 5    87
#> 6    43
#> 7    14
#> 8    82
#> 9    59
#> 10   51
#> 11   97
#> 12   85
#> 13   21
#> 14   54
#> 15   74
#> 16    7
#> 17   73
#> 18   79
#> 19   85
#> 20   37
#> 21   89
#> 22   37
#> 23   34
#> 24   89
#> 25   44
#> 26   79
#> 27   33
#> 28   84
#> 29   35
#> 30   70
#> 31   74
#> 32   42
#> 33   38
#> 34   20
#> 35   28
#> 36   20
#> 37   44
#> 38   87
#> 39   70
#> 40   40
#> 41   44
#> 42   25
#> 43   70
#> 44   39
#> 45   51
#> 46   42
#> 47    6
#> 48   24
#> 49   32
#> 50   14
#> 51    2
#> 52   45
#> 53   18
#> 54   22
#> 55   78
#> 56   65
#> 57   70
#> 58   87
#> 59   70
#> 60   75
#> 61   81
#> 62  100
#> 63   13
#> 64   40
#> 65   89
#> 66   48
#> 67   89
#> 68   23
#> 69   84
#> 70   29
#> 71   13
#> 72   22
#> 73   93
#> 74   28
#> 75   48
#> 76   33
#> 77   45
#> 78   21
#> 79   31
#> 80   17
#> 81   73
#> 82   87
#> 83   83
#> 84   90
#> 85   48
#> 86   64
#> 87   94
#> 88   96
#> 89   60
#> 90   51
#> 91   93
#> 92   34
#> 93   10
#> 94    1
#> 95   43
#> 96   59
#> 97   26
#> 98   15
#> 99   58
#> 100  29
#> 101  24
#> 102  42
#> 103  48
#> 104  76
#> 105  39
#> 106  24
#> 107  53
#> 108  92
#> 109  86
#> 110  40
#> 111  97
#> 112  83
#> 113  90
#> 114  35
#> 115  43
#> 116   1
#> 117  29
#> 118  78
#> 119  22
#> 120  70
#> 121  28
#> 122  37
#> 123  61
#> 124  46
#> 125  67
#> 126  86
#> 127  99
#> 128  71
#> 129  99
#> 130  51
#> 131  44
#> 132  49
#> 133  60
#> 134  56
#> 135  49
#> 136  50
#> 137  91
#> 138   7
#> 139  20
#> 140  24
#> 141  51
#> 142  53
#> 143  16
#> 144  83
#> 145 100
#> 146   2
#> 147  48
#> 148  65
#> 149  44
#> 150  77

# Ten units and three reserves, used in the order of `.reserve`
draw(df, design_simple(n = 10, reserve = 3), seed = 1)
#>    .reserve id
#> 1         0 68
#> 2         0 39
#> 3         0  1
#> 4         0 34
#> 5         0 87
#> 6         0 43
#> 7         0 14
#> 8         0 82
#> 9         0 59
#> 10        0 51
#> 11        1 95
#> 12        2 23
#> 13        3 61
```
