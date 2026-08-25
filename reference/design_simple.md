# Simple random sampling

Draws \`n\` rows uniformly at random.

## Usage

``` r
design_simple(n, replace = FALSE)
```

## Arguments

- n:

  Number of rows to draw. A single non-negative whole number.

- replace:

  Sample with replacement?

## Value

A design object, for use with \[draw()\].

## See also

\[draw()\]

Other designs:
[`design_bootstrap()`](https://elkronos.github.io/dRawn/reference/design_bootstrap.md),
[`design_cluster()`](https://elkronos.github.io/dRawn/reference/design_cluster.md),
[`design_multistage()`](https://elkronos.github.io/dRawn/reference/design_multistage.md),
[`design_reservoir()`](https://elkronos.github.io/dRawn/reference/design_reservoir.md),
[`design_spatial()`](https://elkronos.github.io/dRawn/reference/design_spatial.md),
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
```
