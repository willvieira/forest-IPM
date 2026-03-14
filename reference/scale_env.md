# Scale environmental values to \[0, 1\]

Scales mean annual temperature and precipitation to the \[0, 1\] range
using the observed ranges stored in the package. These scaled values are
what the IPM vital rate functions expect internally.

## Usage

``` r
scale_env(MAT, MAP)
```

## Arguments

- MAT:

  Numeric. Mean annual temperature in degrees Celsius.

- MAP:

  Numeric. Mean annual precipitation in mm/year.

## Value

A named list with elements `MAT` and `MAP`, each scaled to \[0, 1\].

## Examples

``` r
scaled <- scale_env(MAT = 8, MAP = 1200)
```
