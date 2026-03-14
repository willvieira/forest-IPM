# Unscale environmental values from \[0, 1\] to natural units

Reverses the scaling applied by
[`scale_env()`](https://willvieira.github.io/forestIPM/reference/scale_env.md),
converting scaled values back to degrees Celsius and mm/year.

## Usage

``` r
unscale_env(MAT, MAP)
```

## Arguments

- MAT:

  Numeric. Scaled mean annual temperature in \[0, 1\].

- MAP:

  Numeric. Scaled mean annual precipitation in \[0, 1\].

## Value

A named list with elements `MAT` (degrees Celsius) and `MAP` (mm/year).

## Examples

``` r
scaled   <- scale_env(MAT = 8, MAP = 1200)
original <- unscale_env(scaled$MAT, scaled$MAP)
```
