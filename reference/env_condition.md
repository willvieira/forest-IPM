# Create an ipm_env object specifying climate drivers

Create an ipm_env object specifying climate drivers

## Usage

``` r
env_condition(MAT, MAP)
```

## Arguments

- MAT:

  Numeric scalar or `function(t)`. Mean Annual Temperature in Celsius.
  Use a function for time-varying climate scenarios.

- MAP:

  Numeric scalar or `function(t)`. Mean Annual Precipitation in mm/year.

## Value

An object of S3 class `"ipm_env"` with fields `$MAT` and `$MAP`
(unscaled; internal scaling applied by engines).

## Examples

``` r
env <- env_condition(MAT = 8, MAP = 1200)
print(env)
#> <ipm_env>  MAT=8.0°C  MAP=1200 mm/yr
```
