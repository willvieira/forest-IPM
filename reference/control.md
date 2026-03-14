# Configure IPM projection settings

Configure IPM projection settings

## Usage

``` r
control(
  years = 100,
  delta_time = 1,
  store_every = 1,
  bin_width = 1,
  compute_lambda = FALSE,
  progress = TRUE
)
```

## Arguments

- years:

  Positive integer. Number of simulation timesteps. Default 100.

- delta_time:

  Positive numeric. Duration of each timestep in years. Default 1.

- store_every:

  Positive integer. Store stand state every N timesteps. Default 1.

- bin_width:

  Positive integer. Bin width for IPM kernel discretization. Default 1.

- compute_lambda:

  Logical. Whether to compute the asymptotic lambda at each timestep via
  eigendecomposition. Set to FALSE to skip (faster projections when only
  population structure is needed). Default FALSE

- progress:

  Logical. Whether to display a progress bar during projection. Default
  TRUE.

## Value

An object of S3 class `"ipm_control"`.

## Examples

``` r
ctrl <- control(years = 10, compute_lambda = TRUE, progress = FALSE)
print(ctrl)
#> <ipm_control>  10 years | dt=1.0 | store_every=1 | bin_width=1 | lambda=yes | progress=no
```
