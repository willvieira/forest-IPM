# Project population dynamics over time

Project population dynamics over time

## Usage

``` r
project(mod, pars, stand, env, ctrl)
```

## Arguments

- mod:

  An `ipm_spModel` object.

- pars:

  An `ipm_parameters` object.

- stand:

  An `ipm_stand` object. Initial stand state (t = 0).

- env:

  An `ipm_env` object. Climate drivers (may be time-varying).

- ctrl:

  An `ipm_control` object. Required.

## Value

An object of S3 class `"ipm_projection"` with fields: `$species`,
`$years`, `$lambda`, `$stand_series` (list of `ipm_dist_snapshot`
objects, each holding the continuous size distribution per species at
that timestep), `$summary`.

## Examples

``` r
df <- data.frame(size_mm = c(150, 200, 350),
                 species_id = "ABIBAL",
                 plot_size  = 1000)
s    <- stand(df)
mod  <- species_model(s)
pars <- parameters(mod, draw = "mean")
env  <- env_condition(MAT = 8, MAP = 1200)
ctrl <- control(years = 5, compute_lambda = TRUE, progress = FALSE)
proj <- project(mod, pars, s, env, ctrl)
print(proj)
#> <ipm_projection>  1 species | 5 years | store_every=1 | draw=mean
```
