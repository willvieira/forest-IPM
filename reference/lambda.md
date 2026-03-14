# Compute the asymptotic population growth rate (lambda) per species

Compute the asymptotic population growth rate (lambda) per species

## Usage

``` r
lambda(mod, pars, stand, env, ctrl = NULL)
```

## Arguments

- mod:

  An `ipm_spModel` object.

- pars:

  An `ipm_parameters` object. Must contain parameters for at least one
  species in `mod`. Lambda is computed only for species present in both
  `mod` and `pars` - other species in `stand` are used as competitors
  but do not appear in the output.

- stand:

  An `ipm_stand` object. Provides size distributions for all species
  (focal and competitors).

- env:

  An `ipm_env` object. Climate drivers.

- ctrl:

  An `ipm_control` object or NULL. Only `bin_width` and `delta_time` are
  used. If NULL, defaults are used.

## Value

An object of S3 class `"ipm_lambda"` - a named numeric vector with one
element per focal species (species with available parameters).

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
lam  <- lambda(mod, pars, s, env, ctrl)
print(lam)
#> <ipm_lambda>  draw=mean
#>   ABIBAL: 2.0762
```
