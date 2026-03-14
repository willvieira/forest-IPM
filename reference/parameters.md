# Resolve a single parameter realization from Bayesian posteriors

Resolve a single parameter realization from Bayesian posteriors

## Usage

``` r
parameters(mod, draw = "random", seed = NULL)
```

## Arguments

- mod:

  An `ipm_spModel` object.

- draw:

  Character `"mean"` or `"random"`, or a positive integer 1-1000
  selecting a specific posterior draw index.

- seed:

  Integer or NULL. Random seed. When `draw = "random"` and
  `seed = NULL`, a seed is auto-generated so the draw is reproducible.
  Retrieve it from `$seed` on the returned object.

## Value

An object of S3 class `"ipm_parameters"`.

## Examples

``` r
df <- data.frame(size_mm = c(150, 200, 350),
                 species_id = "ABIBAL",
                 plot_size  = 1000)
s    <- stand(df)
mod  <- species_model(s)
pars <- parameters(mod, draw = "mean")
print(pars)
#> <ipm_parameters>  draw=mean | 1 species
```
