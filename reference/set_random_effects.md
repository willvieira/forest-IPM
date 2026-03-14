# Set plot-level random effects for one or more species

Populates the `$random_effects` slot inside each species entry of an
`ipm_parameters` object. Engines
([`lambda()`](https://willvieira.github.io/forestIPM/reference/lambda.md),
[`project()`](https://willvieira.github.io/forestIPM/reference/project.md))
default to `c(0, 0, 0)` for any species whose slot is `NULL`.

## Usage

``` r
set_random_effects(pars, values, species = NULL)
```

## Arguments

- pars:

  An `ipm_parameters` object.

- values:

  A numeric vector of length 3 (growth, survival, recruitment offsets),
  applied to all species named in `species`; *or* a named list keyed by
  species ID, each element a numeric vector of length 3.

- species:

  Character vector of species IDs. Used only when `values` is a numeric
  vector; ignored when `values` is a named list. If `NULL` and `values`
  is a numeric vector, applies to all species.

## Value

The modified `ipm_parameters` object.

## Examples

``` r
df <- data.frame(size_mm = c(150, 200, 350),
                 species_id = "ABIBAL",
                 plot_size  = 1000)
s     <- stand(df)
mod   <- species_model(s)
pars  <- parameters(mod, draw = "mean")
pars2 <- set_random_effects(pars, values = c(0, 0, 0))
```
