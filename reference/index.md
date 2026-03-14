# Package index

## IPM Constructors

Functions to build IPM components

- [`stand()`](https://willvieira.github.io/forestIPM/reference/stand.md)
  : Create an ipm_stand object representing a forest plot
- [`species_model()`](https://willvieira.github.io/forestIPM/reference/species_model.md)
  : Create an ipm_spModel object defining which species to model
- [`env_condition()`](https://willvieira.github.io/forestIPM/reference/env_condition.md)
  : Create an ipm_env object specifying climate drivers
- [`parameters()`](https://willvieira.github.io/forestIPM/reference/parameters.md)
  : Resolve a single parameter realization from Bayesian posteriors
- [`control()`](https://willvieira.github.io/forestIPM/reference/control.md)
  : Configure IPM projection settings

## IPM Engines

Functions to run IPM calculations

- [`lambda()`](https://willvieira.github.io/forestIPM/reference/lambda.md)
  : Compute the asymptotic population growth rate (lambda) per species
- [`project()`](https://willvieira.github.io/forestIPM/reference/project.md)
  : Project population dynamics over time

## Utilities

Helper functions

- [`supported_species()`](https://willvieira.github.io/forestIPM/reference/supported_species.md)
  : Return the tibble of species supported by forestIPM
- [`scale_env()`](https://willvieira.github.io/forestIPM/reference/scale_env.md)
  : Scale environmental values to \[0, 1\]
- [`unscale_env()`](https://willvieira.github.io/forestIPM/reference/unscale_env.md)
  : Unscale environmental values from \[0, 1\] to natural units
- [`set_random_effects()`](https://willvieira.github.io/forestIPM/reference/set_random_effects.md)
  : Set plot-level random effects for one or more species

## Visualization

Plotting and diagnostic functions

- [`plot(`*`<ipm_projection>`*`)`](https://willvieira.github.io/forestIPM/reference/plot.ipm_projection.md)
  : Plot an ipm_projection object
