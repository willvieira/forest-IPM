# Create an ipm_spModel object defining which species to model

Create an ipm_spModel object defining which species to model

## Usage

``` r
species_model(x, on_missing = "error")
```

## Arguments

- x:

  An `ipm_stand` object or a character vector of species IDs.

- on_missing:

  Character. Behavior when a species ID is not in
  [`supported_species()`](https://willvieira.github.io/forestIPM/reference/supported_species.md).
  One of `"error"` (default), `"drop"`, or `"static"`.

## Value

An object of S3 class `"ipm_spModel"` with fields `$species`, `$params`,
and `$on_missing`.

## Examples

``` r
df <- data.frame(size_mm = c(150, 200, 350),
                 species_id = "ABIBAL",
                 plot_size  = 1000)
s   <- stand(df)
mod <- species_model(s)
print(mod)
#> <ipm_spModel>  1 species: ABIBAL
```
