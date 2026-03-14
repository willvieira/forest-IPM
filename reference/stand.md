# Create an ipm_stand object representing a forest plot

Create an ipm_stand object representing a forest plot

## Usage

``` r
stand(data)
```

## Arguments

- data:

  A data.frame or tibble with individual-level tree records. Required
  columns: a size column (named `size` or `dbh`, in mm,= 127 mm), a
  species column (named `sp`, `species`, or`species_id`), and a
  `plot_size` column (plot area in m2,same value for all rows).

## Value

An object of S3 class `"ipm_stand"` with fields `$trees` (data.frame
with standardized columns `size_mm` and `species_id`), `$species`
(character vector of unique species IDs), and `$plot_size` (numeric
scalar).

## Examples

``` r
df <- data.frame(size_mm = c(150, 200, 350),
                 species_id = "ABIBAL",
                 plot_size  = 1000)
s <- stand(df)
print(s)
#> <ipm_stand>  3 trees | 1 species | 1000 m2 plot
```
