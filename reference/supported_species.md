# Return the tibble of species supported by forestIPM

Return the tibble of species supported by forestIPM

## Usage

``` r
supported_species()
```

## Value

A tibble with six columns: species_id, common_name, nom_commun,
growth_model, surv_model, recruit_model. One row per supported species.

## Examples

``` r
sp <- supported_species()
head(sp)
#> # A tibble: 6 × 11
#>   species_id species_name common_name nom_commun MAT_min MAT_max MAP_min MAP_max
#>   <chr>      <chr>        <chr>       <chr>        <dbl>   <dbl>   <dbl>   <dbl>
#> 1 ABIBAL     Abies balsa… Balsam fir  Sapin bau…    -5.2     9      433     1859
#> 2 ACERUB     Acer rubrum  Red maple   Erable ro…    -1.5    24.7    436     2688
#> 3 ACESAC     Acer saccha… Sugar maple Erable a …    -0.5    16.7    443     1944
#> 4 BETALL     Betula alle… Yellow bir… Bouleau j…    -1.4    15.1    529     2177
#> 5 BETPAP     Betula papy… White birch Bouleau b…    -5      13.5    395     1859
#> 6 CARGLA     Carya glabra Pignut hic… Caryer gl…     6.3    23.3    631.    2326
#> # ℹ 3 more variables: growth_model <chr>, surv_model <chr>, recruit_model <chr>
```
