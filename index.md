# forestIPM

[![R-CMD-check](https://github.com/willvieira/forestIPM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/willvieira/forestIPM/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/github/willvieira/forestIPM/graph/badge.svg)](https://app.codecov.io/github/willvieira/forestIPM)

A Bayesian hierarchical Integral Projection Model (IPM) for studying
tree population dynamics in eastern North America. The package
implements growth, survival, and recruitment vital rate models
parameterized from forest inventory data, enabling researchers to
compute species specific population growth rates ($\lambda$) and project
plot-level stand dynamics.

## Documentation

Full documentation, methodology, and worked examples are available in
the companion book:

**[Forest Demography IPM
Book](https://willvieira.github.io/book_forest-demography-IPM/)**

The [function
reference](https://willvieira.github.io/forestIPM/reference/index.html)
is available on the package website.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("willvieira/forestIPM")
```

## Citation

If you use this package in your research, please cite the related
[article](https://willvieira.github.io/ms_forest-ipm-sensitivity/).
