# Convert Individual Tree Sizes to a Mesh-Aligned Size Distribution

Bins individual tree diameter observations into the mesh point grid
defined by an `init_pop()` output object.

## Usage

``` r
dbh_to_sizeDist(dbh, N_intra)
```

## Arguments

- dbh:

  Numeric vector. Individual tree diameters in mm.

- N_intra:

  List. Size distribution object.

## Value

A copy of `N_intra` with `$Nvec` replaced by the count of individuals in
each mesh size class.
