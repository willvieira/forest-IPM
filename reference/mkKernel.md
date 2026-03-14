# Build the IPM Kernel

Assembles the full Integral Projection Model kernel (K = P + F) for a
focal species given its size distribution vectors and environmental
conditions.

## Usage

``` r
mkKernel(
  Nvec_intra,
  Nvec_inter,
  delta_time,
  plotSize,
  Temp,
  Prec,
  pars,
  plot_random
)
```

## Arguments

- Nvec_intra:

  List. Intraspecific size distribution.

- Nvec_inter:

  List. Interspecific size distribution.

- delta_time:

  Numeric. Time step in years.

- plotSize:

  Numeric. Plot area in square meters.

- Temp:

  Numeric. Mean annual temperature (degrees Celsius).

- Prec:

  Numeric. Mean annual precipitation (mm).

- pars:

  Named list. Species-specific parameters with elements `growth`,
  `mort`, `rec`, and `sizeIngrowth`.

- plot_random:

  Numeric vector of length 3. Plot-level random effects for growth (1),
  mortality (2), and recruitment (3).

## Value

A named list with elements `K` (full kernel), `P` (growth x survival
kernel), and `F` (recruitment kernel), each a square matrix of dimension
equal to the number of mesh points.
