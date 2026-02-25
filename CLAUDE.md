# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a R package containing a Bayesian hierarchical Integral Projection Model (IPM) for studying tree population dynamics (population growth rates, coexistence, sensitivity) in eastern North America. Written in R with C++ acceleration via RcppEigen.

## Environment Setup

External data path is stored in `_data.path` (a plain text file pointing to the data directory on the local machine).

## Architecture

### Core Functions (`R/`)

- `vital_rates.R`: Growth (von Bertalanffy), survival (logistic), and ingrowth/recruitment functions parameterized as function of size, competition (basal area), and climate. These are the demographic models fitted with inventory data.
- `kernel.R`: Assembles IPM kernels (P matrix = growth × survival, F matrix = recruitment) using numerical integration over size classes
- `BasalArea_competition.R`: Computes competition metrics (intra/interspecific basal area) from plot data
- `matrix_image.R`: Visualize the IPM kernel discretized matrix
- `src/eigen.cpp`: C++ eigenvalue solver (RcppEigen), compiled via Rcpp on load
