forestIPM — API Architecture Specification
Overview

This package implements a dynamic, species-specific forest Integral Projection Model (IPM) framework capable of:

Computing population growth rate (lambda)

Projecting single-species dynamics

Projecting multi-species community dynamics with competition feedback

Incorporating climate drivers (MAT, MAP)

Incorporating posterior parameter uncertainty (4000 draws per species)

Including plot-level random effects (growth, mortality, recruitment)

Supporting future spatial prediction of random effects

The API is designed to:

Keep a minimal number of object types

Treat single-species models as a special case of multi-species models

Strictly separate structure, ecological state, parameter realization, environmental drivers, and simulation control

Remain extensible to spatial random-effect modeling

Core Object System

The package revolves around five core objects:

stand

species_model

parameters

env_condition

control

Followed by two main engines:

lambda()

project()

1. stand()
Purpose

Defines the ecological state of a plot at time t₀.

This is the foundational object of the workflow.

Responsibilities

Stores all individuals in the plot

Stores species identity per individual

Stores individual sizes

Stores plot size

Defines the initial competition state

Defines which species are present

Example
stand <- stand(
  data = df,              # must contain size + species columns
  plot_size = 400         # m²
)
Internal Structure
class: "ipm_stand"

stand$trees        # individual-level data
stand$species      # unique species vector
stand$plot_size
2. species_model()
Purpose

Defines the structural IPM(s) for one or multiple species.

Contains:

Vital rate functional forms

Kernel constructors

Required parameter names

Scaling rules for MAT/MAP

Random-effect slots (structure only, no values)

Important Design Decision

Single- and multi-species models share the same class.

Usage
Explicit species
mod <- species_model(c("Picea_sp", "Abies_sp"))
Infer from stand
mod <- species_model(stand)

Internally:

mod$species  # character vector
class(mod)   # "ipm_model"
3. parameters()
Purpose

Defines a realization of model parameters for all species in mod.

Handles:

Posterior draw selection

Reproducibility via seed

Random effects (growth, mortality, recruitment)

Key Design Rule

Parameters are always stored as a multi-species structure internally — even if only one species.

No separate community_params class exists.

Usage
pars <- parameters(
  mod,
  draw = "random",   # or numeric index
  seed = 42
)
Internal Structure
class: "ipm_parameters"

pars$species_params
  $Picea_sp
      $fixed
      $random_effects
      $draw_id
  $Abies_sp
      ...
Random Effects

Random effects live inside parameters().

They represent a model realization, not ecological state.

Future extension:

parameters(
  mod,
  draw = 12,
  random_effect = spatial_re(lat, lon)
)
4. env_condition()
Purpose

Defines environmental drivers (climate only).

Contains no ecological state.

Responsibilities

MAT (numeric or function of time)

MAP (numeric or function of time)

Important Design Decision

delta_time does NOT live here.

Example
env <- env_condition(
  MAT = 8,
  MAP = function(t) 1200 + 5 * t
)

Internally:

Accept numeric or function

Always scale MAT/MAP internally to [0,1]

5. control()
Purpose

Defines simulation meta-parameters.

Responsibilities

years (default = 1)

delta_time (default = 1)

Optional verbosity / output controls

Example
ctrl <- control(
  years = 100,
  delta_time = 1
)

If not provided:

years = 1
delta_time = 1
Core Engines
6. lambda()
Purpose

Compute population growth rate(s).

Signature
lambda(mod, pars, stand, env)
Behavior

Returns a named numeric vector

One lambda per species in mod

Example:

Picea_sp  = 1.03
Abies_sp  = 0.97
Betula_sp = 1.11
Important Rule

lambda() does NOT require parameters for heterospecific species.

Heterospecific competition is computed from stand structure only.

No dynamic coupling occurs.

7. project()
Purpose

Project dynamics over time.

Signature
project(mod, pars, stand, env, ctrl)
Behavior Depends on Species Count
Single Species

Independent dynamic

Stand updated through time

Multiple Species

Fully coupled dynamics

At each timestep:

Compute competition from updated stand

Build each species kernel

Project each species

Update stand

Repeat

Parameter Requirements

pars must contain parameters for all species in mod

Species in stand must be included in mod

Validation is automatic.

Full Workflow

Canonical usage:

stand <- stand(data, plot_size = 400)

mod   <- species_model(stand)

pars  <- parameters(mod, draw = "random", seed = 42)

env   <- env_condition(
  MAT = 8,
  MAP = 1200
)

ctrl  <- control(
  years = 100,
  delta_time = 1
)

lambda(mod, pars, stand, env)

project(mod, pars, stand, env, ctrl)
Object Responsibility Separation
Object	Responsibility
stand	Ecological state
species_model	Structural IPM definitions
parameters	Parameter realization + random effects
env_condition	Climate drivers
control	Simulation meta-settings

No object overlaps responsibility.

Internal Execution Flow (Projection)

For each timestep:

Retrieve MAT/MAP from env

Scale climate internally

Compute competition from stand

Retrieve species-specific parameters

Construct kernel per species

Apply projection

Update stand state

Repeat until ctrl$years

Reproducibility Rules

All stochasticity controlled through parameters(draw, seed)

No hidden randomness inside projection

Climate functions are deterministic unless user-defined otherwise

Extensibility

This architecture supports:

Spatial prediction of random effects

Adding trait-based modifiers

Disturbance modules

Alternative competition formulations

Parallel posterior projections

Future GPU acceleration

Without altering the API surface.

Final Design Summary

This API:

Uses a minimal number of object classes

Treats single-species as a special case of multi-species

Separates structure, state, drivers, and parameters cleanly

Makes stand the ecological anchor

Keeps random effects in parameters

Places delta_time only in control

Makes lambda() lightweight

Makes project() dynamically coupled when needed

This is a modeling framework architecture, not just a wrapper.
