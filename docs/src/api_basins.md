# Basins functions

## Basins types
The basins of attraction are often represented as an array or vector. We also provide a convenient extendable structure
that contains the basins themselves, attractors, and the domains on which the basins are defined. All standard basin-related
functions are compatible with this alternate representation.

```@docs
BasinsOfAttraction
ArrayBasinsOfAttraction
SampledBasinsOfAttraction
```

The `map_to_basin` function provides simple interpolation of a point in state space to determine which basin of attraction it
is likely to belong to:

```@docs
map_to_basin
```

## Calculating basins


Calculating basins of attraction, or their state space fractions, can be done with the functions:
- [`basins_fractions`](@ref)
- [`basins_of_attraction`](@ref)

```@docs
basins_fractions
basins_fractions_labels
basins_of_attraction
statespace_sampler
```

## Sampler interface

Sampling initial conditions is done via the following types:

```@docs
InitialConditionsSampler
RandomICsSampler
PrescribedICs
PerParameterICs
```

## Convergence times

```@docs
convergence_and_basins_fractions
convergence_and_basins_of_attraction
```