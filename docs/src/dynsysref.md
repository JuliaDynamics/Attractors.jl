
# `DynamicalSystem` reference

This page is a convenience reference to some of the contents of the DynamicalSystemsBase.jl package (one of the core modules of DynamicalSystems.jl).

## Dynamical systems

The kinds of dynamical systems that can be used in Attractors.jl are listed below for reference

- [`DynamicalSystem`](@ref)
- [`DeterministicIteratedMap`](@ref)
- [`CoupledODEs`](@ref)
- [`StroboscopicMap`](@ref)
- [`PoincareMap`](@ref)
- [`ProjectedDynamicalSystem`](@ref)
- [`ArbitrarySteppable`](@ref)

```@docs
DynamicalSystem
DeterministicIteratedMap
CoupledODEs
StroboscopicMap
PoincareMap
ProjectedDynamicalSystem
ArbitrarySteppable
```

## Relevant dynamical systems API

Here we state only the relevant API functions; others can be found in the DynamicalSystemsBase.jl documentation.

```@docs
step!(::DynamicalSystem, args...; kwargs...)
reinit!(::DynamicalSystem, args...; kwargs...)
set_parameter!
```

## `StateSpaceSet` reference

```@docs
StateSpaceSet
```
