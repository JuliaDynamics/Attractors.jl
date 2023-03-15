# Attractors.jl

```@docs
Attractors
```

```@setup MAIN
using CairoMakie, Attractors
include("../basins_plotting.jl")
```

## Outline of Attractors.jl

![TransitionIndicators.jl](https://github.com/JuliaDynamics/JuliaDynamics/blob/master/videos/attractors/attractorsjl_overview.png?raw=true)


1. First be sure that you are aware of what is a [`DynamicalSystem`](@ref). This is the input to the whole infrastructure of Attractors.jl.
2. The bulk of the work in Attractors.jl is done by the [`AttractorMapper`](@ref) type, that instructs how to find attractors and maps initial conditions to them. It can be used in functions like [`basins_fractions`](@ref).
3. For grouping features, there is a sub-infrastructure for instructing how to group features, which is governed by [`GroupingConfig`](@ref).
4. The infrastructure of finding attractors and their basins fractions is then integrated into a brand new way of doing bifurcation analysis in the [`continuation`](@ref) function.
5. See [Examples for Attractors.jl](@ref) for several applications in real world cases.

## `DynamicalSystem` reference

The kinds of dynamical systems that can be used in Attractors.jl are listed below for reference
```@docs
DynamicalSystem
DeterministicIteratedMap
CoupledODEs
StroboscopicMap
PoincareMap
ProjectedDynamicalSystem
ArbitrarySteppable
```

## `StateSpaceSet` reference

```@docs
StateSpaceSet
```