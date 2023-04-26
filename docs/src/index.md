# Attractors.jl

```@docs
Attractors
```

```@setup MAIN
using CairoMakie, Attractors
```

## Latest news

Our pre-print that discusses the global stability analysis framework offered by Attractors.jl ([`continuation`](@ref)) and the novel continuation offered by [`RecurrencesFindAndMatch`](@ref) is now online: https://arxiv.org/abs/2304.12786 !


## Outline of Attractors.jl

![Attractors.jl flowchart](https://github.com/JuliaDynamics/JuliaDynamics/blob/master/videos/attractors/attractorsjl_overview.png?raw=true)


1. First be sure that you are aware of what is a [`DynamicalSystem`](@ref). This is the input to the whole infrastructure of Attractors.jl.
2. The bulk of the work in Attractors.jl is done by the [`AttractorMapper`](@ref) type, that instructs how to find attractors and maps initial conditions to them. It can be used in functions like [`basins_fractions`](@ref).
3. For grouping features, there is a sub-infrastructure for instructing how to group features, which is governed by [`GroupingConfig`](@ref).
4. The infrastructure of finding attractors and their basins fractions is then integrated into a brand new way of doing bifurcation analysis in the [`continuation`](@ref) function.
5. See [Examples for Attractors.jl](@ref) for several applications in real world cases.
