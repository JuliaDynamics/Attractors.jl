# Finding Attractors
DynamicalSystems.jl defines a generic interface for finding attractors of dynamical systems. One first decides the instance of [`GeneralizedDynamicalSystem`](@ref) they need. Then, an instance of [`AttractorMapper`](@ref) is created from this dynamical system. This `mapper` instance can be used to compute e.g., [`basins_of_attraction`](@ref), and the output can be further analyzed to get e.g., the [`basin_entropy`](@ref).

The example [2D basins of 4D system](@ref) contains most of the functionality documented in this page.

```@docs
AttractorMapper
AttractorsViaProximity
AttractorsViaRecurrences
automatic_Δt_basins
AttractorsViaFeaturizing
```
