# Finding Attractors
Attractors.jl defines a generic interface for finding attractors of dynamical systems. One first decides the instance of [`GeneralizedDynamicalSystem`](@ref) they need. Then, an instance of [`AttractorMapper`](@ref) is created from this dynamical system. This `mapper` instance can be used to compute e.g., [`basins_of_attraction`](@ref), and the output can be further analyzed to get e.g., the [`basin_entropy`](@ref).

The example [2D basins of 4D system](@ref) contains most of the functionality documented in this page.

```@docs
AttractorMapper
```

## Proximity
```@docs
AttractorsViaProximity
```

## Recurrences
```@docs
AttractorsViaRecurrences
AttractorsViaRecurrencesSparse
automatic_Δt_basins
```
## Featurizing
```@docs
AttractorsViaFeaturizing
ClusteringConfig
cluster_features
```
