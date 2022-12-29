# Finding Attractors
Attractors.jl defines a generic interface for finding attractors of dynamical systems. One first decides the instance of [`GeneralizedDynamicalSystem`](@ref) they need. Then, an instance of [`AttractorMapper`](@ref) is created from this dynamical system. This `mapper` instance can be used to compute e.g., [`basins_of_attraction`](@ref), and the output can be further analyzed to get e.g., the [`basin_entropy`](@ref).

```@docs
AttractorMapper
```

## Recurrences
```@docs
AttractorsViaRecurrences
automatic_Δt_basins
```

## Proximity
```@docs
AttractorsViaProximity
```

## Featurizing
```@docs
AttractorsViaFeaturizing
```

## Grouping configurations
### Grouping types
```@docs
GroupingConfig
GroupViaClustering
GroupViaNearestFeature
```

### Grouping utility functions
```@docs
group_features
extract_features
```