# Basin maps

## Finding basins and attractors

Attractors.jl defines a generic and extendable interface for finding mapping initial conditions to their basins of attraction and finding their corresponding attractors. One first decides the instance of [`DynamicalSystem`](@ref) they need. Then, an instance of [`BasinMap`](@ref) is created from this dynamical system. This `bmap` instance can be used to compute e.g., [`basins_of_attraction`](@ref), and the output can be further analyzed to get e.g., the [`basin_entropy`](@ref).

```@docs
BasinMap
extract_attractors
convergence_time
```


## Recurrences-based

```@docs
BasinMapRecurrences
automatic_Δt_recurrences
SubdivisionBasedGrid
subdivision_based_grid
```

## Proximity-based

```@docs
BasinMapProximity
```

## Featurize-and-group

```@docs
BasinMapFeaturizeGroup
```


## Grouping configurations

Grouping configurations that can be given to [`BasinMapFeaturizeGroup`](@ref)
are part of a generic and extendable interface based on the [`group_features`](@ref)
function.
The grouping configuration sets how the features describing the trajectories will be grouped together.
Nevertheless, this grouping infrastructure can also be used and extended completely independently of finding attractors of dynamical systems!

### Grouping interface

```@docs
group_features
GroupingConfig
```

### Grouping types

```@docs
GroupViaClustering
GroupViaHistogram
GroupViaNearestFeature
GroupViaPairwiseComparison
```

### Grouping utils

```@docs
extract_features
```
