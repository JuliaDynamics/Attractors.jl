"""
    GroupingConfig

Supertype for configuration structs on how to group features together.
Used in several occasions such as [`AttractorsViaFeaturizing`](@ref)
or [`aggregate_attractor_fractions`](@ref).

Currently available grouping configurations are:

- [`GroupViaClustering`](@ref)
- [`GroupViaNearestFeature`](@ref)
- [`GroupViaHistogram`](@ref)

`GroupingConfig` defines an extendable interface.
The only thing necessary for a new grouping configuration is to:
1. Make a new type and subtype `GroupingConfig`.
2. Extend the function `group_features(features, config)`.
   If the grouping allows for mapping individual initial conditions to IDs,
   then instead extend the **internal function** `feature_to_group(feature, config)`.
   This will allow doing `id = mapper(u0)` with [`AttractorsViaFeaturizing`](@ref).
   You could still extend `group_features` as well if there are any performance benefits.
3. Include the new grouping file in the `grouping/all_grouping_configs.jl` and list it in
   this documentation string.
"""
abstract type GroupingConfig end
export GroupingConfig

include("cluster_config.jl")
# include("histogram_config.jl")
include("nearest_feature_config.jl")
