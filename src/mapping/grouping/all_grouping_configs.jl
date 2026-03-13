export GroupingConfig, group_features

"""
   GroupingConfig

Supertype for configuration structs on how to group features together.
Used in several occasions such as [`AttractorsViaFeaturizing`](@ref)
or [`aggregate_attractor_fractions`](@ref).

Currently available grouping configurations are:

- [`GroupViaClustering`](@ref)
- [`GroupViaNearestFeature`](@ref)
- [`GroupViaHistogram`](@ref)
- [`GroupViaPairwiseComparison`](@ref)

## For developers

`GroupingConfig` defines an extendable interface.
The only thing necessary for a new grouping configuration is to:

1. Make a new type and subtype `GroupingConfig`.
2. If the grouping allows for mapping individual features to group index,
   then instead extend the **internal function** `feature_to_group(feature, config)`.
   This will also allow doing `id = mapper(u0)` with [`AttractorsViaFeaturizing`](@ref).
3. Else, extend the function `group_features(features, config)`.
   You could still extend `group_features` even if (2.) is satisfied,
   if there are any performance benefits.
4. Include the new grouping file in the `grouping/all_grouping_configs.jl` and list it in
   this documentation string.
"""
abstract type GroupingConfig end


"""
    group_features(features, group_config::GroupingConfig) → labels

Group the given iterable of "features" (anything that can be grouped, typically vectors of real numbers)
according to the configuration and return the labels (vector of equal length as `features`).
See [`GroupingConfig`](@ref) for possible grouping configuration configurations.
"""
function group_features(features, group_config::GroupingConfig)
    return map(f -> feature_to_group(f, group_config), features)
end

"""
    feature_to_group(feature::AbstractVector, group_config::GroupingConfig) → group_label

Map the given feature vector to its group label (integer).
This is an internal function. It is strongly recommended that `feature isa SVector`.
"""
function feature_to_group(feature, group_config::GroupingConfig)
    throw(
        ArgumentError(
            """
            `feature_to_group` not implemented for config $(nameof(typeof(group_config))).
            """
        )
    )
end

include("cluster_config.jl")
include("histogram_config.jl")
include("nearest_feature_config.jl")
include("pairwise_comparison.jl")
