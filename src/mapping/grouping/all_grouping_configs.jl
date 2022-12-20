"""
    GroupingConfig

Supertype for configuration structs on how to group features together.
Used in several occasions such as [`AttractorsViaFeaturizing`](@ref).

`GroupingConfig` defines an extendable interface:
The only thing necessary for a new grouping configuration is to:
1. Make a new type and subtype `GroupingConfig`.
2. Extend the function `group_features(features, config)`.
   If the grouping allows for mapping individual initial conditions to IDs,
   then instead extend the **internal function** `feature_to_group(feature, config)`.
   You can still extend `group_features` as well if there are any performance benefits.
3. Include the new grouping file in the `grouping/all_grouping_configs.jl`.
"""
abstract type GroupingConfig end
export GroupingConfig

include("cluster_config.jl")
# include("histogram_config.jl")
include("nearest_feature_config.jl")
