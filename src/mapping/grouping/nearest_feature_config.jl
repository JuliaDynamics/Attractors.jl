using Neighborhood, Distances
export GroupViaNearestFeature

"""
    GroupViaNearestFeature(templates; kwargs...)

Initialize a struct that contains instructions on how to group features in
[`AttractorsViaFeaturizing`](@ref). `GroupViaNearestFeature` accepts a `template`,
which is a vector of features. Then, generated features from initial conditions in
[`AttractorsViaFeaturizing`](@ref) are labelled according to the feature in
`templates` that is closest (the label is the index of the closest template).

For optimization, `templates` is expected to be a `Vector{<:SVector}`, leading
to hyper-optimized nearest neighbor searches. For the same reason, it is recommended
that the `featurizer` function given in [`AttractorsViaFeaturizing`](@ref) returns
`SVector`.

## Keyword arguments
- `metric = Euclidean()`: metric to be used to quantify distances in the feature space.
- `max_distance = Inf`: Maximum allowed distance between a feature and its nearest
  template for it to be assigned to that template. By default, `Inf` guarantees that a
  feature is assigned to its nearest template regardless of the distance.
  Features that exceed `max_distance` to their nearest template get labelled `-1`.
"""
struct GroupViaNearestFeature{D, T, K} <: GroupingConfig
    templates::Vector{SVector{D,T}}
    max_distance::T
    template_tree::K
    dummy_idxs::Vector{Int}
    dummy_dist::Vector{T}
end
function GroupViaNearestFeature(
        templates::Vector{<:AbstractVector}; metric = Euclidean(), max_distance = Inf
    )
    x = first(templates)
    D = length(x); T = eltype(x)
    t = if templates isa Vector{<:SVector}
        templates
    else
        map(x -> SVector{D,T}(x), templates)
    end
    # The tree performs the nearest neighbor searches efficiently
    tree = searchstructure(KDTree, templates, metric)
    dummy_idxs = [0]; dummy_dist = T[0]
    return GroupViaNearestFeature(templates, tree, T(max_distance), dummy_idxs, dummy_dist)
end

# The following function comes from the source code of the `bulksearch` function
# from Neighborhood.jl. It's the most efficient way to perform one knn search,
# and makes it unecessary to also implement `group_features`. The bulk version
# has the same performance!
@inbounds function feature_to_group(feature, config::GroupViaNearestFeature)
    (; tree, max_distance, dummy_idxs, dummy_dist) = config
    skip = Neighborhood.NearestNeighbors.always_false
    sort_result = false
    Neighborhood.NearestNeighbors.knn_point!(
        tree, feature, sort_result, dummy_dist, dummy_idxs, skip
    )
    label = dummy_idxs[1]
    d = dummy_dist[1]
    if d > max_distance
        return -1
    else
        return label
    end
end
