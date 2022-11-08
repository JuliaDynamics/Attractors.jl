using Neighborhood, Distances
export GroupViaNearestFeature

"""
    GroupViaNearestFeature(templates; kwargs...)

Initialize a struct that contains instructions on how to group features in
[`AttractorsViaFeaturizing`](@ref). `GroupViaNearestFeature` accepts a `template`,
which is a vector of features. Then, generated features from initial conditions in
[`AttractorsViaFeaturizing`](@ref) are labelled according to the feature in
`templates` that is closest (the label is the index of the closest template).

`templates` can be a vector or dictionary mapping keys to templates. Internally
all templates are converted to `SVector` for performance. Hence, it is strongly recommended
that both `templates` and the output of the `featurizer` function in
[`AttractorsViaFeaturizing`](@ref) return `SVector` types.

## Keyword arguments
- `metric = Euclidean()`: metric to be used to quantify distances in the feature space.
- `max_distance = Inf`: Maximum allowed distance between a feature and its nearest
  template for it to be assigned to that template. By default, `Inf` guarantees that a
  feature is assigned to its nearest template regardless of the distance.
  Features that exceed `max_distance` to their nearest template get labelled `-1`.
"""
struct GroupViaNearestFeature{D, T, K <: KDTree, X} <: GroupingConfig
    templates::Vector{SVector{D,T}}
    tree::K # KDTree with templates
    max_distance::T
    dummy_idxs::Vector{Int}
    dummy_dist::Vector{T}
    keys::X
end
function GroupViaNearestFeature(
        templates; metric = Euclidean(), max_distance = Inf
    )
    k, v = collect(keys(templates)), values(templates)
    x = first(v)
    D = length(x); T = eltype(x)
    t = map(x -> SVector{D,T}(x), v)
    # The tree performs the nearest neighbor searches efficiently
    tree = searchstructure(KDTree, t, metric)
    dummy_idxs = [0]; dummy_dist = T[0]
    return GroupViaNearestFeature(t, tree, T(max_distance), dummy_idxs, dummy_dist, k)
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
        return config.keys[label]
    end
end
