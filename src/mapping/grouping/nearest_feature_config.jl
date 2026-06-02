using Neighborhood, Distances
export GroupViaNearestFeature

"""
    GroupViaNearestFeature(templates; kwargs...)

Initialize a struct that contains instructions on how to group features in
[`AttractorsViaFeaturizing`](@ref). `GroupViaNearestFeature` accepts a `template`,
which is a dictionary mapping unique labels to unique feature vectors.
Then, generated feature vectors from initial conditions in
[`AttractorsViaFeaturizing`](@ref) are labelled according to the feature vector in
`templates` that is closest (the label is the key of the closest template).

`templates` can also be just a vector of feature vectors.

## Keyword arguments
- `metric = Euclidean()`: metric to be used to quantify distances in the feature space.
- `max_distance = Inf`: Maximum allowed distance between a feature and its nearest
  template for it to be assigned to that template. By default, `Inf` guarantees that a
  feature is assigned to its nearest template regardless of the distance.
  Features that exceed `max_distance` to their nearest template get labelled `-1`.
- `use_svector = true`: Convert template vectors to `SVector`.
- `use_kdtree = true`: If `true`, use a KDTree for nearest-neighbor lookup.
  If `false`, use brute-force search over all template vectors.
"""
struct GroupViaNearestFeature{D, T, K, X, M} <: GroupingConfig
    templates::Vector{SVector{D, T}}
    tree::K
    max_distance::T
    dummy_idxs::Vector{Int}
    dummy_dist::Vector{T}
    keys::X
    metric::M
end
function GroupViaNearestFeature(
        templates; metric = Euclidean(), max_distance = Inf,
        use_svector = true, use_kdtree = true,
    )
    k, v = collect(keys(templates)), values(templates)
    x = first(v)
    D = length(x); T = eltype(x)
    t = use_svector ? map(x -> SVector{D, T}(x), v) : v
    tree = use_kdtree ? searchstructure(KDTree, t, metric) : nothing
    dummy_idxs = [0]; dummy_dist = T[0]
    return GroupViaNearestFeature(t, tree, T(max_distance), dummy_idxs, dummy_dist, k, metric)
end

# The following function comes from the source code of the `bulksearch` function
# from Neighborhood.jl. It's the most efficient way to perform one knn search,
# and makes it unnecessary to also implement `group_features`. The bulk version
# has the same performance!
@inbounds function feature_to_group(feature, config::GroupViaNearestFeature)
    (; tree, max_distance, dummy_idxs, dummy_dist, templates, metric) = config
    if isnothing(tree)
        d, label = findmin(i -> metric(feature, templates[i]), eachindex(templates))
    else
        skip = Returns(false)
        sort_result = false
        Neighborhood.NearestNeighbors.knn_point!(
            tree, feature, sort_result, dummy_dist, dummy_idxs, skip
        )
        label = dummy_idxs[1]
        d = dummy_dist[1]
    end
    if d > max_distance
        return -1
    else
        return config.keys[label]
    end
end
