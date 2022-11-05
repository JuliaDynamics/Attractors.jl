using Neighborhood, Distances

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
  Features that exceed `max_distance` to all `templates` get labelled `-1` (divergent).
"""
struct GroupViaNearestFeature{D, T, M <: Metric}
    templates::Vector{SVector{D,T}}
    metric::M
    max_distance::T
end
function GroupViaNearestFeature(templates; metric = Euclidean(), max_distance = Inf)
    t = if templates isa Vector{<:AbstractVector}
        templates
    else
        x = first(templates)
        D = length(x); T = eltype(x)
        map(x -> SVector{D,T}(x), templates)
    end
    return GroupViaNearestFeature(templates, metric, max_distance)
end

function group_features(features::Vector{<:AbstractVector}, config::GroupViaNearestFeature)
   # prepare for nearest-neighbors algorithm (kNN with k=1)
   template_tree = searchstructure(KDTree, config.templates, config.metric)
   labels_nested, distances_nested = bulksearch(template_tree, features, NeighborNumber(1))
   labels = reduce(vcat, labels_nested) # make it a vector
   distances = reduce(vcat, distances_nested)
   # Make label -1 if error strictly bigger than threshold
   labels[distances .> config.max_distance_template] .= -1
   # labels[i] is the label of template nearest to feature i
   # (i-th column of features matrix)
   return labels
end
