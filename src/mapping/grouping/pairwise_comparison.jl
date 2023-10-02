export GroupViaPairwiseComparison


"""
    GroupViaPairwiseComparison(; distance_threshold::Real, kwargs...)

Initialize a struct that contains instructions on how to group features in
[`AttractorsViaFeaturizing`](@ref). `GroupViaPairwiseComparison` groups features and
identifies clusters by considering the pairwise distance (using `distance_metric`) between
features. 

## Keyword arguments
    
* `distance_threshold`: a real number defining the maximum distance two features can be to
  be considered in the same cluster - above the threshold, features are different. This
  value simply needs to be large enough to differentiate clusters.
* `distance_metric = Euclidean()`: A metric to be used in the clustering. It can be any
  function `f(a, b)` that returns the distance between any type of data structure (usually
  vectors or matrices of reals). Needs to be consistent with the `featurizer` function. All
  metrics from Distances.jl can be used here.
* `rescale_features = true`: if true, rescale each dimension of the extracted features
  separately into the range `[0,1]`. This typically leads to more accurate clustering.

## Description
This algorithm assumes that the the features well-separated into distinct clouds, with the
maximum radius of the cloud controlled by `distance_threshold`. Since the systems are
deterministic, this is achievable with a good-enough `featurizer` function, by removing
transients, and running the trajectories for sufficiently long. It then considers that
features belong to the same attractor when their pairwise distance, computed using
`distance_metric`, is smaller or equal than `distance_threshold`, and belong to different
attractors when the distance is bigger than the threshold. Attractors correspond to each
grouping of similar features. In this way, the key parameter `distance_threshold` is
simply the amount of variations permissible in the features. If they are well-chosen, the
value can be relatively small and does not need to be fine tuned.
"""
struct GroupViaPairwiseComparison{R<:Real, M} <: GroupingConfig
    distance_threshold::R
    distance_metric::M
    rescale_features::Bool 
end

function GroupViaPairwiseComparison(;
        distance_threshold, #impossible to set a good default value, depends on the features
        distance_metric=Euclidean(), rescale_features=false, 
    )
    return GroupViaPairwiseComparison(
        distance_threshold,
        distance_metric, rescale_features,
    )
end

function group_features(
    features, config::GroupViaPairwiseComparison;
    kwargs...
)
    if config.rescale_features
        features = _rescale_to_01(features)
    end
    
    labels = _cluster_features_into_labels(features, config, config.distance_threshold; kwargs...)
    return labels
end

# TODO: add support for par_weight,plength and spp in the computation of the distance metric?
function _cluster_features_into_labels(features, config::GroupViaPairwiseComparison, distance_threshold::Real; kwargs...)
    labels_features = Vector{Int64}(undef, length(features)) #labels of all features
    metric = config.distance_metric
    
    # Assign feature 1 as a new attractor
    labels_features[1] = 1
    cluster_idxs = [1] # idxs of the features that define each cluster
    cluster_labels = [1] # labels for the clusters, going from 1 : num_clusters
    next_cluster_label = 2
    
    for idx_feature = 2:length(features)
        feature = features[idx_feature]
        dist_to_clusters = _distance_dict(feature, features, cluster_idxs, cluster_labels, metric; kwargs...)
        min_dist, closest_cluster_label = findmin(dist_to_clusters)
        
        if min_dist > distance_threshold #bigger than threshold => new attractor
            feature_label = next_cluster_label
            push!(cluster_idxs, idx_feature)
            push!(cluster_labels, next_cluster_label)
            # @info "New attractor $next_cluster_label, min dist was $min_dist > $distance_threshold" #TODO: allow this when debugging verbose mode on!
            next_cluster_label += 1
        else #smaller than threshold => assign to closest cluster 
            feature_label = closest_cluster_label
        end
        
        labels_features[idx_feature] = feature_label 
    end 
    return labels_features
end

function _distance_dict(feature, features, cluster_idxs, cluster_labels, metric; kwargs...)
    if metric isa Metric
        dist_to_clusters = Dict(cluster_label => evaluate(metric, feature, features[cluster_idxs[idx_cluster]]) for (idx_cluster, cluster_label) in enumerate(cluster_labels))
    else
        dist_to_clusters = Dict(cluster_label => metric(feature, features[cluster_idxs[idx_cluster]]) for (idx_cluster, cluster_label) in enumerate(cluster_labels))
    end
    
    return dist_to_clusters 
end