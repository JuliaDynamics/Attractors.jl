export GroupViaPairwiseComparison

"""
    GroupViaPairwiseComparison(; threshold::Real, kwargs...)

Initialize a struct that contains instructions on how to group features in
[`AttractorsViaFeaturizing`](@ref). `GroupViaPairwiseComparison` groups features and
identifies clusters by considering the pairwise distance between features. It can be used
as an alternative to the clustering method in `GroupViaClustering`, having the
advantage that it is simpler, typically faster and uses less memory. 

## Keyword arguments
    
* `threshold`: a real number defining the maximum distance two features can be to
  be considered in the same cluster - above the threshold, features are different. This
  value simply needs to be large enough to differentiate clusters.
* `metric = Euclidean()`: A metric to be used in the clustering. It can be any
  function `f(a, b)` that returns the distance between any type of data structure (usually
  vectors or matrices of reals). Needs to be consistent with the `featurizer` function. All
  metrics from Distances.jl can be used here.
* `rescale_features = true`: if true, rescale each dimension of the extracted features
  separately into the range `[0,1]`. This typically leads to more accurate clustering.

## Description
This algorithm assumes that the features are well-separated into distinct clouds, with the
maximum radius of the cloud controlled by `threshold`. Since the systems are
deterministic, this is achievable with a good-enough `featurizer` function, by removing
transients, and running the trajectories for sufficiently long. It then considers that
features belong to the same attractor when their pairwise distance, computed using
`metric`, is smaller than or equal to `threshold`, and that they belong
to different attractors when the distance is bigger. Attractors correspond to each
grouping of similar features. In this way, the key parameter `threshold` is
basically the amount of variation permissible in the features belonging to the same
attractor. If they are well-chosen, the value can be relatively small and does not need to
be fine tuned. 

The `threshold` should achieve a balance: one one hand, it should be large enough
to account for variations in the features from the same attractor - if it's not large
enough, the algorithm will find duplicate attractors. On the other hand, it should be
small enough to not group together features from distinct attractors. This requires some
knowledge of how spread the features are. If it's too big, the algorithm will miss some
attractors, as it groups 2+ distinct attractors together. Therefore, as a rule of thumb,
one can repeat the procedure a few times, starting with a relatively large value and
reducing it until no more attractors are found and no duplicates appear.

The method uses relatively little memory, as it only stores vectors whose size is on order
of the number of attractors of the system. 
"""
struct GroupViaPairwiseComparison{R<:Real, M} <: GroupingConfig
    threshold::R
    metric::M
    rescale_features::Bool 
end

function GroupViaPairwiseComparison(;
        threshold, #impossible to set a good default value, depends on the features
        metric=Euclidean(), rescale_features=false, 
    )
    return GroupViaPairwiseComparison(
        threshold,
        metric, rescale_features,
    )
end

function group_features(
    features, config::GroupViaPairwiseComparison;
    kwargs...
)
    if config.rescale_features
        features = _rescale_to_01(features)
    end
    
    labels = _cluster_features_into_labels(features, config, config.threshold; kwargs...)
    return labels
end

# TODO: add support for par_weight,plength and spp in the computation of the distance metric?
function _cluster_features_into_labels(features, config::GroupViaPairwiseComparison, threshold::Real; kwargs...)
    labels_features = Vector{Int64}(undef, length(features)) #labels of all features
    metric = config.metric
    
    # Assign feature 1 as a new attractor
    labels_features[1] = 1
    cluster_idxs = [1] # idxs of the features that define each cluster
    cluster_labels = [1] # labels for the clusters, going from 1 : num_clusters
    next_cluster_label = 2
    
    for idx_feature = 2:length(features)
        feature = features[idx_feature]
        dist_to_clusters = _distance_dict(feature, features, cluster_idxs, cluster_labels, metric; kwargs...)
        min_dist, closest_cluster_label = findmin(dist_to_clusters)
        
        if min_dist > threshold #bigger than threshold => new attractor
            feature_label = next_cluster_label
            push!(cluster_idxs, idx_feature)
            push!(cluster_labels, next_cluster_label)
            # @info "New attractor $next_cluster_label, min dist was $min_dist > $threshold" #TODO: allow this when debugging verbose mode on!
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