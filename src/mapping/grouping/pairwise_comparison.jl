export GroupViaPairwiseComparison

"""
    GroupViaPairwiseComparison <: GroupingConfig
    GroupViaPairwiseComparison(; threshold::Real, metric...)

Initialize a struct that contains instructions on how to group features in
[`AttractorsViaFeaturizing`](@ref). `GroupViaPairwiseComparison` groups features and
identifies clusters by considering the pairwise distance between features. It can be used
as an alternative to the clustering method in `GroupViaClustering`, having the
advantage that it is simpler, typically faster and uses less memory.

## Keyword arguments

* `threshold = 0.1`: A real number defining the maximum distance two features can have to
  be considered in the same cluster - above the threshold, features are different. This
  value simply needs to be large enough to differentiate clusters.
  A good value for `threshold` depends on the feature variability within a cluster
  the chosen metric, and whether features are rescaled. See description below for more.
* `metric = Euclidean()`: A function `metric(a, b)` that returns the distance between two
  features `a` and `b`, outputs of `featurizer`. Any `Metric` from Distances.jl can be used
  here.
* `rescale_features = true`: if true, rescale each dimension of the extracted features
  separately into the range `[0, 1]`. This typically leads to more accurate grouping
  for the default `metric. threshold`, however, it should be avoid for when
  there is only one attractor for the system, because it leads to it being
  wrongly classified as many different attractors at the same location).

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

The method scales as O(N) in memory and performance with N the number of features.
This is a huge difference versus the O(N^2) of [`GroupViaClustering`](@ref).
"""
@kwdef struct GroupViaPairwiseComparison{R<:Real, M} <: GroupingConfig
    threshold::R = 0.1
    metric::M = Euclidean()
    rescale_features::Bool = true
end

function group_features(features, config::GroupViaPairwiseComparison)
    if config.rescale_features
        features = _rescale_to_01(features)
    end
    labels = _cluster_features_into_labels(features, config)
    return labels
end

function _cluster_features_into_labels(features, config::GroupViaPairwiseComparison)
    labels = ones(Int, length(features))
    metric = config.metric
    threshold = config.threshold

    # Assign feature 1 as a new cluster center
    labels[1] = 1
    cluster_idxs = [1] # idxs of the features that define each cluster
    cluster_labels = [1] # labels for the clusters, going from 1 : num_clusters
    next_cluster_label = 2

    for idx_feature = 2:length(features)
        feature = features[idx_feature]
        # because here the cluster labels are always the positive integers,
        # we don't need a dictionary. `Vector` keys are the IDs we need.
        dist_to_clusters = [metric(feature, features[idx]) for idx in cluster_idxs]
        min_dist, closest_cluster_label = findmin(dist_to_clusters)

        if min_dist > threshold # bigger than threshold: new attractor
            # assign feature as new cluster center
            feature_label = next_cluster_label
            push!(cluster_idxs, idx_feature)
            push!(cluster_labels, next_cluster_label)
            next_cluster_label += 1
        else # smaller than threshold: assign to closest cluster
            feature_label = closest_cluster_label
        end

        labels[idx_feature] = feature_label
    end
    return labels
end
