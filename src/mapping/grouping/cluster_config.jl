using Distances, Clustering, Distributions
export ClusteringConfig, cluster_features

"""
    ClusteringConfig(; kwargs...)

Initialize a struct that contains information used to cluster "features".
These features are typically extracted from trajectories/datasets in
[`AttractorsViaFeaturizing`](@ref), or manualy by the user.

The clustering is done in the function [`cluster_features`](@ref).

The default clustering method is an improvement over existing literature, see Description.

## Keyword arguments
* `clust_distance_metric = Euclidean()`: metric to be used in the clustering. Used in both versions.
### Supervised method
* `templates = nothing`: Enables supervised version, see below. If given (i.e. different
  than `nothing`), `templates` must be a `Dict` of cluster labels to cluster features. The
  labels must be of `Int` type, and the features are `Vector`s representing a cluster
  (which can be an attractor, for instance). The label `-1` is reserved for invalid
  trajectories, which either diverge or whose clustering failed.
* `max_distance_template = Inf`: Maximum allowed distance between a feature and its nearest
  template for it to be assigned to that template. By default, `Inf` guarantees that a
  feature is assigned to its nearest template regardless of the distance.
### Unsupervised method
* `rescale_features = true`: if true, rescale each dimension of the extracted features
  separately into the range `[0,1]`. This typically leads to more accurate clustering.
* `min_neighbors = 10`: minimum number of neighbors (i.e. of similar features) each
  feature needs to have in order to be considered in a cluster (fewer than this, it is
  labeled as an outlier, `-1`).
* `optimal_radius_method::Union{Real, String} = "silhouettes_optim"`: if a real number, it
  is the radius used to cluster features in the unsupervised method. Otherwise, it
  determines the method used to automatically determine that radius. Possible
  values are:
  - `"silhouettes"`: Performs a linear (sequential) search for the radius that maximizes a
    statistic of the silhouette values of clusters (typically the mean). This can be chosen
    with `silhouette_statistic`. The linear search may take some time to finish. To
    increase speed, the number of radii iterated through can be reduced by decreasing
    `num_attempts_radius` (see its entry below).
  - `"silhouettes_optim"`: Same as `"silhouettes"` but performs an optimized search via
    Optim.jl. It's faster than `"silhouettes"`, with typically the same accuracy (the
    search here is not guaranteed to always find the global maximum, though it typically
    gets close).
  - `"knee"`: chooses the the radius according to the knee (a.k.a. elbow,
    highest-derivative method) and is quicker, though generally leading to much worse
    clustering. It requires that `min_neighbors` > 1.
#### Keywords for optimal radius estimation
* `num_attempts_radius = 100`: number of radii that the `optimal_radius_method` will try
  out in its iterative procedure. Higher values increase the accuracy of clustering,
  though not necessarily much, while always reducing speed.
* `silhouette_statistic::Function = mean`: statistic (e.g. mean or minimum) of the
  silhouettes that is maximized in the "optimal" clustering. The original implementation
  in [^Stender2021] used the `minimum` of the silhouettes, and typically performs less
  accurately than the `mean`.
* `max_used_features = 0`: if not `0`, it should be an `Int` denoting the max amount of
  features to be used when finding the optimal radius. Useful when clustering a very large
  number of features (e.g., high accuracy estimation of fractions of basins of
  attraction).

## Description
The trajectory `X`, which may for instance be an attractor, is transformed into a vector
of features. Each feature is a number useful in _characterizing the attractor_, and
distinguishing it from other attrators. Example features are the mean or standard
deviation of one of the of the timeseries of the trajectory, the entropy of the first two
dimensions, the fractal dimension of `X`, or anything else you may fancy. The vectors of
features are then used to identify clusters of the trajectories. An application thus is to
cluster attractors. Currently two methods are offered to achieve this.

### Unsupervised method
The **unsupervised method** does not rely on templates, and instead uses the DBSCAN
clustering algorithm to automatically identify clusters of similar features. To achieve
this, each feature is considered a point in feature space. Each cluster then basically
groups points that are closely packed together. Closely packed means that the points have
at least `min_neighbors` inside a ball of radius `optimal_radius` centered on them. This
method typically works well if the radius is chosen well, which is not necessarily an easy
task. Currently, three methods are implemented to automatically estimate an "optimal"
radius.
#### Estimating the optimal radius
The default method is the **silhouettes method**, which includes keywords `silhouette` and
`silhouette_optim`. Both of them search for the radius that optimizes the clustering,
meaning the one that maximizes a statistic `silhouette_statistic` (e.g. mean value) of a
quantifier for the quality of each cluster. This quantifier is the silhouette value of
each identified cluster. A silhouette value measures how similar a point is to the cluster
it currently belongs to, compared to the other clusters, and ranges from -1 (worst
matching) to +1 (ideal matching). If only one cluster is found, the assigned silhouette is
0. So for each attempted radius in the search the clusters are computed, their silhouettes
calculated, and the statistic of these silhouettes computed. The algorithm then finds the
radius that leads to the maximum such statistic. For `optimal_radius_method =
"silhouettes"`, the search is done linearly, from a minimum to a maximum candidate radius
for `optimal_radius_method = "silhouettes"`; `optimal_radius_method = silhouettes_optim`,
it is done via an optimized search performed by Optim.jl which is typically faster and
with similar accuracy. A third alternative is the`"elbow"` method, which works by
calculating the distance of each point to its k-nearest-neighbors (with `k=min_neighbors`)
and finding the distance corresponding to the highest derivative in the curve of the
distances, sorted in ascending order. This distance is chosen as the optimal radius. It is
described in [^Kriegel1996] and [^Schubert2017]. It typically performs considerably worse
than the `"silhouette"` methods.

### Supervised method
In the alternative, **supervised version**, the user provides features to be used as
templates guiding the clustering via the `templates` keyword. Features are clustered to
their nearest template if the distance to that template is below or equal to the threshold
`max_distance_template`. The distances are computed via the metric defined by
`clust_distance_metric`. So the algorithm first finds the template nearest to each feature
(in that metric) via the `k`-nearest-neighbors algorithm for `k=1` and then compares their
distance to the threshold. If the distance is below or equal to the threshold, the feature
is assigned the template's label, which is given in `templates`. If the distance is
larger, the feature is assigned label `-1`.

[^Stender2021]:
    Stender & Hoffmann, [bSTAB: an open-source software for computing the basin
    stability of multi-stable dynamical systems](https://doi.org/10.1007/s11071-021-06786-5)
[^Kriegel1996]: Ester, Kriegel, Sander and Xu: A Density-Based Algorithm for Discovering
    Clusters in Large Spatial Databases with Noise
[^Schubert2017]:
    Schubert, Sander, Ester, Kriegel and Xu: DBSCAN Revisited, Revisited: Why and How You
    Should (Still) Use DBSCAN
"""
mutable struct ClusteringConfig{A, M}
    templates::A
    clust_distance_metric::M
    max_distance_template::Float64
    min_neighbors::Int
    rescale_features::Bool
    optimal_radius_method::Union{Real, String}
    num_attempts_radius::Int
    silhouette_statistic::Function
    max_used_features::Int
end

function ClusteringConfig(; templates::Union{Nothing, Dict} = nothing,
        clust_distance_metric=Euclidean(), max_distance_template = Inf, min_neighbors = 10,
        rescale_features=true, optimal_radius_method::Union{Real, String}="silhouettes_optim",
        num_attempts_radius=100, silhouette_statistic = mean, max_used_features = 0,
    )
    return ClusteringConfig(
        templates, clust_distance_metric,
        Float64(max_distance_template), min_neighbors,
        rescale_features, optimal_radius_method,
        num_attempts_radius, silhouette_statistic, max_used_features
    )
end

include("cluster_utils.jl")

#####################################################################################
# Clustering classification functions
#####################################################################################
"""
    cluster_features(features, cc::ClusteringConfig)
Cluster the given `features::Vector{<:AbstractVector}`, according to given
[`ClusteringConfig`](@ref). Return `cluster_labels`, which contains, for each feature, the
labels (indices) of the corresponding cluster.
"""
function cluster_features(features::Vector{<:AbstractVector}, cc::ClusteringConfig)
    # All methods require the features in a matrix format
    f = reduce(hcat, features) # Convert to Matrix from Vector{Vector}
    f = float.(f)
    if !isnothing(cc.templates)
        cluster_features_templates(f, cc)
    else
        cluster_features_clustering(
            f, cc.min_neighbors, cc.clust_distance_metric,
            cc.rescale_features, cc.optimal_radius_method,
            cc.num_attempts_radius, cc.silhouette_statistic, cc.max_used_features,
        )
    end
end

# Supervised method: for each template, return the label of its nearest template (cluster)
# if the distance is smaller than the threshold.
function cluster_features_templates(features, cc::ClusteringConfig)
    #puts each vector into a column, with the ordering based on the order given in keys(d)
    templates = float.(reduce(hcat, [cc.templates[i] for i ∈ keys(cc.templates)]))

    #prepare for nearest-neighbors algorithm (kNN with k=1)
    template_tree = searchstructure(KDTree, templates, cc.clust_distance_metric)
    cluster_labels, dists_to_templates = bulksearch(template_tree, features, NeighborNumber(1))
    cluster_labels = reduce(vcat, cluster_labels) # make it a vector
    dists_to_templates = reduce(vcat, dists_to_templates)

    # Make label -1 if error strictly bigger than threshold
    cluster_labels[dists_to_templates .> cc.max_distance_template] .= -1
    matlabels_to_dictlabels = Dict(1:length(keys(cc.templates)) .=> keys(cc.templates))

    # cluster_user_labels[i] is the label of template nearest to feature i
    # (i-th column of features matrix)
    cluster_user_labels = replace(cluster_labels, matlabels_to_dictlabels...)
    return cluster_user_labels
end

"""
Do "min-max" rescaling of vector `vec`: rescale it such that its values span `[0,1]`.
"""
function _rescale!(vec::Vector{T}) where T
    vec .-= minimum(vec)
    max = maximum(vec)
    if max == 0 return zeros(T, length(vec)) end
    vec ./= maximum(vec)
end

# Unsupervised method: clustering in feature space
function cluster_features_clustering(
        features, min_neighbors, metric, rescale_features, optimal_radius_method,
        num_attempts_radius, silhouette_statistic, max_used_features
    )
    # needed because dbscan, as implemented, needs to receive as input a matrix D x N
    # such that D < N
    dimfeats, nfeats = size(features)
    if dimfeats ≥ nfeats @warn "Not enough features. The algorithm needs the number of features
         $nfeats to be greater or equal than the number of dimensions $dimfeats";
           return 1:nfeats, zeros(nfeats) end

    if rescale_features
        features = mapslices(_rescale!, features; dims=2)
    end

    dists = pairwise(metric, features)

    # These functions are called from cluster_utils.jl
    if optimal_radius_method isa String
      features_for_optimal = if max_used_features == 0
          features
      else
          StatsBase.sample(features, minimum(length(features), max_used_features); replace = false)
      end
      ϵ_optimal = optimal_radius_dbscan(
          features_for_optimal, min_neighbors, metric, optimal_radius_method,
          num_attempts_radius, silhouette_statistic
      )
    elseif optimal_radius_method isa Real
      ϵ_optimal = optimal_radius_method
    else
      error("Specified optimal_radius_method is incorrect. Please specify the radius
       directly as a Real number or the method to compute it as a String")
    end

    dbscanresult = dbscan(dists, ϵ_optimal, min_neighbors)
    cluster_labels = cluster_assignment(dbscanresult)
    return cluster_labels
end
