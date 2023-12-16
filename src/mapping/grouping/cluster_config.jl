using Distances, Clustering, Distributions
import Optim
export GroupViaClustering

"""
    GroupViaClustering(; kwargs...)

Initialize a struct that contains instructions on how to group features in
[`AttractorsViaFeaturizing`](@ref). `GroupViaClustering` clusters features into
groups using DBSCAN, similar to the original work by bSTAB [Stender2021](@cite) and
MCBB [Gelbrecht2020](@cite). Several options on clustering are available, see keywords below.

The defaults are a significant improvement over existing literature, see Description.

## Keyword arguments

* `clust_distance_metric = Euclidean()`: A metric to be used in the clustering.
  It can be any function `f(a, b)` that returns the distance between real-valued vectors
  `a, b`. All metrics from Distances.jl can be used here.
* `rescale_features = true`: if true, rescale each dimension of the extracted features
  separately into the range `[0,1]`. This typically leads to more accurate clustering.
* `min_neighbors = 10`: minimum number of neighbors (i.e. of similar features) each
  feature needs to have, including counting its own self,
  in order to be considered in a cluster (fewer than this, it is
  labeled as an outlier, `-1`).
* `use_mmap = false`: whether to use an on-disk map for creating the distance matrix
  of the features. Useful when the features are so many where a matrix with side their
  length would not fit to memory.

### Keywords for optimal radius estimation

* `optimal_radius_method::Union{Real, String} = "silhouettes_optim"`: if a real number, it
  is the radius used to cluster features. Otherwise, it determines the method used to
  automatically determine that radius. Possible values are:
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
* `num_attempts_radius = 100`: number of radii that the `optimal_radius_method` will try
  out in its iterative procedure. Higher values increase the accuracy of clustering,
  though not necessarily much, while always reducing speed.
* `silhouette_statistic::Function = mean`: statistic (e.g. mean or minimum) of the
  silhouettes that is maximized in the "optimal" clustering. The original implementation
  in [Stender2021](@cite) used the `minimum` of the silhouettes, and typically performs less
  accurately than the `mean`.
* `max_used_features = 0`: if not `0`, it should be an `Int` denoting the max amount of
  features to be used when finding the optimal radius. Useful when clustering a very large
  number of features (e.g., high accuracy estimation of fractions of basins of
  attraction).

## Description

The DBSCAN clustering algorithm is used to automatically identify clusters of similar
features. Each feature vector is a point in a feature space. Each cluster then basically
groups points that are closely packed together. Closely packed means that the points have
at least `min_neighbors` inside a ball of radius `optimal_radius` centered on them. This
method typically works well if the radius is chosen well, which is not necessarily an easy
task. Currently, three methods are implemented to automatically estimate an "optimal"
radius.

### Estimating the optimal radius

The default method is the **silhouettes method**, which includes keywords `silhouette` and
`silhouette_optim`. Both of them search for the radius that optimizes the clustering,
meaning the one that maximizes a statistic `silhouette_statistic` (e.g. mean value) of a
quantifier for the quality of each cluster. This quantifier is the silhouette value of
each identified cluster. A silhouette value measures how similar a point is to the cluster
it currently belongs to, compared to the other clusters, and ranges from -1 (worst
matching) to +1 (ideal matching). If only one cluster is found, the assigned silhouette is
zero. So for each attempted radius in the search the clusters are computed, their silhouettes
calculated, and the statistic of these silhouettes computed. The algorithm then finds the
radius that leads to the maximum such statistic. For `optimal_radius_method =
"silhouettes"`, the search is done linearly, from a minimum to a maximum candidate radius
for `optimal_radius_method = "silhouettes"`; `optimal_radius_method = silhouettes_optim`,
it is done via an optimized search performed by Optim.jl which is typically faster and
with similar accuracy. A third alternative is the`"elbow"` method, which works by
calculating the distance of each point to its k-nearest-neighbors (with `k=min_neighbors`)
and finding the distance corresponding to the highest derivative in the curve of the
distances, sorted in ascending order. This distance is chosen as the optimal radius. It is
described in [Ester1996](@cite) and [Schubert2017](@cite). It typically performs considerably worse
than the `"silhouette"` methods.
"""
struct GroupViaClustering{R<:Union{Real, String}, M, F<:Function} <: GroupingConfig
    clust_distance_metric::M
    min_neighbors::Int
    rescale_features::Bool
    optimal_radius_method::R
    num_attempts_radius::Int
    silhouette_statistic::F
    max_used_features::Int
    use_mmap::Bool
end

function GroupViaClustering(;
        clust_distance_metric=Euclidean(), min_neighbors = 10,
        rescale_features=true, num_attempts_radius=100,
        optimal_radius_method::Union{Real, String} = "silhouettes_optim",
        silhouette_statistic = mean, max_used_features = 0,
        use_mmap = false,
    )
    return GroupViaClustering(
        clust_distance_metric, min_neighbors,
        rescale_features, optimal_radius_method,
        num_attempts_radius, silhouette_statistic, max_used_features,
        use_mmap,
    )
end

#####################################################################################
# API function (group features)
#####################################################################################
# The keywords `par_weight, plength, spp` enable the "for-free" implementation of the
# MCBB algorithm (weighting the distance matrix by parameter value as well).
# The keyword version of this function is only called in
# `GroupingAcrossParametersContinuation` and is not part of public API!
function group_features(
        features, config::GroupViaClustering;
        par_weight::Real = 0, plength::Int = 1, spp::Int = 1,
    )
    nfeats = length(features); dimfeats = length(features[1])
    if dimfeats ≥ nfeats
        throw(ArgumentError("""
        Not enough features. The algorithm needs the number of features
        $nfeats to be greater or equal than the number of dimensions $dimfeats
        """))
    end
    if config.rescale_features
        features = _rescale_to_01(features)
    end
    ϵ_optimal = _extract_ϵ_optimal(features, config)
    distances = _distance_matrix(features, config; par_weight, plength, spp)
    labels = _cluster_distances_into_labels(distances, ϵ_optimal, config.min_neighbors)
    return labels
end

function _distance_matrix(features, config::GroupViaClustering;
        par_weight::Real = 0, plength::Int = 1, spp::Int = 1
    )
    metric = config.clust_distance_metric
    L = length(features)
    if config.use_mmap
        pth, s = mktemp()
        dists = Mmap.mmap(s, Matrix{Float32}, (L, L))
    else
        dists = zeros(L, L)
    end
    if metric isa Metric # then the `pairwise` function is valid
        # ensure that we give the vector of static vectors to pairwise!
        pairwise!(metric, dists, vec(features); symmetric = true)
    else # it is any arbitrary distance function, e.g., used in aggregating attractors
        @inbounds for i in eachindex(features)
            Threads.@threads for j in i:length(features)
                v = metric(features[i], features[j])
                dists[i, j] = dists[j, i] = v # utilize symmetry
            end
        end
    end

    if par_weight ≠ 0 # weight distance matrix by parameter value
        par_vector = kron(range(0, 1, plength), ones(spp))
        length(par_vector) ≠ size(dists, 1) && error("Feature size doesn't match.")
        @inbounds for k in eachindex(par_vector)
            # We can optimize the loop here due to symmetry of the metric.
            # Instead of going over all `j` we go over `(k+1)` to end,
            # and also add value to transpose. (also assume that if j=k, distance is 0)
            for j in (k+1):size(dists, 1)
                pdist = par_weight*abs(par_vector[k] - par_vector[j])
                dists[k,j] += pdist
                dists[j,k] += pdist
            end
        end
    end
    return dists
end

function _extract_ϵ_optimal(features, config::GroupViaClustering)
    (; min_neighbors, clust_distance_metric, optimal_radius_method,
    num_attempts_radius, silhouette_statistic, max_used_features) = config

    if optimal_radius_method isa String
        # subsample features to accelerate optimal radius search
        if max_used_features == 0 || max_used_features > length(features)
            features_for_optimal = vec(features)
        else
            features_for_optimal = sample(vec(features), max_used_features; replace = false)
        end
        # get optimal radius (function dispatches on the radius method)
        ϵ_optimal, v_optimal = optimal_radius_dbscan(
            features_for_optimal, min_neighbors, clust_distance_metric,
            optimal_radius_method, num_attempts_radius, silhouette_statistic
        )
    elseif optimal_radius_method isa Real
        ϵ_optimal = optimal_radius_method
    else
        error("Specified `optimal_radius_method` is incorrect. Please specify the radius
        directly as a Real number or the method to compute it as a String")
    end
    return ϵ_optimal
end

# Already expecting the distance matrix, the output of `pairwise`
function _cluster_distances_into_labels(distances, ϵ_optimal, min_neighbors)
    dbscanresult = dbscan(distances, ϵ_optimal; min_neighbors, metric=nothing)
    cluster_labels = cluster_assignment(dbscanresult)
    return cluster_labels
end

"""
Do "min-max" rescaling of vector of feature vectors so that its values span `[0,1]`.
"""
_rescale_to_01(features::Vector{<:AbstractVector}) = _rescale_to_01(StateSpaceSet(features))
function _rescale_to_01(features::AbstractStateSpaceSet)
    mini, maxi = minmaxima(features)
    rescaled = map(f -> (f .- mini) ./ (maxi .- mini), features)
    return typeof(features)(rescaled) # ensure it stays the same type
end

#####################################################################################
# Utilities
#####################################################################################
"""
Util function for `classify_features`. Returns the assignment vector, in which the i-th
component is the cluster index of the i-th feature.
"""
function cluster_assignment(clusters, data; include_boundary=true)
    assign = zeros(Int, size(data)[2])
    for (idx, cluster) in enumerate(clusters)
        assign[cluster.core_indices] .= idx
        if cluster.boundary_indices != []
            if include_boundary
                assign[cluster.boundary_indices] .= idx
            else
                assign[cluster.boundary_indices] .= -1
            end
        end
    end
    return assign
end
function cluster_assignment(dbscanresult::Clustering.DbscanResult)
    labels = dbscanresult.assignments
    return replace!(labels, 0=>-1)
end

#####################################################################################
# Finding optimal ϵ
#####################################################################################
include("cluster_optimal_ϵ.jl")
