export ClusteringAcrossParametersContinuation
import ProgressMeter
import Mmap

struct ClusteringAcrossParametersContinuation{A, E} <: BasinsFractionContinuation
    mapper::A
    info_extraction::E
    par_weight::Float64
    use_mmap::Bool
end

"""
    ClusteringAcrossParametersContinuation(mapper::AttractorsViaFeaturizing; kwargs...)
A method for [`basins_fractions_continuation`](@ref).
It uses clustering across of features across a parameter range, potentially weighted
by the distance in parameter space. Its input `mapper` must have
a `GroupByClustering` as its grouping configuration.
This is the original continuation method described in [^Gelbrecht2020].

## Keyword Arguments
- `info_extraction::Function` a function that takes as an input a vector of features
  (corresponding to a cluster) and returns a description of the cluster.
  By default, the centroid of the cluster is used.
- `par_weight = 0` The distance matrix between features has a special extra weight that
  is proportional to the distance `|p[i] - p[j]|` between the parameters used when
  extracting features. This keyword argument is the weight coeficient that ponderates
  the distance matrix. Notice that the range of parameters is normalized from 0 to 1
  such that the largest distance in the parameter space is 1. The normalization is done
  because the feature space is also (by default) normalized to 0-1.
- `use_mmap = false` this parameter whether the feature distance matrix should be computed
  in memory or on the disk using memory map. Should be used if a matrix with side
  `length(prange)*samples_per_parameter` exceeds available memory.

## Description

The method first integrates and computes a set of statistics on the trajectories, for example
the mean and standard deviation of the trajectory data points. Once all the featured statistics
have been  computed, the algorithm computes the distance between each feature. A special weight
term is added so that the distance between features with different parameters is increased.
It helps to discriminate between attractors with very different parameters. At last, the
distance matrix is clustered with the DBSCAN algorithm.

[^Gelbrecht2021]:
    Maximilian Gelbrecht et al 2021, Monte Carlo basin bifurcation analysis,
    [New J. Phys.22 03303](http://dx.doi.org/10.1088/1367-2630/ab7a05)
"""
function ClusteringAcrossParametersContinuation(
        mapper::AttractorsViaFeaturizing;
        info_extraction = mean_across_features,
        par_weight = 0.0,
        use_mmap = false,
    )
    if !(mapper.group_config isa GroupViaClustering)
        throw(ArgumentError(
            "This method needs `GroupViaClustering` as the grouping configuration."
        ))
    end
    return ClusteringAcrossParametersContinuation(
        mapper, info_extraction, par_weight, use_mmap
    )
end

function mean_across_features(fs)
    means = zeros(length(first(fs)))
    N = length(fs)
    for f in fs
        for i in eachindex(f)
            means[i] += f[i]
        end
    end
    return means ./ N
end

function basins_fractions_continuation(
        continuation::ClusteringAcrossParametersContinuation, prange, pidx, ics;
        show_progress = true, samples_per_parameter = 100
    )
    (; mapper, info_extraction, par_weight) = continuation
    spp, n = samples_per_parameter, length(prange)
    features = _get_features_prange(mapper, ics, n, spp, prange, pidx, show_progress)
    cluster_labels = cluster_all_features(
        features, mapper.group_config, par_weight;
        prange, use_mmap = continuation.use_mmap, show_progress,
        samples_per_parameter
    )
    fractions_curves, attractors_info = _label_fractions(cluster_labels, n, spp, features[1], info_extraction)
    return fractions_curves, attractors_info
end

# TODO: This whole business makes me feel like we should
# just make the `mmap` and then `par_weight` parameters part of `GroupViaClustering`,
# and then the below function is a literaly call to `group_features`...
# Although it does become a pain in the ass with all this keeping track of the
# samples per parameter nad the length of the paramter space... GRRRR
# But I guess we can define a single 5 argument method that calls the 2 argument
# metod only for the `GroupViaClustering` type.
function cluster_all_features(
        features, clust_config, par_weight;
        prange = nothing, use_mmap = false, show_progress = true,
        samples_per_parameter = 100, cluster_in_slice = 0. 
    )
    # The distance matrix can get very large so we allow  memory map based array
    if use_mmap
        pth, s = mktemp()
        dists = Mmap.mmap(s, Matrix{Float32}, (length(features), length(features)))
    else
        dists = zeros(length(features), length(features))
    end
    progress = ProgressMeter.ProgressUnknown(desc="Clustering", enabled=show_progress)
    ProgressMeter.next!(progress)
    _get_dist_matrix!(features, dists, prange, samples_per_parameter, par_weight, clust_config, cluster_in_slice)
    ProgressMeter.next!(progress)
    cluster_labels = _cluster_across_parameters(dists, features, clust_config)
    ProgressMeter.finish!(progress)
    return cluster_labels
end


function _get_features_prange(mapper::AttractorsViaFeaturizing, ics, n, spp, prange, pidx, show_progress)
    progress = ProgressMeter.Progress(n;
        desc="Generating features", enabled=show_progress
    )
    # Extract the first possible feature to initialize the features container
    feature = extract_features(mapper, ics; N = 1)
    features = Vector{typeof(feature[1])}(undef, n*spp)
    # Collect features
    for (i, p) in enumerate(prange)
        set_parameter!(mapper.integ, pidx, p)
        current_features = extract_features(mapper, ics; show_progress, N = spp)
        features[((i - 1)*spp + 1):i*spp] .= current_features
        ProgressMeter.next!(progress)
    end
    return features
end


function _get_dist_matrix!(features, dists, prange, spp, par_weight, clust_config, cluster_in_slice)
    # Construct distance matrix
    metric = clust_config.clust_distance_metric
    pairwise!(metric, dists, features; symmetric = true)
    # use parameter distance weight (w is the weight for one parameter only)
    # Parameter range is normalized from 0 to 1.
    if par_weight ≠ 0
        isnothing(prange) && error("`par_weight` isn't 0, but `prange` isn't given!")
        par_weight = par_weight/abs(prange[end]-prange[1])
        par_array = kron(range(0,1,length(prange)), ones(spp))
        @inbounds for k in eachindex(par_array)
            for j in eachindex(par_array)
                dists[k,j] += par_weight*abs(par_array[k] - par_array[j])
                # Add a weight to the dist if we do not want to cluster for the same parameter slice
                if par_array[k] == par_array[j]
                    dists[k,j] += cluster_in_slice
                end
            end
        end
    end
    @show dists
end

function _cluster_across_parameters(dists, features, cc::GroupViaClustering)
    ϵ_optimal =  _extract_ϵ_optimal(features, cc)
    cluster_labels = _cluster_distances_into_labels(dists, ϵ_optimal, cc.min_neighbors)
    return cluster_labels
end

function _label_fractions(cluster_labels, n, spp, feature, info_extraction)
    # finally we collect/group stuff into their dictionaries
    fractions_curves = Vector{Dict{Int, Float64}}(undef, n)
    dummy_info = info_extraction(feature)
    attractors_info = Vector{Dict{Int, typeof(dummy_info)}}(undef, n)
    for i in 1:n
        # Here we know which indices correspond to which parameter value
        # because they are sequentially increased every `spp`
        # (steps per parameter)
        current_labels = view(cluster_labels, ((i - 1)*spp + 1):i*spp)
        current_ids = unique(current_labels)
        # getting fractions is easy; use API function that takes in arrays
        fractions_curves[i] = basins_fractions(current_labels, current_ids)
        attractors_info[i] = Dict(id => info_extraction(
            view(current_labels, findall(isequal(id), current_labels)))
            for id in current_ids
        )
    end
    return fractions_curves, attractors_info
end
