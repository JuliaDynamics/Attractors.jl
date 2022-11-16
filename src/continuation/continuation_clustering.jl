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

    # The distance matrix can get very large. The use of memory map based array is
    # necessary.
    if continuation.use_mmap
        # Create temp file
        pth, s = mktemp()
        dists = Mmap.mmap(s, Matrix{Float32}, (length(features), length(features)))
    else
        dists = zeros(length(features), length(features))
    end

    progress = ProgressMeter.ProgressUnknown(;
        desc="Clustering: ", enabled=show_progress
    )
    _get_dist_matrix!(features, dists, prange, spp, par_weight, mapper)
    cluster_labels = _cluster_across_parameters(dists, features, mapper)
    ProgressMeter.finish!(progress)

    fractions_curves, attractors_info = _label_fractions(cluster_labels, n, spp, features[1], info_extraction)

    return fractions_curves, attractors_info
end


function _get_features_prange(mapper::AttractorsViaFeaturizing, ics, n, spp, prange, pidx, show_progress)
    progress = ProgressMeter.Progress(n;
        desc="Generating features: ", enabled=show_progress
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


function _get_dist_matrix!(features, dists, prange, spp, par_weight,
    mapper::AttractorsViaFeaturizing)
    # Construct distance matrix
    metric = mapper.group_config.clust_distance_metric
    pairwise!(metric, dists, features; symmetric = true)
    # use parameter distance weight (w is the weight for one parameter only)
    # Parameter range is normalized from 0 to 1.
    par_array = kron(range(0,1,length(prange)), ones(spp))
    @inbounds for k in eachindex(par_array)
        for j in eachindex(par_array)
            dists[k,j] += par_weight*abs(par_array[k] - par_array[j])
        end
    end
end

function _cluster_across_parameters(dists, features, mapper::AttractorsViaFeaturizing)
    # Cluster the values accross parameters
    cc = mapper.group_config
    ϵ_optimal =  _extract_ϵ_optimal(features, cc)
    cluster_labels = _cluster_distances_into_labels(dists, ϵ_optimal, cc.min_neighbors)
    return cluster_labels
end

function _label_fractions(cluster_labels, n, spp, feature, info_extraction)
    # And finally collect/group stuff into their dictionaries
    fractions_curves = Vector{Dict{Int, Float64}}(undef, n)
    dummy_info = info_extraction(feature)
    attractors_info = Vector{Dict{Int, typeof(dummy_info)}}(undef, n)
    for i in 1:n
        current_labels = view(cluster_labels, ((i - 1)*spp + 1):i*spp)
        current_ids = unique(current_labels)
        # getting fractions is easy; use predefined function
        fractions_curves[i] = basins_fractions(current_labels, current_ids)
        attractors_info[i] = Dict(id => info_extraction(
            view(current_labels, findall(isequal(id), current_labels)))
            for id in current_ids
        )
    end
    return fractions_curves, attractors_info
end

