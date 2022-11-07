export ClusteringAcrossParametersContinuation
import ProgressMeter
import Mmap

struct ClusteringAcrossParametersContinuation{A, M, I, E, C} <: BasinsFractionContinuation
    mapper::A
    info_extraction::E
    samples_per_parameter::I 
    par_weight::M
    mmap_limit::I
    cluster_config::C
end

"""
    ClusteringAcrossParametersContinuation(mapper::AttractorsViaFeaturizing; kwargs...)
 A method to compute the continuation of a basin when a parameter changes, see
  [`basins_fractions_continuation`](@ref). This method computes the Features accross a 
  range of parameters before performing a clustering into classes using a special metric, 
  see [^Gelbrecht2020]. 

The function takes as an input a `mapper` that maps initial conditions to attractors 
  using the featurizing method [^Stender2021]. See [`AttractorMapper`](@ref) for how to  
  use the `mapper`.

## Keyword Arguments
- `prange` Range of parameter to analyze.
- `pidx` Number or symbol of the parameter to change in the array of parameters of the dynamical system. 
- `ics` Sampler function to generate initial conditions to sample. 
- `samples_per_parameter = 100` Number of samples per parameter.
- `par_weight = 1` The distance matrix between features has a special extra  weight that 
  is proportional to the distance |pi - pj| between the parameters of each features. This 
  keyword argument is the weight coeficient that ponderates the distance matrix.
- `mmap_limit = 20000` this parameter sets the limit of features that should be clustered using only the 
  RAM memory. Above this limit the program uses a memory map based array to store the distance 
  matrix on the disk. 

## Description

The method first simulate and compute a set of statistics on the trajectories, for example the mean 
 and standard deviation of the trajectory data points. Once all the featured statistics have been 
 computed, the algorithm computes the distance between each feature. A special weitgh term is added 
 so that the distance between features with different parameters is increased. It helps to discriminate
 between attractors with very different parameters. At last, the distance matrix is clustered with the 
 DBSCAN algorithm

[^Gelbrecht2020]:
    Gelbrecht, M., Kurths, J., & Hellmann, F. (2020). Monte Carlo basin bifurcation
    analysis. New Journal of Physics, 22(3), 033032.

[^Stender2021]:
    Stender & Hoffmann, [bSTAB: an open-source software for computing the basin
    stability of multi-stable dynamical systems](https://doi.org/10.1007/s11071-021-06786-5)
"""
function ClusteringAcrossParametersContinuation(
        mapper::AttractorMapper;
        info_extraction = mean_across_features,
        samples_per_parameter = 100, 
        par_weight = 1, 
        mmap_limit = 20000,
        cluster_config = ClusterConfig()
    )
    return ClusteringAcrossParametersContinuation(mapper, info_extraction, samples_per_parameter, par_weight, mmap_limit, cluster_config)
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
        continuation::ClusteringAcrossParametersContinuation, prange, pidx, ics::Function;
        show_progress = true
    )
    (; mapper, info_extraction, samples_per_parameter, par_weight, mmap_limit, cluster_config) = continuation
    spp, n = samples_per_parameter, length(prange)

    features = _get_features_prange(mapper, ics, n, spp, prange, pidx, show_progress)

    # The distance matrix can get very large. The use of memory map based array is 
    # necessary.
    if length(features) > mmap_limit
        # Create temp file
        pth, s = mktemp()  
        dists = Mmap.mmap(s, Matrix{Float32}, ( length(features), length(features)))  
    else
        dists = zeros(length(features), length(features))
    end
    
    _get_dist_matrix!(features, dists, prange, spp, par_weight, cluster_config)

    cluster_labels = _cluster_across_parameters(dists, features, cluster_config)

    fractions_curves, attractors_info = _label_fractions(cluster_labels, n, spp, features[1], info_extraction)

    return fractions_curves, attractors_info
end
 

function _get_features_prange(mapper::AttractorsViaFeaturizing, ics, n, spp, prange, pidx, show_progress)
    progress = ProgressMeter.Progress(n;
        desc="Continuating basins fractions:", enabled=show_progress
    )
    # Extract the first possible feature to initialize the features container
    feature = extract_features(mapper, ics; N = 1)
    features = Vector{typeof(feature[1])}(undef, n*spp)
    # Collect features
    for (i, p) in enumerate(prange)
        set_parameter!(mapper.integ, pidx, p)
        current_features = extract_features_threaded(mapper, ics; show_progress, N = spp)
        features[((i - 1)*spp + 1):i*spp] .= current_features
        next!(progress)
    end
    return features
end


function _get_dist_matrix!(features, dists, prange, spp, par_weight, cluster_config)
    # Construct distance matrix
    metric = cluster_config.clust_distance_metric
    pairwise!(metric, dists, features; symmetric = true)
    # use parameter distance weight (w is the weight for one parameter only)
    # Parameter range is rescaled from 0 to 1.
    par_array = kron(range(0,1,length(prange)), ones(spp))
    for k in 1:length(par_array)
        for j in 1:length(par_array)
            dists[k,j] += par_weight*abs(par_array[k]-par_array[j])
        end
    end
end


function _cluster_across_parameters(dists, features, cluster_config)
    # Cluster the values accross parameters
    cc = cluster_config
    ftrs = reduce(hcat, features) # Convert to Matrix from Vector{Vector}
    cluster_labels = cluster_features_clustering(ftrs, cc.min_neighbors, cc.clust_distance_metric, 
        false, cc.optimal_radius_method, cc.num_attempts_radius, cc.silhouette_statistic, 
        cc.max_used_features; dists
    )
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

function _get_features_prange(mapper::AttractorsViaRecurrences, ics, n, spp, prange, pidx, show_progress)
    progress = ProgressMeter.Progress(n;
        desc="Continuating basins fractions:", enabled=show_progress
    )

    set_parameter!(mapper.integ, pidx, prange[1])
    fs = basins_fractions(mapper, ics; show_progress = true, N = spp)
    fractions_curves = [fs]
    current_attractors = deepcopy(mapper.bsn_nfo.attractors)
    attractors_info = [current_attractors]

    for (i, p) in enumerate(prange)
        set_parameter!(mapper.integ, pidx, p)
        reset!(mapper)
        fs = basins_fractions(mapper, ics; show_progress = true, N = spp)
        push!(fractions_curves, fs)
        push!(attractors_info, deepcopy(mapper.bsn_nfo.attractors))
        next!(progress)
    end
    # collect Datasets
    vec_att = Dataset[]
    for att in attractors_info
        for a in att
            push!(vec_att, a[2])
        end
    end 
    return vec_att
end
