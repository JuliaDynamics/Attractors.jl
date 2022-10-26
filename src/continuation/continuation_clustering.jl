export ClusteringAcrossParametersContinuation
import ProgressMeter

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
- `samples_per_parameter` Number of samples per parameter.
- `par_weight` The distance matrix between features has a special extra  weight that is 
 proportional to the distance |p_i - p_j| between the parameters of each features. This 
 keyword argument is the weight coeficient that ponderates the distance matrix.

## Description

[^Gelbrecht2020]:
    Gelbrecht, M., Kurths, J., & Hellmann, F. (2020). Monte Carlo basin bifurcation
    analysis. New Journal of Physics, 22(3), 033032.

[^Stender2021]:
    Stender & Hoffmann, [bSTAB: an open-source software for computing the basin
    stability of multi-stable dynamical systems](https://doi.org/10.1007/s11071-021-06786-5)
"""
function ClusteringAcrossParametersContinuation(
        mapper::AttractorsViaFeaturizing;
        info_extraction = mean_across_features
        # TODO: Here we can add more keywords regarding how to cluster across parameters.
    )
    return (; mapper, info_extraction)
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
        continuation::NamedTuple, prange, pidx, ics::Function;
        samples_per_parameter = 100, show_progress = true, par_weight = 1
    )
    spp, n = samples_per_parameter, length(prange)
    (; mapper, info_extraction) = continuation
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

    # Construct distance matrix
    cc = mapper.cluster_config
    metric = cc.clust_method_norm
    dists = pairwise(metric, features)

    # use parameter distance weight (w is the weight for one parameter only)
    # Parameter range is rescaled from 0 to 1.
    par_array = kron(range(0,1,length(prange)), ones(spp))
    for k in 1:length(par_array)
        for j in 1:length(par_array)
            dists[k,j] += par_weight*metric(par_array[k],par_array[j])
        end
    end

    cluster_labels =  _cluster_features_across_parameters(features, dists, cc)

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


function _cluster_features_across_parameters(features, dists, cc)
    # Cluster them
    metric = cc.clust_method_norm
    # @show features
    f = reduce(hcat, features) # Convert to Matrix from Vector{Vector}
    f = float.(f)
    features_for_optimal = if cc.max_used_features == 0
        f
    else
        StatsBase.sample(f, minimum(length(features), cc.max_used_features); replace = false)
    end 
    ϵ_optimal = optimal_radius_dbscan(
        features_for_optimal, cc.min_neighbors, metric, cc.optimal_radius_method,
        cc.num_attempts_radius, cc.silhouette_statistic
    )
    if ϵ_optimal > 0 
        dbscanresult = dbscan(dists, ϵ_optimal, cc.min_neighbors)
    else
        @warn "Optimal radius ϵ_optimal is 0, using 5 instead" 
        dbscanresult = dbscan(dists, 5, cc.min_neighbors)
    end
    cluster_labels = cluster_assignment(dbscanresult)
    return cluster_labels
end
