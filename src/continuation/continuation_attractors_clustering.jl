export ClusteringAttractorsContinuation
import ProgressMeter
import Mmap

struct ClusteringAttractorsContinuation{A, G, F} <: BasinsFractionContinuation
    mapper::A
    featurizer::F
    cluster_in_slice::Float64
    par_weight::Float64
    group_config::G
end

"""
    ClusteringAttractorsContinuation(mapper::AttractorsViaFeaturizing; kwargs...)
A method for [`basins_fractions_continuation`](@ref).
It uses clustering across of features across a parameter range, potentially weighted
by the distance in parameter space. Its input `mapper` must have
a `GroupByClustering` as its grouping configuration.

## Keyword Arguments
## Description

"""
function ClusteringAttractorsContinuation(
        mapper::AttractorsViaRecurrences;
        featurizer = identity,
        cluster_in_slice = 100., 
        par_weight = 0.0,
        group_config = GroupViaClustering()
    )
    return ClusteringAttractorsContinuation(
        mapper, featurizer, cluster_in_slice, par_weight, group_config
    )
end


function basins_fractions_continuation(
        continuation::ClusteringAttractorsContinuation, prange, pidx, ics;
        show_progress = true, samples_per_parameter = 100
    )
    (; mapper, featurizer, cluster_in_slice, par_weight, group_config) = continuation
    spp, n = samples_per_parameter, length(prange)

    features, pindex, f_curves, att_info = _get_attractors_prange(mapper, ics, n, spp, prange, pidx, featurizer, show_progress)

    cluster_labels = cluster_all_features(features, group_config, par_weight; prange = pindex, samples_per_parameter = 1, cluster_in_slice)

    j_curves = _label_fractions(cluster_labels, f_curves, att_info, pindex)

    return f_curves, att_info, j_curves
end


function _get_attractors_prange(mapper::AttractorsViaRecurrences, ics, n, spp, prange, pidx, featurizer, show_progress)
    progress = ProgressMeter.Progress(n;
        desc="Generating features: ", enabled=show_progress
    )
    # Do the first parameter to build the arrays
    set_parameter!(mapper.integ, pidx, prange[1])
    fs = basins_fractions(mapper, ics; show_progress, N = spp)
    fractions_curves = [fs]
    current_attractors = deepcopy(mapper.bsn_nfo.attractors)
    attractors_info = [current_attractors]

    for (i, p) in enumerate(prange[2:end])
        set_parameter!(mapper.integ, pidx, p)
        reset!(mapper)
        fs = basins_fractions(mapper, ics; show_progress, N = spp)
        push!(fractions_curves, fs)
        push!(attractors_info, deepcopy(mapper.bsn_nfo.attractors))
        ProgressMeter.next!(progress)
    end

    # Set up containers
    example_feature = featurizer(first(values(attractors_info[1])))
    features = typeof(example_feature)[]
    parameter_idxs = Int[]
    # Transform original data into sequential vectors
    for i in eachindex(fractions_curves)
        ai = attractors_info[i]
        A = length(ai)
        append!(parameter_idxs, (i for _ in 1:A))
        for k in keys(ai)
            push!(features, featurizer(ai[k]))
        end
     end
    return features, parameter_idxs, fractions_curves, attractors_info
end


function _label_fractions(clustered_labels, fractions_curves, attractors_info, pindex)
    P = length(fractions_curves)
    original_labels = keytype(first(fractions_curves))[]
    unlabeled_fractions = zeros(P)
    # Transform original data into sequential vectors
    for i in eachindex(fractions_curves)
        fs = fractions_curves[i]
        ai = attractors_info[i]
        unlabeled_fractions[i] = get(fs, -1, 0.0)
        for k in keys(ai)
            push!(original_labels, k)
        end
    end

    # Anyways, time to reconstruct the joint fractions
    joint_fractions = [Dict{Int,Float64}() for _ in 1:P]
    current_p_idx = 0
    for j in eachindex(clustered_labels)
        new_label = clustered_labels[j]
        p_idx = pindex[j]
        if p_idx > current_p_idx
            current_p_idx += 1
            joint_fractions[current_p_idx][-1] = unlabeled_fractions[current_p_idx]
        end
        d = joint_fractions[current_p_idx]
        orig_frac = get(fractions_curves[current_p_idx], original_labels[j], 0)
        d[new_label] = get(d, new_label, 0) + orig_frac
    end

    return joint_fractions
end

