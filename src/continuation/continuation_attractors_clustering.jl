export ClusteringAttractorsContinuation
import ProgressMeter
import Mmap

struct ClusteringAttractorsContinuation{A, G} <: BasinsFractionContinuation
    mapper::A
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
        cluster_in_slice = 100., 
        par_weight = 0.0,
        group_config = GroupViaClustering()
    )
    return ClusteringAttractorsContinuation(
        mapper, cluster_in_slice, par_weight, group_config
    )
end


function basins_fractions_continuation(
        continuation::ClusteringAttractorsContinuation, prange, pidx, ics;
        show_progress = true, samples_per_parameter = 100
    )
    (; mapper, cluster_in_slice, par_weight, group_config) = continuation
    spp, n = samples_per_parameter, length(prange)

    keys,  pindex, vecs, f_curves, att_info = _get_attractors_prange(mapper, ics, n, spp, prange, pidx, show_progress)

    dists = _get_dist_matrix(vecs, pindex, prange, par_weight, cluster_in_slice, group_config)

    cluster_labels = _cluster_across_parameters(dists, vecs, group_config)

    _label_fractions!(cluster_labels, n, keys, pindex, prange, f_curves, att_info)

    return f_curves, att_info
end


function _get_attractors_prange(mapper::AttractorsViaRecurrences, ics, n, spp, prange, pidx, show_progress)
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
    # collect Datasets
    vec_att = Dataset[]
    par_index = Int64[]
    key_array = Int64[]
    for (k,att) in enumerate(attractors_info)
        for a in att
            push!(vec_att, a[2])
            push!(key_array, a[1])
            push!(par_index, k)
        end
    end
    return key_array, par_index, vec_att, fractions_curves, attractors_info
end


function _get_dist_matrix(vec_att, p_index, p_range, par_weight, cluster_in_slice, group_config::GroupViaClustering)
    metric = group_config.clust_distance_metric
    # Construct distance matrix
    dists = pairwise(metric, vec_att; symmetric = true)
    par_weight = 1/norm(p_range[end]-p_range[1])
    for (k,pk) in enumerate(p_index)
        for (j,pj) in enumerate(p_index)
            dists[k,j] += par_weight*norm(p_range[pk]-p_range[pj])
            # Add a weight to the dist if we do not want to cluster for the same parameter slice
            if p_range[pk] == p_range[pj]
                dists[k,j] += cluster_in_slice
            end
        end
    end
    return dists
end

function _cluster_across_parameters(dists, vecs, group_config)
    gc = group_config
    # DO THE OPTIMIZATION HERE
    # cluster with dbscan
    # ϵ_optimal =  _extract_ϵ_optimal(vecs, gc)
    ϵ_optimal =  gc.optimal_radius_method
    cluster_labels = _cluster_distances_into_labels(dists, ϵ_optimal, gc.min_neighbors)
    # dbscanresult = dbscan(dists, 0.3, 1)
    # cluster_labels = cluster_assignment(dbscanresult)
    return cluster_labels
end

function _label_fractions!(cluster_labels, n, key_array, pindex, prange, fractions_curves, attractors_info)
    # And finally collect/group stuff into their dictionaries
    c = 0; next_label = maximum(cluster_labels) + 1
    for i in 1:n
        l = sum(pindex .== i)
        current_labels = view(cluster_labels, c+1:c+l)
        d = Dict()
        for k in c+1:c+l
            kk = key_array[k] # Fetch original key number 
            if cluster_labels[k] == -1
                # if the attractor is in the junk cluster we give it a proper number
                push!(d, kk => next_label) 
                next_label += 1
            else
                push!(d, kk => cluster_labels[k])
            end
        end
        swap_dict_keys!(fractions_curves[i], d)
        swap_dict_keys!(attractors_info[i], d)
        c += l
    end

    return fractions_curves, attractors_info
end

