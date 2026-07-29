export FeaturizeGroupAcrossParameter
import ProgressMeter
import Mmap

struct FeaturizeGroupAcrossParameter{B <: BasinMapFeaturizeGroup, E} <: GlobalContinuationAlgorithm
    bmap::B
    info_extraction::E
    store_features::Bool
end

"""
    FeaturizeGroupAcrossParameter <: GlobalContinuationAlgorithm
    FeaturizeGroupAcrossParameter(bmap::BasinMapFeaturizeGroup [, info_extraction]; kw...)

A method for [`global_continuation`](@ref) that featurizes and groups
trajectories across the whole parameter axis to establish a continuation of the groups.

`info_extraction::Function` is a function that takes as an input a vector of feature-vectors
(corresponding to a cluster) and returns a description of the cluster.
By default, the centroid of the cluster is used as its description.
This output becomes the `attractors` representation in [`GlobalContinuationOutput`](@ref).

If the keyword `store_features` is `true` (default), then
this method adds some information into the `other` field of `GlobalContinuationOutput`:
- `"features"` contains the features per parameter step
- `"labels"` contains their corresponding labels (both containers are sorted)
which allows subsequent analysis of the grouped features.

## Description

It uses the featurizing approach discussed in [`BasinMapFeaturizeGroup`](@ref)
and hence requires an instance of that bmap as an input.
When used in [`global_continuation`](@ref), features are extracted
and then grouped across a parameter range. Said differently, all features
of all initial conditions across all parameter values are put into the same "pool"
and then grouped as dictated by the `group_config` of the bmap.
After the grouping is finished the feature label fractions are distributed
to each parameter value they came from.

This continuation method is based on, but strongly generalizes, the approaches
in the papers [Gelbrecht2020](@cite) and [Stender2021](@cite).
"""
function FeaturizeGroupAcrossParameter(bmap::BasinMapFeaturizeGroup; store_features = true)
    return FeaturizeGroupAcrossParameter(bmap, mean_across_features, store_features)
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

function global_continuation(
        continuation::FeaturizeGroupAcrossParameter, pcurve::Vector, sampler::InitialConditionSampler;
        show_progress = true,
    )
    (; bmap, info_extraction) = continuation
    # spp means 'samples per parameter'
    spp, n = length(sampler), length(pcurve)
    features = _get_features_pcurve(bmap, sampler, n, spp, pcurve, show_progress)
    labels = group_features(features, bmap.group_config)
    fractions_cont, attractors_cont, feat_cont, label_cont = label_fractions_across_parameter(
        labels, features, n, spp, info_extraction
    )
    # wrap output into the designated type:
    if continuation.store_features
        other = Dict("features" => feat_cont, "labels" => label_cont)
    else
        other = Dict{String, Nothing}()
    end
    out = GlobalContinuationOutput(attractors_cont, fractions_cont, Dict{String,Vector}(), other, pcurve)
    return out
end

function _get_features_pcurve(bmap::BasinMapFeaturizeGroup, sampler, n, spp, pcurve, show_progress)
    progress = ProgressMeter.Progress(
        n; desc = "Generating features", enabled = show_progress, offset = 2,
    )
    # Extract the first possible feature to initialize the features container
    u0s = generate_ics(sampler)
    set_parameters!(referenced_dynamical_system(bmap), first(pcurve))
    current_features = extract_features(bmap, u0s; show_progress, N = length(sampler))
    features = Vector{typeof(current_features[1])}(undef, n * spp)
    features[1:spp] .= current_features
    # Collect features across parameter axis
    for i in 2:length(pcurve)
        set_parameters!(referenced_dynamical_system(bmap), pcurve[i])
        u0s = generate_ics(sampler)
        current_features = extract_features(bmap, u0s; show_progress, N = length(sampler))
        features[((i - 1) * spp + 1):(i * spp)] .= current_features
        ProgressMeter.next!(progress)
    end
    return features
end

function label_fractions_across_parameter(labels, features, n, spp, info_extraction)
    # finally we collect/group stuff into their dictionaries
    fractions_cont = Vector{Dict{Int, Float64}}(undef, n)
    dummy_info = info_extraction([first(features)])
    attractors_cont = Vector{Dict{Int, typeof(dummy_info)}}(undef, n)
    features_cont = []
    labels_cont = []
    for i in 1:n
        # Here we know which indices correspond to which parameter value
        # because they are sequentially increased every `spp`
        # (steps per parameter)
        current_labels = view(labels, ((i - 1) * spp + 1):(i * spp))
        current_features = view(features, ((i - 1) * spp + 1):(i * spp))
        push!(labels_cont, current_labels)
        push!(features_cont, current_features)
        current_ids = unique(current_labels)
        # getting fractions is easy; use API function that takes in arrays
        fractions_cont[i] = basins_fractions(current_labels, current_ids)
        attractors_cont[i] = Dict(
            id => info_extraction(
                    view(current_features, findall(isequal(id), current_labels))
                ) for id in current_ids
        )
    end
    return fractions_cont, attractors_cont, features_cont, labels_cont
end
