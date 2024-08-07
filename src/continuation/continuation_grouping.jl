export FeaturizeGroupAcrossParameter
import ProgressMeter
import Mmap

struct FeaturizeGroupAcrossParameter{A<:AttractorsViaFeaturizing, E} <: GlobalContinuationAlgorithm
    mapper::A
    info_extraction::E
    par_weight::Float64
end

"""
    FeaturizeGroupAcrossParameter <: GlobalContinuationAlgorithm
    FeaturizeGroupAcrossParameter(mapper::AttractorsViaFeaturizing; kwargs...)

A method for [`global_continuation`](@ref).
It uses the featurizing approach discussed in [`AttractorsViaFeaturizing`](@ref)
and hence requires an instance of that mapper as an input.
When used in [`global_continuation`](@ref), features are extracted
and then grouped across a parameter range. Said differently, all features
of all initial conditions across all parameter values are put into the same "pool"
and then grouped as dictated by the `group_config` of the mapper.
After the grouping is finished the feature label fractions are distributed
to each parameter value they came from.

This continuation method is based on, but strongly generalizes, the approaches
in the papers [Gelbrecht2020](@cite) and [Stender2021](@cite).

## Keyword arguments

- `info_extraction::Function` a function that takes as an input a vector of feature-vectors
  (corresponding to a cluster) and returns a description of the cluster.
  By default, the centroid of the cluster is used.
  This is what the `attractors_cont` contains in the return of `global_continuation`.
"""
function FeaturizeGroupAcrossParameter(
        mapper::AttractorsViaFeaturizing;
        info_extraction = mean_across_features,
        par_weight = 0.0,
    )
    return FeaturizeGroupAcrossParameter(
        mapper, info_extraction, par_weight
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

function global_continuation(
        continuation::FeaturizeGroupAcrossParameter, pcurve, ics;
        show_progress = true, samples_per_parameter = 100
    )
    (; mapper, info_extraction, par_weight) = continuation
    spp, n = samples_per_parameter, length(pcurve)
    features = _get_features_pcurve(mapper, ics, n, spp, pcurve, show_progress)

    # This is a special clause for implementing the MCBB algorithm (weighting
    # also by parameter value, i.e., making the parameter value a feature)
    # It calls a special `group_features` function that also incorporates the
    # parameter value (see below). Otherwise, we call normal `group_features`.
    # TODO: We have deprecated this special clause. In the next version we need to cleanup
    # the source code and remove the `par_weight` and its special treatment in `group_features`.
    if mapper.group_config isa GroupViaClustering && par_weight ≠ 0
        labels = group_features(features, mapper.group_config; par_weight, plength = n, spp)
    else
        labels = group_features(features, mapper.group_config)
    end
    fractions_cont, attractors_cont =
    label_fractions_across_parameter(labels, 1features, n, spp, info_extraction)
    return fractions_cont, attractors_cont
end

function _get_features_pcurve(mapper::AttractorsViaFeaturizing, ics, n, spp, pcurve, show_progress)
    progress = ProgressMeter.Progress(n;
        desc="Generating features", enabled=show_progress, offset = 2,
    )
    # Extract the first possible feature to initialize the features container
    feature = extract_features(mapper, ics; N = 1)
    features = Vector{typeof(feature[1])}(undef, n*spp)
    # Collect features
    for (i, p) in enumerate(pcurve)
        set_parameters!(mapper.ds, p)
        current_features = extract_features(mapper, ics; show_progress, N = spp)
        features[((i - 1)*spp + 1):i*spp] .= current_features
        ProgressMeter.next!(progress)
    end
    return features
end

function label_fractions_across_parameter(labels, features, n, spp, info_extraction)
    # finally we collect/group stuff into their dictionaries
    fractions_cont = Vector{Dict{Int, Float64}}(undef, n)
    dummy_info = info_extraction([first(features)])
    attractors_cont = Vector{Dict{Int, typeof(dummy_info)}}(undef, n)
    for i in 1:n
        # Here we know which indices correspond to which parameter value
        # because they are sequentially increased every `spp`
        # (steps per parameter)
        current_labels = view(labels, ((i - 1)*spp + 1):i*spp)
        current_features = view(features, ((i - 1)*spp + 1):i*spp)
        current_ids = unique(current_labels)
        # getting fractions is easy; use API function that takes in arrays
        fractions_cont[i] = basins_fractions(current_labels, current_ids)
        attractors_cont[i] = Dict(
            id => info_extraction(
                view(current_features, findall(isequal(id), current_labels))
            ) for id in current_ids
        )
    end
    return fractions_cont, attractors_cont
end
