export GroupAcrossParameterContinuation
import ProgressMeter
import Mmap

struct GroupAcrossParameterContinuation{A<:AttractorsViaFeaturizing, E} <: BasinsFractionContinuation
    mapper::A
    info_extraction::E
    par_weight::Float64
end

"""
    GroupAcrossParameterContinuation(mapper::AttractorsViaFeaturizing; kwargs...)

A method for [`basins_fractions_continuation`](@ref).
It uses the featurizing approach discussed in [`AttractorsViaFeaturizing`](@ref)
and hence requires an instance of that mapper as an input.
When used in [`basins_fractions_continuation`](@ref), features are extracted
and then grouped across a parameter range. Said differently, all features
of all initial conditions across all parameter values are put into the same "pool"
and then grouped as dictated by the `group_config` of the mapper.
After the grouping is finished the feature label fractions are distributed
to each parameter value they came from.

## Keyword arguments
- `info_extraction::Function` a function that takes as an input a vector of features
  (corresponding to a cluster) and returns a description of the cluster.
  By default, the centroid of the cluster is used.
- `par_weight = 0`: See below the section on MCBB.

## MCBB special version
If the chosen grouping method is [`GroupViaClustering`](@ref), the additional keyword
`par_weight::Real` can be used. If it is ≠ 0, the distance matrix between features
obtains an extra weight that is proportional to the distance `par_weight*|p[i] - p[j]|`
between the parameters used when extracting features.
The range of parameters is normalized to 0-1
such that the largest distance in the parameter space is 1. The normalization is done
because the feature space is also (by default) normalized to 0-1.

This version of the algorithm is the original "MCBB" continuation method described
in [^Gelbrecht2020], besides the improvements of clustering accuracy and performance
done by the developer team of Attractors.jl.

[^Gelbrecht2021]:
    Maximilian Gelbrecht et al 2021, Monte Carlo basin bifurcation analysis,
    [New J. Phys.22 03303](http://dx.doi.org/10.1088/1367-2630/ab7a05)
"""
function GroupAcrossParameterContinuation(
        mapper::AttractorsViaFeaturizing;
        info_extraction = mean_across_features,
        par_weight = 0.0,
    )
    return GroupAcrossParameterContinuation(
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

function basins_fractions_continuation(
        continuation::GroupAcrossParameterContinuation, prange, pidx, ics;
        show_progress = true, samples_per_parameter = 100
    )
    (; mapper, info_extraction, par_weight) = continuation
    spp, n = samples_per_parameter, length(prange)
    features = _get_features_prange(mapper, ics, n, spp, prange, pidx, show_progress)

    # This is a special clause for implementing the MCBB algorithm (weighting
    # also by parameter value, i.e., making the parameter value a feature)
    # It calls a special `group_features` function that also incorporates the
    # parameter value (see below). Otherwise, we call normal `group_features`.
    if mapper.group_config isa GroupViaClustering && par_weight ≠ 0
        labels = group_features(features, mapper.group_config; par_weight, plength = n, spp)

    else
        labels = group_features(features, mapper.group_config)
    end
    fractions_curves, attractors_info =
    label_fractions_across_parameter(labels, n, spp, features[1], info_extraction)
    return fractions_curves, attractors_info
end

function _get_features_prange(mapper::AttractorsViaFeaturizing, ics, n, spp, prange, pidx, show_progress)
    progress = ProgressMeter.Progress(n;
        desc="Generating features", enabled=show_progress, offset = 2,
    )
    # Extract the first possible feature to initialize the features container
    feature = extract_features(mapper, ics; N = 1)
    features = Vector{typeof(feature[1])}(undef, n*spp)
    # Collect features
    for (i, p) in enumerate(prange)
        set_parameter!(mapper.ds, pidx, p)
        current_features = extract_features(mapper, ics; show_progress, N = spp)
        features[((i - 1)*spp + 1):i*spp] .= current_features
        ProgressMeter.next!(progress)
    end
    return features
end

function label_fractions_across_parameter(labels, n, spp, feature, info_extraction)
    # finally we collect/group stuff into their dictionaries
    fractions_curves = Vector{Dict{Int, Float64}}(undef, n)
    dummy_info = info_extraction(feature)
    attractors_info = Vector{Dict{Int, typeof(dummy_info)}}(undef, n)
    for i in 1:n
        # Here we know which indices correspond to which parameter value
        # because they are sequentially increased every `spp`
        # (steps per parameter)
        current_labels = view(labels, ((i - 1)*spp + 1):i*spp)
        current_ids = unique(current_labels)
        # getting fractions is easy; use API function that takes in arrays
        fractions_curves[i] = basins_fractions(current_labels, current_ids)
        attractors_info[i] = Dict(
            id => info_extraction(
                view(current_labels, findall(isequal(id), current_labels))
            ) for id in current_ids
        )
    end
    return fractions_curves, attractors_info
end
