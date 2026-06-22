export stability_measures_along_continuation

function global_continuation(
        ascm::AttractorSeedContinueMatch{<:StabilityMeasuresAccumulator}, pcurve, ics;
        samples_per_parameter = 100, show_progress = true,
        featurizer = nothing, group_config = nothing,
    )
    N = samples_per_parameter
    progress = ProgressMeter.Progress(
        length(pcurve);
        desc = "Global continuation (accumulator):", PMKWARGS..., enabled = show_progress
    )
    mapper = ascm.mapper
    prev_attractors = empty(extract_attractors(mapper))
    additional_ics = typeof(current_state(referenced_dynamical_system(mapper)))[]
    attractors_cont = Dict[]
    # a dict of stability measures (not basin fractions) is stored per parameter
    measures_cont = []
    centroids_cont = Dict[]
    for (i, p) in enumerate(pcurve)
        set_parameters!(referenced_dynamical_system(mapper), p)
        reset_mapper!(mapper)
        empty!(additional_ics)
        for att in values(prev_attractors)
            for u0 in ascm.seeding(att)
                push!(additional_ics, u0)
            end
        end
        if ics isa PerParameterInitialConditions
            pics = ics.generator(p, N)
        else
            pics = ics
        end
        # the return is unused; the accumulator records every mapped initial condition itself
        basins_fractions(mapper, pics; N, additional_ics, show_progress, offset = 2)
        # finalize with optional aggregation. When grouping, group once and reuse the same
        # labels for both the merged measures and the centroids that match groups across
        # parameters.
        if !isnothing(featurizer) && !isnothing(group_config)
            group_map, centroids = _group_and_centroids(
                extract_attractors(mapper), featurizer, group_config
            )
            measures, agg_attractors = finalize_accumulator(mapper; group_map)
            push!(centroids_cont, centroids)
        else
            measures, agg_attractors = finalize_accumulator(mapper)
        end
        prev_attractors = deepcopy(extract_attractors(mapper))
        push!(attractors_cont, agg_attractors)
        push!(measures_cont, measures)
        showvalues = i < length(pcurve) ? [("pcurve index", i + 1)] : []
        ProgressMeter.next!(progress; showvalues)
    end

    # match by feature centroid when grouping, else use the user's matcher
    if !isnothing(featurizer) && !isnothing(group_config)
        rmaps = match_sequentially!(centroids_cont, MatchByFeatureDistance())
        match_sequentially!(attractors_cont, rmaps)
        transposed = accumulator_continuation_output(measures_cont, rmaps)
    else
        rmaps = match_sequentially!(
            attractors_cont, ascm.matcher; pcurve, ds = referenced_dynamical_system(mapper)
        )
        transposed = accumulator_continuation_output(measures_cont, rmaps)
    end
    return transposed, attractors_cont
end

function accumulator_continuation_output(measures_cont, rmaps)
    # match
    for (i, rmap) in enumerate(rmaps)
        for dict in values(measures_cont[i + 1])
            swap_dict_keys!(dict, rmap)
        end
    end
    # "transpose" (i.e., swap nesting order)
    transposed = Dict{String, Vector{Dict{Int64, Any}}}()
    for measure in measures_cont[1]
        measure_name = measure[1]
        transposed[measure_name] = Vector{Dict{Int64, Any}}()
    end
    for measures in measures_cont
        for (measure_name, measure_dict) in measures
            push!(transposed[measure_name], measure_dict)
        end
    end
    return transposed
end


# make sure to allow the possiblity that the proximity options can also be
# vectors of same length as `pcurve`; Same for the distributions
"""
    stability_measures_along_continuation(
        ds::DynamicalSystem, attractors_cont, pcurve, ics;
        kw...
    )

Perform a global continuation of all stability measures estimated by
[`StabilityMeasuresAccumulator`](@ref) using the found attractors of
a previous call to [`global_continuation`](@ref) using the `ds`.

This method is special because it always creates an [`AttractorsViaProximity`](@ref)
mapper for the attractors at a given point along the global continuation,
and then estimates the stability measures using [`StabilityMeasuresAccumulator`](@ref)
and the proximity mapper.

There are two reasons to use this method:

1. You are interested in measures related to the convergence time, which is defined
   more rirogously and is estimated more accurately for a proximity mapper.
2. You want more control over the values of `ε, finite_time, weighting_distribution`,
   all of which are allowed to be `Vector`s with the same length as `pcurve`.
   (they can always be functions)

## Keyword arguments

- `ε = nothing`: given to [`AttractorsViaProximity`](@ref).
- `proximity_mapper_options = NamedTuple()`: extra keywords for `AttractorsViaProximity`.
- `distance, finite_time, weighting_distribution`: given to [`StabilityMeasuresAccumulator`](@ref).
- `samples_per_parameter = 1000`: how many samples to use when estimating stability measures
  via [`StabilityMeasuresAccumulator`](@ref). Ignored when `ics` is not a function.
- `featurizer = nothing`, `group_config = nothing`: passed to [`finalize_accumulator`](@ref).
  When both are provided, attractors with similar features are merged before computing
  stability measures at each parameter step. See [`aggregate_attractor_fractions`](@ref)
  for the same interface.
  Additionally, when both are given, [`MatchByFeatureDistance`](@ref) is automatically
  applied after the loop to ensure that the grouped attractor labels are consistent
  across the parameter curve: groups at consecutive parameter steps are matched by
  minimising the Euclidean distance between their `featurizer` outputs.
"""
function stability_measures_along_continuation(
        ds::DynamicalSystem,
        attractors_cont,
        pcurve,
        ics;
        ε = nothing,
        weighting_distribution = EverywhereUniform(),
        finite_time = 1.0,
        samples_per_parameter = 1000,
        distance = Centroid(),
        proximity_mapper_options = NamedTuple(),
        featurizer = nothing,
        group_config = nothing,
        show_progress = true
    )
    progress = ProgressMeter.Progress(
        length(pcurve); desc = "Continuing accumulator quantifiers:", enabled = show_progress
    )
    N = samples_per_parameter
    measures_cont = []
    centroids_cont = Dict[]
    measure_names = nothing
    for (i, p) in enumerate(pcurve)
        set_parameters!(ds, p)
        attractors = attractors_cont[i]

        if isempty(attractors)
            push!(measures_cont, Dict{String, Dict{Int64, Float64}}())
            push!(centroids_cont, Dict{Int, Any}())
            ProgressMeter.next!(progress)
            continue
        end

        if ε isa AbstractVector
            ε_ = ε[i]
        elseif ε isa Function
            ε_ = ε(p, attractors)
        else
            ε_ = ε
        end
        if weighting_distribution isa AbstractVector
            wd = weighting_distribution[i]
        elseif weighting_distribution isa Function
            wd = weighting_distribution(p, attractors)
        else
            wd = weighting_distribution
        end
        if finite_time isa AbstractVector
            ft = finite_time[i]
        elseif finite_time isa Function
            ft = finite_time(p, attractors)
        else
            ft = finite_time
        end
        if distance isa AbstractVector
            d = distance[i]
        elseif distance isa Function
            d = distance(p, attractors)
        else
            d = distance
        end

        accumulator = StabilityMeasuresAccumulator(
            AttractorsViaProximity(ds, attractors; ε = ε_, proximity_mapper_options...);
            weighting_distribution = wd, finite_time = ft,
            distance = d
        )

        if ics isa PerParameterInitialConditions
            pics = ics.generator(p, N)
        else
            pics = ics
        end
        basins_fractions(accumulator, pics; N, show_progress = false)
        # When grouping, group once and reuse the same labels for both the merged measures
        # and the centroids that match groups across parameters.
        if !isnothing(featurizer) && !isnothing(group_config)
            group_map, centroids = _group_and_centroids(attractors, featurizer, group_config)
            measures, _ = finalize_accumulator(accumulator; group_map)
            push!(centroids_cont, centroids)
        else
            measures, _ = finalize_accumulator(accumulator)
        end
        if measure_names === nothing
            measure_names = collect(keys(measures))
        end
        push!(measures_cont, measures)
        ProgressMeter.next!(progress)
    end

    if !isnothing(featurizer) && !isnothing(group_config) && length(centroids_cont) > 1
        rmaps = match_sequentially!(centroids_cont, MatchByFeatureDistance())
        for (i, rmap) in enumerate(rmaps)
            for dict in values(measures_cont[i + 1])
                swap_dict_keys!(dict, rmap)
            end
        end
    end

    # change the measures format to the expected output
    transposed = Dict{String, Vector{Dict{Int64, Float64}}}()
    if measure_names === nothing
        return transposed
    end
    for measure_name in measure_names
        transposed[measure_name] = Vector{Dict{Int64, Float64}}()
    end
    for measures in measures_cont
        for measure_name in measure_names
            push!(transposed[measure_name], get(measures, measure_name, Dict{Int64, Float64}()))
        end
    end
    return transposed
end
