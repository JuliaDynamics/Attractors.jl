export stability_quantifiers_along_continuation

# This function is practically identical with the original one in
# `continuation_ascm_generic.jl`; only small differences exist
function global_continuation(
        ascm::AttractorSeedContinueMatch{<:StabilityQuantifiersAccumulator}, pcurve, ics;
        samples_per_parameter = 100, show_progress = true,
    )
    N = samples_per_parameter
    progress = ProgressMeter.Progress(
        length(pcurve);
        desc = "Global continuation (accumulator):", PMKWARGS..., enabled = show_progress
    )
    bmap = ascm.bmap
    prev_attractors = empty(extract_attractors(bmap))
    additional_ics = typeof(current_state(referenced_dynamical_system(bmap)))[]
    attractors_cont = Dict[]
    # difference one: this isn't fractions
    quantifiers_cont = []
    for (i, p) in enumerate(pcurve)
        set_parameters!(referenced_dynamical_system(bmap), p)
        reset_mapper!(bmap)
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
        # difference two: we don't care about the return of basins_fractions
        # as initial condition mapping is accumulated anyways
        basins_fractions(bmap, pics; N, additional_ics, show_progress, offset = 2)
        quantifiers = finalize_accumulator(bmap)
        prev_attractors = deepcopy(extract_attractors(bmap))
        push!(attractors_cont, prev_attractors)
        push!(quantifiers_cont, quantifiers)
        showvalues = i < length(pcurve) ? [("pcurve index", i + 1)] : []
        ProgressMeter.next!(progress; showvalues)
    end

    rmaps = match_sequentially!(
        attractors_cont, ascm.matcher; pcurve, ds = referenced_dynamical_system(bmap)
    )
    # and difference four, a bit more involved matching for quantifiers:
    transposed = accumulator_continuation_output(quantifiers_cont, rmaps)
    return transposed, attractors_cont
end

function accumulator_continuation_output(quantifiers_cont, rmaps)
    # match
    for (i, rmap) in enumerate(rmaps)
        for dict in values(quantifiers_cont[i + 1])
            swap_dict_keys!(dict, rmap)
        end
    end
    # "transpose" (i.e., swap nesting order)
    transposed = Dict{String, Vector{Dict{Int64, Any}}}()
    for quantifier in quantifiers_cont[1]
        quantifier_name = quantifier[1]
        transposed[quantifier_name] = Vector{Dict{Int64, Any}}()
    end
    for quantifiers in quantifiers_cont
        for (quantifier_name, quantifier_dict) in quantifiers
            push!(transposed[quantifier_name], quantifier_dict)
        end
    end
    return transposed
end


# make sure to allow the possiblity that the proximity options can also be
# vectors of same length as `pcurve`; Same for the distributions
"""
    stability_quantifiers_along_continuation(
        ds::DynamicalSystem, attractors_cont, pcurve, ics;
        kw...
    )

Perform a global continuation of all stability quantifiers estimated by
[`StabilityQuantifiersAccumulator`](@ref) using the found attractors of
a previous call to [`global_continuation`](@ref) using the `ds`.

This method is special because it always creates an [`BasinMapProximity`](@ref)
bmap for the attractors at a given point along the global continuation,
and then estimates the stability quantifiers using [`StabilityQuantifiersAccumulator`](@ref)
and the proximity bmap.

There are two reasons to use this method:

1. You are interested in quantifiers related to the convergence time, which is defined
   more rirogously and is estimated more accurately for a proximity bmap.
2. You want more control over the values of `ε, finite_time, weighting_distribution`,
   all of which are allowed to be `Vector`s with the same length as `pcurve`.
   (they can always be functions)

## Keyword arguments

- `ε = nothing`: given to [`BasinMapProximity`](@ref).
- `proximity_mapper_options = NamedTuple()`: extra keywords for `BasinMapProximity`.
- `distance, finite_time, weighting_distribution`: given to [`StabilityQuantifiersAccumulator`](@ref).
- `samples_per_parameter = 1000`: how many samples to use when estimating stability quantifiers
  via [`StabilityQuantifiersAccumulator`](@ref). Ignored when `ics` is not a function.

## Aggregating attractors

This function computes stability quantifiers for whatever attractors it is given. To obtain
quantifiers for *aggregated* groups of attractors, first merge them with
[`aggregate_continuation`](@ref) and pass the resulting `agg_attractors_cont` here:
each merged group is then treated as a single attractor, so all quantifiers (including those
that need the raw basin data, like medians and critical shock magnitudes) are computed
correctly for the group, with IDs that stay consistent along the parameter axis.
"""
function stability_quantifiers_along_continuation(
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
        show_progress = true
    )
    progress = ProgressMeter.Progress(
        length(pcurve); desc = "Continuing accumulator quantifiers:", enabled = show_progress
    )
    N = samples_per_parameter
    quantifiers_cont = []
    quantifier_names = nothing
    for (i, p) in enumerate(pcurve)
        set_parameters!(ds, p)
        attractors = attractors_cont[i]

        if isempty(attractors)
            push!(quantifiers_cont, Dict{String, Dict{Int64, Float64}}())
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

        accumulator = StabilityQuantifiersAccumulator(
            BasinMapProximity(ds, attractors; ε = ε_, proximity_mapper_options...);
            weighting_distribution = wd, finite_time = ft,
            distance = d
        )

        if ics isa PerParameterInitialConditions
            pics = ics.generator(p, N)
        else
            pics = ics
        end
        basins_fractions(accumulator, pics; N, show_progress = false)
        quantifiers = finalize_accumulator(accumulator)
        if quantifier_names === nothing
            quantifier_names = collect(keys(quantifiers))
        end
        push!(quantifiers_cont, quantifiers)
        ProgressMeter.next!(progress)
    end

    # change the quantifiers format to the expected output
    transposed = Dict{String, Vector{Dict{Int64, Float64}}}()
    if quantifier_names === nothing
        return transposed
    end
    for quantifier_name in quantifier_names
        transposed[quantifier_name] = Vector{Dict{Int64, Float64}}()
    end
    for quantifiers in quantifiers_cont
        for quantifier_name in quantifier_names
            push!(transposed[quantifier_name], get(quantifiers, quantifier_name, Dict{Int64, Float64}()))
        end
    end
    return transposed
end
