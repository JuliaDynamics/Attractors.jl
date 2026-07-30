export stability_quantifiers_along_continuation

# this fuction is called in `global_continuation` when the bmap is the accumulator
function match_and_generate_output(bmap::StabilityQuantifiersAccumulator, quantifiers_cont, rmaps)
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
    return Vector{Dict{Int64, Float64}}(transposed["basin_fraction"]), transposed
end


"""
    stability_quantifiers_along_continuation(
        ds::DynamicalSystem, attractors_cont, pcurve, sampler;
        kw...
    )

Perform a global continuation of all stability quantifiers estimated by
[`StabilityQuantifiersAccumulator`](@ref) using the found attractors of
a previous call to [`global_continuation`](@ref) using the `ds`.
The inputs `pcurve, sampler` are the same as in [`global_continuation`](@ref)
(i.e., a vector and an `InitialConditionSampler`).

This method is special because it always creates a [`BasinMapProximity`](@ref)
for the attractors at a given point along the global continuation,
and then estimates the stability quantifiers using [`StabilityQuantifiersAccumulator`](@ref)
and the proximity map.

There are two reasons to use this method:

1. You are interested in quantifiers related to the convergence time, which is defined
   more rirogously and is estimated more accurately for a proximity map.
2. You want more control over the values of `ε, finite_time, weighting_distribution`,
   all of which are allowed to be `Vector`s with the same length as `pcurve`.
   (they can always be functions)

## Keyword arguments

- `ε = nothing`: given to [`BasinMapProximity`](@ref).
- `proximity_mapper_options = NamedTuple()`: extra keywords for `BasinMapProximity`.
- `distance, finite_time, weighting_distribution`: given to [`StabilityQuantifiersAccumulator`](@ref).
- `basin_entropy = true`: given to [`StabilityQuantifiersAccumulator`](@ref) to enable/disable
    sampled-basin entropy quantifiers.
- `n_basin_entropy = 100`: given to [`StabilityQuantifiersAccumulator`](@ref) as the number
    of nearest neighbors used by sampled-basin entropy quantifiers.

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
        sampler;
        ε = nothing,
        weighting_distribution = EverywhereUniform(),
        finite_time = 1.0,
        basin_entropy = true,
        n_basin_entropy = 100,
        distance = Centroid(),
        proximity_mapper_options = NamedTuple(),
        show_progress = true
    )
    progress = ProgressMeter.Progress(
        length(pcurve); desc = "Continuing accumulator quantifiers:", enabled = show_progress
    )

    # Deprecation: remove this at later version.
    if !(sampler isa InitialConditionSampler)
        @warn "You must now pass a subtype of `InitialConditionSampler` explicitly."
        if sampler isa AbstractVector
            sampler = PrescribedICs(sampler)
        else
            sampler = RandomICSampler(sampler, 1000)
        end
    end

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
            distance = d,
            basin_entropy = basin_entropy,
            n_basin_entropy = n_basin_entropy,
        )

        pics = generate_ics(sampler, p)
        basins_fractions(accumulator, pics; N = length(sampler), show_progress = false)
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
