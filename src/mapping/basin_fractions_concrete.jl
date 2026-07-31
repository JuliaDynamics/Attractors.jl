"""
    basins_fractions(bmap::BasinMap, ics::InitialConditionsSampler; kwargs...)

Approximate the state space fractions `fs` of the basins of attraction of a dynamical
system by mapping initial conditions to attractors using `bmap`
(which contains a reference to a [`DynamicalSystem`](@ref)).
The fractions are simply the ratios of how many initial conditions ended up
at each attractor. Return the fractions as a dictionary mapping
basin IDs (integers) to their fractions.
If you also want to obtain a vector of labels corresponding to each initial condition,
use the function [`basins_fractions_labels`](@ref) instead.

Initial conditions are sampled using a concrete subtype of [`InitialConditionsSampler`](@ref)
and typically will use [`RandomICsSampler`](@ref) or [`PrescribedICs`](@ref).

## Return

The function will always return `fractions`, which is
a dictionary whose keys are the labels given to each attractor
(always integers enumerating the different attractors), and whose
values are the respective basins fractions. The label `-1` is given to any initial condition
where `bmap` could not match to an attractor (this depends on the `bmap` type).

If `ics` is a `StateSpaceSet` the function will also return `labels`, which is a
_vector_, of equal length to `ics`, that contains the label each initial
condition was mapped to.

See [`BasinMap`](@ref) for all possible `bmap` types, and use
[`extract_attractors`](@ref) (after calling `basins_fractions`) to extract
the stored attractors from the `bmap`.
See also [`convergence_and_basins_fractions`](@ref).

## Keyword arguments

* `show_progress = true`: Display a progress bar of the process.
"""
function basins_fractions(bmap::BasinMap, sampler::InitialConditionsSampler; kw...)
    fs, labels = basins_fractions_labels(bmap, sampler; fill_labels = false, kw...)
    return fs
end

"""
    basins_fractions_labels(bmap::BasinMap, sampler::InitialConditionsSampler; kw...)

Same as [`basins_fractions`](@ref) but return two outputs: the fractions dictionary
and a vector of integers which contains the labels of the provided initial conditions.
"""
function basins_fractions_labels(bmap::BasinMap, sampler::InitialConditionsSampler;
        show_progress = true,
        # This is an internal keyword used by `basin_fractions`
        fill_labels = true,
        # this is an internal keyword used in the ASCM global conitnuation
        additional_ics = [],
        # and this is another internal keyword for the offset of the progress bar
        # for when this function is called in global continuation
        offset = 0,
        # and yet another internal keyword used in `global_continuation`
        # to notify which parameters to use in the sampler
        params = current_parameters(referenced_dynamical_system(bmap))
    )
    progress = ProgressMeter.Progress(
        length(sampler);
        desc = "Running basin map:", PMKWARGS..., offset, enabled = show_progress
    )
    labels = Vector{Int}(undef, fill_labels ? length(sampler) : 0)
    ffs = if can_map_individual_ic(bmap)
        basins_fractions_individual(bmap, sampler, progress, labels, additional_ics)
    else
        # collect all initial conditions
        ics = generate_ics(sampler, params)
        icscol = collect(ics) # thankfully this also copies for `Vector`
        append!(icscol, additional_ics)
        basins_fractions_grouped(bmap, icscol, progress, labels)
    end
    return ffs, labels
end

function basins_fractions_individual(bmap, sampler, progress, labels, additional_ics)
    ics = generate_ics(sampler, current_parameters(referenced_dynamical_system(bmap)))
    fs = Dict{Int, Int}()
    for u0 in additional_ics
        label = bmap(u0)
        fs[label] = get(fs, label, 0) + 1
    end
    for (i, ic) in enumerate(ics)
        label = bmap(ic)
        fs[label] = get(fs, label, 0) + 1
        !isempty(labels) && (labels[i] = label)
        ProgressMeter.next!(progress)
    end
    ffs = Dict(k => v / (length(sampler) + length(additional_ics)) for (k, v) in fs)
    return ffs
end

"""
    basins_fractions_grouped(bmap, ics, progress, labels)

Internal function called by `basins_fractions` for `BasinMap`s that
cannot map individual initial conditions to attractor labels.
Must be extended for new mappers that fall under this category.
`ics` is the initial conditions already collected into vector,
`progress` is an initialized progress bar of appropriate size,
`labels` is an initialized container of labels that should be overwritten
in the function call if it is not `empty`.
See the implementation of `BasinMapFeaturizeGroup` for an example.
"""
function basins_fractions_grouped(bmap, ics, progress, labels)
    error("Must be implemented for bmap of type $(nameof(typeof(bmap)))")
end

# Generator dispatch is needed because that's what `Sampler` types return
_get_ic(ics::Base.Generator, i) = first(iterate(ics))
_get_ic(ics::Function, i) = ics()
_get_ic(ics::AbstractVector, i) = ics[i]

"""
    extract_attractors(bmap::AttractorsMapper) → attractors

Return a dictionary mapping label IDs to attractors found by the `bmap`.
This function should be called after mapping initial conditions with `bmap`
(e.g., calling [`basins_fractions`](@ref))
so that the attractors have actually been found first.

For developing a new bmap: extend the internal function `_extract_attractors`.
"""
function extract_attractors(bmap::BasinMap)
    attractors = _extract_attractors(bmap)
    ds = referenced_dynamical_system(bmap)
    # name attractor variables if possible
    # add prior compat entry to fix some incompat version bugs
    if isdefined(DynamicalSystemsBase, :referenced_sciml_model)
        check = isnothing(DynamicalSystemsBase.referenced_sciml_model(ds))
    else
        check = isnothing(DynamicalSystemsBase.referrenced_sciml_model(ds))
    end
    check && return attractors
    names = named_variables(ds)
    for (k, A) in attractors
        attractors[k] = StateSpaceSet(A; names)
    end
    return attractors
end

"internal function for whether the bmap can map individual i.c."
can_map_individual_ic(::BasinMap) = true
