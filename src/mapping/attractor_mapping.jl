# Definition of the attracting mapping API and exporting
# At the end it also includes all files related to mapping

export AttractorMapper,
    AttractorsViaRecurrences,
    AttractorsViaProximity,
    AttractorsViaFeaturizing,
    ClusteringConfig,
    basins_fractions,
    convergence_and_basins_of_attraction,
    convergence_and_basins_fractions,
    convergence_time,
    basins_of_attraction,
    automatic_Δt_basins,
    extract_attractors,
    subdivision_based_grid,
    SubdivisionBasedGrid,
    reset_mapper!,
    StabilityMeasuresAccumulator,
    finalize_accumulator

ValidICS = Union{AbstractVector, Function}

#########################################################################################
# AttractorMapper structure definition
#########################################################################################
"""
    AttractorMapper(ds::DynamicalSystem, args...; kwargs...) → mapper

Subtypes of `AttractorMapper` are structures that map initial conditions of `ds` to
attractors. The found attractors are stored inside the `mapper`, and can be obtained
by calling `attractors = extract_attractors(mapper)`.

Currently available mapping methods:

* [`AttractorsViaProximity`](@ref)
* [`AttractorsViaRecurrences`](@ref)
* [`AttractorsViaFeaturizing`](@ref)

All `AttractorMapper` subtypes can be used with [`basins_fractions`](@ref)
or [`basins_of_attraction`](@ref).

In addition, some mappers can be called as a function of an initial condition:
```julia
label = mapper(u0)
```
and this will on the fly compute and return the label of the attractor `u0` converges at.
The mappers that can do this are:
* [`AttractorsViaProximity`](@ref)
* [`AttractorsViaRecurrences`](@ref)
* [`AttractorsViaFeaturizing`](@ref) with the [`GroupViaHistogram`](@ref) configuration.

See also [`StabilityMeasuresAccumulator`](@ref) that extends this interface
to accelerate estimation of stability measures.

## For developers

`AttractorMapper` defines an extendable interface. A new type needs to subtype
`AttractorMapper` and implement [`extract_attractors`](@ref), `id = mapper(u0)`
and the internal function `Attractors.referenced_dynamical_system(mapper)`.
From these, everything else in the entire rest of the library just works!
If it is not possible to implement `id = mapper(u0)`, then instead extend
`basins_fractions(mapper, ics)` with `ics` a vector of initial conditions.
"""
abstract type AttractorMapper end

referenced_dynamical_system(mapper::AttractorMapper) = mapper.ds

# Generic pretty printing
function generic_mapper_print(io, mapper)
    ps = 14
    println(io, "$(nameof(typeof(mapper)))")
    println(io, rpad(" system: ", ps), nameof(typeof(referenced_dynamical_system(mapper))))
    return ps
end
Base.show(io::IO, mapper::AttractorMapper) = generic_mapper_print(io, mapper)

#########################################################################################
# Generic basin fractions method structure definition
#########################################################################################
# It works for all mappers that define the function-like-object behavior
"""
    basins_fractions(
        mapper::AttractorMapper,
        ics::Union{AbstractVector, Function};
        kwargs...
    )

Approximate the state space fractions `fs` of the basins of attraction of a dynamical
system by mapping initial conditions to attractors using `mapper`
(which contains a reference to a [`DynamicalSystem`](@ref)).
The fractions are simply the ratios of how many initial conditions ended up
at each attractor.

Initial conditions to use are defined by `ics`. It can be:
* an `AbstractVector` of initial conditions, in which case all are used.
  Typically this is a `StateSpaceSet`.
* a 0-argument function `ics()` that spits out random initial conditions.
  Then `N` random initial conditions are chosen.
  See [`statespace_sampler`](@ref) to generate such functions.

## Return

The function will always return `fractions`, which is
a dictionary whose keys are the labels given to each attractor
(always integers enumerating the different attractors), and whose
values are the respective basins fractions. The label `-1` is given to any initial condition
where `mapper` could not match to an attractor (this depends on the `mapper` type).

If `ics` is a `StateSpaceSet` the function will also return `labels`, which is a
_vector_, of equal length to `ics`, that contains the label each initial
condition was mapped to.

See [`AttractorMapper`](@ref) for all possible `mapper` types, and use
[`extract_attractors`](@ref) (after calling `basins_fractions`) to extract
the stored attractors from the `mapper`.
See also [`convergence_and_basins_fractions`](@ref).

## Keyword arguments

* `N = 1000`: Number of random initial conditions to generate in case `ics` is a function.
* `show_progress = true`: Display a progress bar of the process.
"""
function basins_fractions(mapper::AttractorMapper, ics::ValidICS;
        show_progress = true, N = 1000, additional_fs::Dict = Dict(),
    )
    used_container = ics isa AbstractVector
    N = used_container ? length(ics) : N
    progress = ProgressMeter.Progress(N;
        desc="Mapping initial conditions to attractors:", enabled = show_progress
    )
    fs = Dict{Int, Int}()
    used_container && (labels = Vector{Int}(undef, N))

    for i ∈ 1:N
        ic = _get_ic(ics, i)
        label = mapper(ic; show_progress)
        fs[label] = get(fs, label, 0) + 1
        used_container && (labels[i] = label)
        show_progress && ProgressMeter.next!(progress)
    end
    # the non-public-API `additional_fs` i s used in the continuation methods
    additive_dict_merge!(fs, additional_fs)
    N = N + (isempty(additional_fs) ? 0 : sum(values(additional_fs)))
    # Transform count into fraction
    ffs = Dict(k => v/N for (k, v) in fs)
    if used_container
        return ffs, labels
    else
        return ffs
    end
end

_get_ic(ics::Function, i) = ics()
_get_ic(ics::AbstractVector, i) = ics[i]

"""
    extract_attractors(mapper::AttractorsMapper) → attractors

Return a dictionary mapping label IDs to attractors found by the `mapper`.
This function should be called after calling [`basins_fractions`](@ref)
with the given `mapper` so that the attractors have actually been found first.

For developing a new mapper: extend the internal function `_extract_attractors`.
"""
function extract_attractors(mapper::AttractorMapper)
    attractors = _extract_attractors(mapper)
    ds = referenced_dynamical_system(mapper)
    # name attractor variables if possible
    isnothing(referrenced_sciml_model(ds)) && return attractors
    names = named_variables(ds)
    for (k, A) in attractors
        attractors[k] = StateSpaceSet(A; names)
    end
    return attractors
end


"""
    convergence_time(mapper::AttractorMapper) → t

Return the approximate time the `mapper` took to converge to an attractor.
This function should be called just right after `mapper(u0)` was called with
`u0` the initial condition of interest. Hence it is only valid with `AttractorMapper`
subtypes that support this syntax.

Obtaining the convergence time is computationally free,
so that [`convergence_and_basins_fractions`](@ref) can always
be used instead of [`basins_fractions`](@ref).
"""
function convergence_time end



#########################################################################################
# Includes
#########################################################################################
include("attractor_mapping_proximity.jl")
include("recurrences/attractor_mapping_recurrences.jl")
include("grouping/attractor_mapping_featurizing.jl")
include("stability_measures_accumulator.jl")