# Definition of the attracting mapping API and exporting
# At the end it also includes all files related to mapping

export BasinMap,
    BasinMapRecurrences,
    BasinMapProximity,
    BasinMapFeaturizeGroup,
    ClusteringConfig,
    basins_fractions,
    basins_fractions_labels,
    convergence_and_basins_of_attraction,
    convergence_and_basins_fractions,
    convergence_time,
    basins_of_attraction,
    automatic_Δt_recurrences,
    extract_attractors,
    subdivision_based_grid,
    SubdivisionBasedGrid,
    reset_mapper!,
    StabilityQuantifiersAccumulator,
    finalize_accumulator

#########################################################################################
# BasinMap structure definition
#########################################################################################
"""
    BasinMap(ds::DynamicalSystem, args...; kwargs...) → bmap

Subtypes of `BasinMap` are structures that map initial conditions of `ds` to unique IDs.
These IDs (typically) correspond to the basins of attraction and corresponding attractors of `ds`.
The found attractors are stored inside the `bmap`, and can be obtained
by calling `attractors = extract_attractors(bmap)` when applicable.

Currently available mapping methods:

* [`BasinMapProximity`](@ref)
* [`BasinMapRecurrences`](@ref)
* [`BasinMapFeaturizeGroup`](@ref)

All `BasinMap` subtypes can be used with [`basins_fractions`](@ref)
or [`basins_of_attraction`](@ref).

In addition, most basin maps can be called as a function of an initial condition:
```julia
label = bmap(u0)
```
and this will on the fly compute and return the label of the attractor `u0` converges at.
The mappers that can do this are:

* [`BasinMapProximity`](@ref)
* [`BasinMapRecurrences`](@ref)
* [`BasinMapFeaturizeGroup`](@ref) with the [`GroupViaHistogram`](@ref)
  or [`GroupViaNearestFeature`](@ref) configurations.

See also [`StabilityQuantifiersAccumulator`](@ref) that extends this interface
to accelerate estimation of stability quantifiers.

## For developers

`BasinMap` defines an extendable interface. A new type needs to subtype
`BasinMap` and implement the following:

- [`extract_attractors`](@ref)
- the internal function `Attractors.referenced_dynamical_system(bmap)`.
- the internal function `Attractors.can_map_individual_ic(bmap)::Bool`.
- and if possible, `id = bmap(u0)`

From these, everything else in the rest of the library just works!

If it is not possible to implement `id = bmap(u0)`, then instead extend
the function `basins_fractions_grouped(bmap, ics, progress, labels)`,
where `ics` is always a `Vector` of initial conditions, `progress` is a preinitialized
progress bar, and `labels` is a preinitialized container of labels.
If `!isempty(labels)`, then its full length must be filled with the ids corresponding
to the first N entries of `ics`.
"""
abstract type BasinMap end

referenced_dynamical_system(bmap::BasinMap) = bmap.ds

# Generic pretty printing
function generic_mapper_print(io, bmap)
    ps = 14
    println(io, "$(nameof(typeof(bmap)))")
    println(io, rpad(" system: ", ps), nameof(typeof(referenced_dynamical_system(bmap))))
    return ps
end
Base.show(io::IO, bmap::BasinMap) = generic_mapper_print(io, bmap)

"""
    convergence_time(bmap::BasinMap) → t

Return the approximate time the `bmap` took to converge to an attractor.
This function should be called just right after `bmap(u0)` was called with
`u0` the initial condition of interest. Hence it is only valid with `BasinMap`
subtypes that support this syntax.

Obtaining the convergence time is computationally free,
so that [`convergence_and_basins_fractions`](@ref) can always
be used instead of [`basins_fractions`](@ref).
"""
function convergence_time end

#########################################################################################
# Includes
#########################################################################################

include("sampler_api.jl")
include("basins_types.jl")
include("basin_fractions_concrete.jl")
include("attractor_mapping_proximity.jl")
include("recurrences/attractor_mapping_recurrences.jl")
include("grouping/attractor_mapping_featurizing.jl")
