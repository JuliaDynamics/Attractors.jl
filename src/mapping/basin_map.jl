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
- `id = bmap(u0)`
- the internal function `Attractors.referenced_dynamical_system(bmap)`.

From these, everything else in the entire rest of the library just works!

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

#########################################################################################
# Includes
#########################################################################################
# Instantiate Grid type so that BasinsOfAttraction subtype ArrayBasinsOfAttraction
# may be loaded which allows attractor_mapping_recurrences.jl to be loaded
abstract type Grid end

include("sampler_api.jl")
include("basin_fractions_concrete.jl")
include("attractor_mapping_proximity.jl")
include("recurrences/attractor_mapping_recurrences.jl")
include("grouping/attractor_mapping_featurizing.jl")

"internal function for whether the bmap can map individual i.c."
allows_mapper_u0(::BasinMap) = true
function allows_mapper_u0(bmap::BasinMapFeaturizeGroup)
    if bmap.group_config isa GroupViaClustering
        return false
    elseif bmap.group_config isa GroupViaPairwiseComparison
        return false
    else
        return true
    end
end
