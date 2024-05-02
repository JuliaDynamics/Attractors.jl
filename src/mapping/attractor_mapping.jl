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
    SubdivisionBasedGrid

#########################################################################################
# AttractorMapper structure definition
#########################################################################################
"""
    AttractorMapper(ds::DynamicalSystem, args...; kwargs...) → mapper

Subtypes of `AttractorMapper` are structures that map initial conditions of `ds` to
attractors. Currently available mapping methods:

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
"""
abstract type AttractorMapper end

# Generic pretty printing
function generic_mapper_print(io, mapper)
    ps = 14
    println(io, "$(nameof(typeof(mapper)))")
    println(io, rpad(" system: ", ps), nameof(typeof(mapper.ds)))
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
        ics::Union{StateSpaceSet, Function};
        kwargs...
    )

Approximate the state space fractions `fs` of the basins of attraction of a dynamical
system by mapping initial conditions to attractors using `mapper`
(which contains a reference to a [`DynamicalSystem`](@ref)).
The fractions are simply the ratios of how many initial conditions ended up
at each attractor.

Initial conditions to use are defined by `ics`. It can be:
* a `StateSpaceSet` of initial conditions, in which case all are used.
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
function basins_fractions(mapper::AttractorMapper, ics::Union{AbstractStateSpaceSet, Function};
        show_progress = true, N = 1000, additional_fs::Dict = Dict(),
    )
    used_StateSpaceSet = ics isa AbstractStateSpaceSet
    N = used_StateSpaceSet ? size(ics, 1) : N
    progress = ProgressMeter.Progress(N;
        desc="Mapping initial conditions to attractors:", enabled = show_progress
    )
    fs = Dict{Int, Int}()
    used_StateSpaceSet && (labels = Vector{Int}(undef, N))

    for i ∈ 1:N
        ic = _get_ic(ics, i)
        label = mapper(ic; show_progress)
        fs[label] = get(fs, label, 0) + 1
        used_StateSpaceSet && (labels[i] = label)
        show_progress && ProgressMeter.next!(progress)
    end
    # the non-public-API `additional_fs` i s used in the continuation methods
    additive_dict_merge!(fs, additional_fs)
    N = N + (isempty(additional_fs) ? 0 : sum(values(additional_fs)))
    # Transform count into fraction
    ffs = Dict(k => v/N for (k, v) in fs)
    if used_StateSpaceSet
        return ffs, labels
    else
        return ffs
    end
end

_get_ic(ics::Function, i) = ics()
_get_ic(ics::AbstractStateSpaceSet, i) = ics[i]

"""
    extract_attractors(mapper::AttractorsMapper) → attractors

Return a dictionary mapping label IDs to attractors found by the `mapper`.
This function should be called after calling [`basins_fractions`](@ref)
with the given `mapper` so that the attractors have actually been found first.

For `AttractorsViaFeaturizing`, the attractors are only stored if
the mapper was called with pre-defined initial conditions rather than
a sampler (function returning initial conditions).
"""
extract_attractors(::AttractorMapper) = error("not implemented")


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
# Generic basins of attraction method structure definition
#########################################################################################
# It works for all mappers that define a `basins_fractions` method.
"""
    basins_of_attraction(mapper::AttractorMapper, grid::Tuple) → basins, attractors

Compute the full basins of attraction as identified by the given `mapper`,
which includes a reference to a [`DynamicalSystem`](@ref) and return them
along with (perhaps approximated) found attractors.

`grid` is a tuple of ranges defining the grid of initial conditions that partition
the state space into boxes with size the step size of each range.
For example, `grid = (xg, yg)` where `xg = yg = range(-5, 5; length = 100)`.
The grid has to be the same dimensionality as the state space expected by the
integrator/system used in `mapper`. E.g., a [`ProjectedDynamicalSystem`](@ref)
could be used for lower dimensional projections, etc. A special case here is
a [`PoincareMap`](@ref) with `plane` being `Tuple{Int, <: Real}`. In this special
scenario the grid can be one dimension smaller than the state space, in which case
the partitioning happens directly on the hyperplane the Poincaré map operates on.

`basins_of_attraction` function is a convenience 5-lines-of-code wrapper which uses the
`labels` returned by [`basins_fractions`](@ref) and simply assigns them to a full array
corresponding to the state space partitioning indicated by `grid`.

See also [`convergence_and_basins_of_attraction`](@ref).
"""
function basins_of_attraction(mapper::AttractorMapper, grid::Tuple; kwargs...)
    basins = zeros(Int32, map(length, grid))
    I = CartesianIndices(basins)
    A = StateSpaceSet([generate_ic_on_grid(grid, i) for i in vec(I)])
    fs, labels = basins_fractions(mapper, A; kwargs...)
    attractors = extract_attractors(mapper)
    vec(basins) .= vec(labels)
    return basins, attractors
end

# Type-stable generation of an initial condition given a grid array index
@generated function generate_ic_on_grid(grid::NTuple{B, T}, ind) where {B, T}
    gens = [:(grid[$k][ind[$k]]) for k=1:B]
    quote
        Base.@_inline_meta
        @inbounds return SVector{$B, Float64}($(gens...))
    end
end

#########################################################################################
# Includes
#########################################################################################
include("attractor_mapping_proximity.jl")
include("recurrences/attractor_mapping_recurrences.jl")
include("grouping/attractor_mapping_featurizing.jl")
