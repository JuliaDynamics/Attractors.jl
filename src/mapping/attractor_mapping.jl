# Definition of the attracting mapping API and exporting
# At the end it also includes all files related to mapping

const mx_dimension_sparse = 30

export AttractorMapper,
    AttractorsViaRecurrences,
    AttractorsViaProximity,
    AttractorsViaFeaturizing,
    ClusteringConfig,
    basins_fractions,
    basins_of_attraction,
    automatic_Δt_basins,
    extract_attractors,
    subdivision_based_grid,
    SubdivisionBasedGrid,
    generate_ics_on_grid

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
    text = "$(nameof(typeof(mapper)))"
    println(io, text)
    println(io, rpad(" rule f: ", ps), DynamicalSystemsBase.rulestring(mapper))
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
stystem by mapping initial conditions to attractors using `mapper`
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

If `ics` is a `StateSpaceSet` the function will also return `labels`, which is
_vector_, of equal length to `ics`, that contains the label each initial
condition was mapped to.

See [`AttractorMapper`](@ref) for all possible `mapper` types, and use
[`extract_attractors`](@ref) (after calling `basins_fractions`) to extract
the stored attractors from the `mapper`.

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
extract_attractors(::AttractorMapper) = error("not imlemented")

#########################################################################################
# Generic basins of attraction method structure definition
#########################################################################################
# It works for all mappers that define a `basins_fractions` method.
"""
    basins_of_attraction(mapper::AttractorMapper, grid::Tuple) → basins, attractors
Compute the full basins of attraction as identified by the given `mapper`,
which includes a reference to a [`GeneralizedDynamicalSystem`](@ref) and return them
along with (perhaps approximated) found attractors.

`grid` is a tuple of ranges defining the grid of initial conditions that partition
the state space into boxes with size the step size of each range.
For example, `grid = (xg, yg)` where `xg = yg = range(-5, 5; length = 100)`.
The grid has to be the same dimensionality as the state space expected by the
integrator/system used in `mapper`. E.g., a [`projected_integrator`](@ref)
could be used for lower dimensional projections, etc. A special case here is
a [`poincaremap`](@ref) with `plane` being `Tuple{Int, <: Real}`. In this special
scenario the grid can be one dimension smaller than the state space, in which case
the partitioning happens directly on the hyperplane the Poincaré map operates on.

`basins_of_attraction` function is a convenience 5-lines-of-code wrapper which uses the
`labels` returned by [`basins_fractions`](@ref) and simply assings them to a full array
corresponding to the state space partitioning indicated by `grid`.
"""
function basins_of_attraction(mapper::AttractorMapper, grid::Tuple; kwargs...)
    basins = zeros(Int32, map(length, grid))
    A = generate_ics_on_grid(grid, basins)
    fs, labels = basins_fractions(mapper, A; kwargs...)
    attractors = extract_attractors(mapper)
    vec(basins) .= vec(labels)
    return basins, attractors
end

"""
Generates initial conditions on `grid`, used typically in `basins_of_attraction`. A common
use case is studying basins in small-dimensional systems. In this case,
`generate_ic_on_grid`, using CartesianIndices on all dimensions, works well. A particular
case occurs when studying basins in high-dimensional systems, in which `grid` is
high-dimensional but only has varying values in two dimensions. The first function will
not work, so a special function `_generate_ics_on_grid_vary_few_dims` is used to avoid the
problem. This separation is controlled by:

* `min_num_fixed_dims`: minimum number of fixed dimensions above which it becomes worth it
  to use the special function to generate ics.
"""
function generate_ics_on_grid(grid, basins; min_num_fixed_dims=4)
    grid_lengths = length.(grid)
    number_fixed_dims = length(findall(len->len<=1, grid_lengths))
    if number_fixed_dims >= min_num_fixed_dims
        return _generate_ics_on_grid_vary_few_dims(grid)
    else
        return _generate_ics_on_grid(grid, basins)
    end
end

function _generate_ics_on_grid(grid, basins)
    I = CartesianIndices(basins)
    A = StateSpaceSet([generate_ic_on_grid(grid, i) for i in vec(I)])
    return A
end

# Type-stable generation of an initial condition given a grid array index
@generated function generate_ic_on_grid(grid::NTuple{B, T}, ind) where {B, T}
    gens = [:(grid[$k][ind[$k]]) for k=1:B]
    quote
        Base.@_inline_meta
        @inbounds return SVector{$B, Float64}($(gens...))
    end
end

"""
    Generates ics on grid for a reduced grid, which contains only dimensions that are varying.
    Then construct the full-dimension vector from that. 
"""
function _generate_ics_on_grid_vary_few_dims(grid::NTuple{B, T}) where {B, T} #only worth it if only a few varyign dims
    grid_lengths = length.(grid)
    idxs_varying_dims = findall(len->len>1, grid_lengths)
    reduced_grid = grid[idxs_varying_dims]
    I_reduced = CartesianIndices(length.(reduced_grid))
    A_reduced = [generate_ic_on_grid(reduced_grid, i) for i in vec(I_reduced)]
    
    ics_fixed = [grid_dim[1] for grid_dim in grid]
    A = _expand_A(A_reduced, ics_fixed, idxs_varying_dims)
    return Dataset(A)
end

function _expand_A(vec_reduced, ic_fixed, idxs_varying_dims) where {A}
    vec = [deepcopy(ic_fixed) for _ in vec_reduced]
    for i in eachindex(vec)
        vec[i][idxs_varying_dims] .= vec_reduced[i]
    end
    return vec 
end

#########################################################################################
# Includes
#########################################################################################
include("attractor_mapping_proximity.jl")
include("attractor_mapping_recurrences.jl")
include("attractor_mapping_featurizing.jl")
