using Neighborhood
export ics_from_grid, map_to_basin

# It works for all mappers that define a `basins_fractions` method.
"""
    basins_of_attraction(mapper::AttractorMapper, grid::Tuple) → array_basins_of_attraction

Compute the full basins of attraction as identified by the given `mapper`,
which includes a reference to a [`DynamicalSystem`](@ref) and return them
along with (perhaps approximated) found attractors contained within the 
[`ArrayBasinsOfAttraction`](@ref) structure. 

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

Note that, to ensure backwards compatibility the return type can be decomposed such that 
`basins, attractors = array_basins_of_attraction`.

See also [`convergence_and_basins_of_attraction`](@ref).
"""
function basins_of_attraction(mapper::AttractorMapper, grid::Tuple; kwargs...)
    basins = zeros(Int32, map(length, grid))
    A = ics_from_grid(grid)
    fs, labels = basins_fractions(mapper, A; kwargs...)
    attractors = extract_attractors(mapper)
    vec(basins) .= vec(labels)
    return ArrayBasinsOfAttraction(basins, attractors, grid)
end

"""
    basins_of_attraction(mapper::AttractorMapper, ics; kwargs...) → sampled_basin_of_attraction

Compute the full basins of attraction as identified by the given `mapper`,
which includes a reference to a [`DynamicalSystem`](@ref) and return them
along with (perhaps approximated) found attractors contained within the 
[`SampledBasinsOfAttraction`](@ref) structure.

The initial conditions `ics`, and the keyword arguments `kwargs` are the same
as in [`basins_fractions`](@ref) with the same function signature. This function
is a small convenience wrapper which uses the sampled initial conditions and their 
corresponding labels, from `basins_fractions`, to construct a `SampledBasinsOfAttraction`

Note that, as with the other `basins_of_attraction` function, the return can be decomposed: 
`basins, attractors = array_basins_of_attraction`.

"""
function basins_of_attraction(mapper::AttractorMapper, ics::ValidICS; show_progress = true, N = 1000)
    used_container = ics isa AbstractVector
    ics_vec = used_container ? ics : [ics() for _ in 1:N]
    _, labels = basins_fractions(mapper, ics_vec, show_progress=show_progress)
    attractors = extract_attractors(mapper)
    return SampledBasinsOfAttraction(labels, attractors, ics_vec)
end

"""
    ics_from_grid(grid)

Generate all initial conditions corresponding to the given `grid`
(a state space tessellation used e.g. in [`basins_of_attraction`](@ref)).
"""
function ics_from_grid(grid)
    I = CartesianIndices(length.(grid))
    A = [generate_ic_on_grid(grid, i) for i in vec(I)]
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
    basins_fractions(basins::AbstractArray [,ids]) → fs::Dict
    basins_fractions(BoA::BasinsOfAttraction [,ids]) → fs::Dict

Calculate the state space fraction of the basins of attraction encoded in `basins`.
The elements of `basins` are integers, enumerating the attractor that the entry of
`basins` converges to (i.e., like the output of [`basins_of_attraction`](@ref)).
Return a dictionary that maps attractor IDs to their relative fractions.
Optionally you may give a vector of `ids` to calculate the fractions of only
the chosen ids (by default `ids = unique(basins)`).

The second function signature is a simple wrapper of the first signature 
allowing compatibility with the [`BasinsOfAttraction`](@ref) type

In [Menck2013](@cite) the authors use these fractions to quantify the stability of a basin of
attraction, and specifically how it changes when a parameter is changed.
For this, see [`global_continuation`](@ref).
"""
function basins_fractions(basins::AbstractArray, ids = unique(basins))
    fs = Dict{eltype(basins), Float64}()
    N = length(basins)
    for ξ in ids
        B = count(isequal(ξ), basins)
        fs[ξ] = B/N
    end
    return fs
end

basins_fractions(BoA::BasinsOfAttraction, ids = unique(extract_basins(BoA))) = basins_fractions(extract_basins(BoA), ids)

"""
    convergence_and_basins_of_attraction(mapper::AttractorMapper, grid)

An extension of [`basins_of_attraction`](@ref).
Return `basins, attractors, convergence`, with `basins, attractors` as in
`basins_of_attraction`, and `convergence` being an array with same shape
as `basins`. It contains the time each initial condition took
to converge to its attractor.
It is useful to give to [`shaded_basins_heatmap`](@ref).

See also [`convergence_time`](@ref). Note that this function is not able to
be updated with a [`BasinsOfAttraction`](@ref) return type as the return is
an iterable also contains `iterations`, but `array_basins, iterations` breaks
backwards compatibility.

# Keyword arguments

- `show_progress = true`: show progress bar.
"""
function convergence_and_basins_of_attraction(mapper::AttractorMapper, grid; show_progress = true)
    if length(grid) != dimension(referenced_dynamical_system(mapper))
        @error "The mapper and the grid must have the same dimension"
    end
    basins = zeros(length.(grid))
    ds = referenced_dynamical_system(mapper)
    iterations = zeros(typeof(current_time(ds)), length.(grid))
    I = CartesianIndices(basins)
    progress = ProgressMeter.Progress(
        length(basins); desc = "Basins and convergence: ", dt = 1.0
    )

    for (k, ind) in enumerate(I)
        show_progress && ProgressMeter.update!(progress, k)
        u0 = Attractors.generate_ic_on_grid(grid, ind)
        basins[ind] = mapper(u0)
        iterations[ind] = convergence_time(mapper)
    end
    attractors = extract_attractors(mapper)
    return basins, attractors, iterations
end

"""
    convergence_and_basins_fractions(mapper::AttractorMapper, ics)

An extension of [`basins_fractions`](@ref).
Return `fs, labels, convergence`. The first two are as in `basins_fractions`,
and `convergence` is a vector containing the time each initial condition took
to converge to its attractor.
Only usable with mappers that support `id = mapper(u0)`.

See also [`convergence_time`](@ref).

# Keyword arguments

- `show_progress = true`: show progress bar.
"""
function convergence_and_basins_fractions(mapper::AttractorMapper, ics::ValidICS;
        show_progress = true, N = 1000,
    )
    N = ics isa Function ? N : length(ics)
    progress = ProgressMeter.Progress(N;
        desc="Mapping initial conditions to attractors:", enabled = show_progress
    )
    fs = Dict{Int, Int}()
    labels = Vector{Int}(undef, N)
    iterations = Vector{typeof(current_time(mapper.ds))}(undef, N)

    for i ∈ 1:N
        ic = _get_ic(ics, i)
        label = mapper(ic; show_progress = false)
        fs[label] = get(fs, label, 0) + 1
        labels[i] = label
        iterations[i] = convergence_time(mapper)
        show_progress && ProgressMeter.next!(progress)
    end
    # Transform count into fraction
    ffs = Dict(k => v/N for (k, v) in fs)
    return ffs, labels, iterations
end

"""
    map_to_basin(BoA::BasinOfAttraction, point; kwargs...) → id

Given a `point`, `map_to_basin` interpolates to which basin it should belong.

For `ArrayBasinsOfAttraction` this finds the label of the closest grid cell. 

For `SampledBasinsOfAttraction` this finds the label corresponding to the 
nearest neighbor using `Neighborhood.jl`, in which case the keyword arguments are:
* `tree`: search tree constructor (e.g. `KDTree`, `BallTree`)
* `metric`: distance metric (e.g. `Euclidean()`, `Chebyshev()`)
* `searchstructure_kwargs...`: additional keyword arguments passed to `searchstructure`

For developing a new `BasinOfAttraction` subtype extend the internal function `map_to_domain`

"""
function map_to_basin(BoA::BasinsOfAttraction, point; kwargs...)
    n = map_to_domain(BoA, point; kwargs...)
    return BoA.basins[n]
end

# Map point to index of the nearest member of the domain of the basin of attraction
map_to_domain(BoA::ArrayBasinsOfAttraction, point) = basin_cell_index(point, BoA.grid) # Nearest grid cell
function map_to_domain(BoA::SampledBasinsOfAttraction, point; tree = KDTree, metric = Euclidean(), ss_kwargs...) 
    ss = searchstructure(tree, BoA.sampled_points, metric, ss_kwargs...) # Nearest sampled point
    idx, _ = nn(ss, point, args...; kwargs...)
    return idx
end