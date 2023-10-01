#####################################################################################
# Grid construction
#####################################################################################
Base.@kwdef struct RegularGrid{D, R <: AbstractRange} <: Grid
    grid_steps::SVector{D, Float64}
    grid_minima::SVector{D, Float64}
    grid_maxima::SVector{D, Float64}
    grid::NTuple{D, R}
end

minmax_grid_extent(g::RegularGrid) = g.grid_minima, g.grid_maxima
mean_cell_diagonal(g::RegularGrid{D}) where {D} = norm(g.grid_steps)

struct IrregularGrid{D} <: Grid
    grid::NTuple{D,Vector{Float64}}
end
minmax_grid_extent(g::IrregularGrid) = minmax_grid_extent(g.grid)

minmax_grid_extent(g::NTuple) = map(minimum, g), map(maximum, g)
function mean_cell_diagonal(g::NTuple)
    steps = map(r -> mean(diff(r)), g)
    return norm(steps)
end

"""
    SubdivisionBasedGrid(grid::NTuple{D, <:AbstractRange}, lvl_array::Array{Int, D})

Given a coarse `grid` tesselating the state space, construct a `SubdivisionBasedGrid`
based on the given level array `lvl_array`.
The level array is has non-negative integer values.
Value 0 means that the corresponding cell of the coarse `grid` is not subdivided
any further. Value `n > 0` means that the corresponding cell will be subdivided
in total `2^n` times (along each dimension), resulting in finer cells
within the original coarse cell.
"""
struct SubdivisionBasedGrid{D, R <: AbstractRange} <:Grid
    grid_steps::Dict{Int, Vector{Int}}
    grid_minima::SVector{D, Float64}
    grid_maxima::SVector{D, Float64}
    lvl_array::Array{Int, D}
    grid::NTuple{D, R}
    max_grid::NTuple{D, R}
end
minmax_grid_extent(g::SubdivisionBasedGrid) = g.grid_minima, g.grid_maxima
mean_cell_diagonal(g::SubdivisionBasedGrid) = mean_cell_diagonal(g.grid)

"""
    subdivision_based_grid(ds::DynamicalSystem, grid; maxlevel = 4)

Construct a grid structure `SubdivisionBasedGrid` that can be directly passed
as a grid to [`AttractorsViaRecurrences`](@ref). The input `grid` is an
orginally coarse grid (a tuple of `AbstractRange`s).

This approach is designed for _continuous time_ systems in which different areas of
the state space flow may have significantly different velocity. In case of
originally coarse grids, this may lead [`AttractorsViaRecurrences`](@ref)
being stuck in some state space regions with
a small motion speed and false identification of attractors. To prevent this from
happening we provide an algorithm expansion to dynamically evaluate different regions speed
of motion to handle areas of the grid which should be more coarse or dense than others.

To achieve this, function make use of `make_irregular_array` which automatically constructs
an array of discretization levels indices for a `grid`
(a tuple of `AbstractRange`s) originally specified by user.
Upon construction function automatically stores necessary parameters
to further adopt mapping of initial conditions to specific grid density levels.
"""
function subdivision_based_grid(ds::DynamicalSystem, grid; maxlevel = 4)
    lvl_array = make_irregular_array(ds, grid, maxlevel)
    return SubdivisionBasedGrid(grid, lvl_array)
end

function SubdivisionBasedGrid(grid::NTuple{D, <:AbstractRange}, lvl_array::Array{Int, D}) where {D}
    unique_lvls = unique(lvl_array)
    any(<(0), unique_lvls) && error("Level array cannot contain negative values!")
    grid_steps = Dict{Int, Vector{Int}}()
    for i in unique_lvls
        grid_steps[i] = [length(axis)*2^i for axis in grid]
    end
    grid_maxima = SVector{D, Float64}(maximum.(grid))
    grid_minima = SVector{D, Float64}(minimum.(grid))

    function scale_axis(axis, multiplier)
        new_length = length(axis) * (2^multiplier)
        return range(first(axis), last(axis), length=new_length)
    end
    multiplier = maximum(keys(grid_steps))
    scaled_axis = [scale_axis(axis, multiplier) for axis in grid]
    max_grid = Tuple(scaled_axis)

    return SubdivisionBasedGrid(grid_steps, grid_minima, grid_maxima, lvl_array, grid, max_grid)
end

function make_irregular_array(ds::DynamicalSystem, grid, maxlevel = 4)
    isdiscretetime(ds) && error("Dynamical system must be continuous time.")
    indices = CartesianIndices(length.(grid))
    f, p = dynamic_rule(ds), current_parameters(ds)
    udummy = copy(current_state(ds))
    velocities = zeros(length.(grid))
    for ind in indices
        u0 = Attractors.generate_ic_on_grid(grid, ind)
        velocity = if !isinplace(ds)
            f(u0, p, 0.0)
        else
            f(udummy, u0, p, 0.0)
            udummy
        end
        if (isequal(norm(velocity),NaN))
            velocities[ind] = Inf
        else
            velocities[ind] = norm(velocity)
        end
    end

    maxvel = maximum(filter(x -> x != Inf, velocities))
    velratios = maxvel./velocities
    result = [round(Int,log2(clamp(x, 1, 2^maxlevel))) for x in velratios]
    return result
end


#####################################################################################
# Mapping a state to a cartesian index according to the grid
#####################################################################################
function basin_cell_index(u, grid_nfo::RegularGrid)
    D = length(grid_nfo.grid_minima) # compile-type deduction
    @inbounds for i in eachindex(grid_nfo.grid_minima)
        if !(grid_nfo.grid_minima[i] ≤ u[i] ≤ grid_nfo.grid_maxima[i])
            return CartesianIndex{D}(-1)
        end
    end
    # Snap point to grid
    ind = @. round(Int, (u - grid_nfo.grid_minima)/grid_nfo.grid_steps) + 1
    return CartesianIndex{D}(ind...)
end

function basin_cell_index(u, grid_nfo::IrregularGrid)
    D = length(grid_nfo.grid) # compile-type deduction
    for (axis, coord) in zip(grid_nfo.grid, u)
        if coord < first(axis) || coord > last(axis)
            return CartesianIndex{D}(-1)
        end
    end
    cell_indices = map((x, y) -> searchsortedlast(x, y), grid_nfo.grid, u)
    return CartesianIndex{D}(cell_indices...)
end

function basin_cell_index(u, grid_nfo::SubdivisionBasedGrid)
    D = length(grid_nfo.grid) # compile-type deduction
    initial_index = basin_cell_index(u, RegularGrid(SVector{D,Float64}(step.(grid_nfo.grid)), grid_nfo.grid_minima, grid_nfo.grid_maxima, grid_nfo.grid))
    if initial_index == CartesianIndex{D}(-1)
        return initial_index
    end
    cell_area = grid_nfo.lvl_array[initial_index]
    grid_maxima = grid_nfo.grid_maxima
    grid_minima = grid_nfo.grid_minima
    grid_steps = grid_nfo.grid_steps
    max_level = maximum(keys(grid_steps))
    grid_step = (grid_maxima - grid_minima .+ 1) ./ grid_steps[cell_area]
    ind = @. round(Int, (u - grid_minima)/grid_step, RoundDown) * (2^(max_level-cell_area)) + 1
    return CartesianIndex{D}(ind...)
end
