#####################################################################################
# Grid construction
#####################################################################################
struct RegularGrid{D, R <: AbstractRange} <: Grid
    grid_steps::SVector{D, Float64}
    grid_minima::SVector{D, Float64}
    grid_maxima::SVector{D, Float64}
    grid::NTuple{D, R}
end
function RegularGrid(grid::NTuple)
    D = length(grid)
    grid_steps = SVector{D,Float64}(step.(grid))
    grid_maxima = SVector{D,Float64}(maximum.(grid))
    grid_minima = SVector{D,Float64}(minimum.(grid))
    return RegularGrid(grid_steps, grid_minima, grid_maxima, grid)
end

minmax_grid_extent(g::RegularGrid) = g.grid_minima, g.grid_maxima
mean_cell_diagonal(g::RegularGrid{D}) where {D} = norm(g.grid_steps)

struct IrregularGrid{D} <: Grid
    grid::NTuple{D, Vector{Float64}}
end
minmax_grid_extent(g::IrregularGrid) = minmax_grid_extent(g.grid)

minmax_grid_extent(g::NTuple) = map(minimum, g), map(maximum, g)
function mean_cell_diagonal(g::NTuple)
    steps = map(r -> mean(diff(r)), g)
    return norm(steps)
end
mean_cell_diagonal(g::IrregularGrid) = mean_cell_diagonal(g.grid)

"""
    SubdivisionBasedGrid(grid::NTuple{D, <:AbstractRange}, lvl_array::Array{Int, D})

Given a coarse `grid` tesselating the state space, construct a `SubdivisionBasedGrid`
based on the given level array `lvl_array` that should have the same dimension as `grid`.
The level array has non-negative integer values, with 0 meaning that the
corresponding cell of the coarse `grid` should not be subdivided any further.
Value `n > 0` means that the corresponding cell will be subdivided
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
    subdivision_based_grid(ds::DynamicalSystem, grid; maxlevel = 4, q = 0.99)

Construct a grid structure [`SubdivisionBasedGrid`](@ref) that can be directly passed
as a grid to [`AttractorsViaRecurrences`](@ref). The input `grid` is an
orginally coarse grid (a tuple of `AbstractRange`s).
The state space speed is evaluate in all cells of the `grid`. Cells with small speed
(when compared to the "max" speed) resultin in this cell being subdivided more.
To avoid problems with spikes in the speed, the `q`-th quantile of the velocities
is used as the "max" speed (use `q = 1` for true maximum).
The subdivisions in the resulting grid are clamped to at most value `maxlevel`.

This approach is designed for _continuous time_ systems in which different areas of
the state space flow may have significantly different velocity. In case of
originally coarse grids, this may lead [`AttractorsViaRecurrences`](@ref)
being stuck in some state space regions with
a small motion speed and false identification of attractors.
"""
function subdivision_based_grid(ds::DynamicalSystem, grid; maxlevel = 4, q = 0.99)
    lvl_array = make_lvl_array(ds, grid, maxlevel, q)
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

function make_lvl_array(ds::DynamicalSystem, grid, maxlevel, q)
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
        speed = norm(velocity)
        if isnan(speed)
            velocities[ind] = Inf
        else
            velocities[ind] = speed
        end
    end

    maxvel = quantile(filter(x -> x != Inf, velocities), q)
    # large ratio means small velocity means high subdivision
    ratios = maxvel ./ velocities
    # subdivision is just the log2 of the ratio. We do this fancy
    # computation because this way zeros are handled correctly
    # (and we also clamp the values of the level array correctly in 0-maxlevel)
    result = [round(Int,log2(clamp(x, 1, 2^maxlevel))) for x in ratios]
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
