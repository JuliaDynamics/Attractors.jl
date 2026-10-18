abstract type Grid end

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
    grid_steps = SVector{D, Float64}(step.(grid))
    grid_maxima = SVector{D, Float64}(maximum.(grid))
    grid_minima = SVector{D, Float64}(minimum.(grid))
    return RegularGrid(grid_steps, grid_minima, grid_maxima, grid)
end

minmax_grid_extent(g::RegularGrid) = g.grid_minima, g.grid_maxima
mean_cell_diagonal(g::RegularGrid) = norm(g.grid_steps)

struct IrregularGrid{D, T <: Tuple} <: Grid
    grid::T       # cell centers, sorted ascending
end

function IrregularGrid(grid::Tuple)
    all(issorted, grid) || error("IrregularGrid dimensions must be sorted ascending.")
    return IrregularGrid{length(grid), typeof(grid)}(grid)
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
struct SubdivisionBasedGrid{D, R <: AbstractRange} <: Grid
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
as a grid to [`BasinMapRecurrences`](@ref). The input `grid` is an
originally coarse grid (a tuple of `AbstractRange`s).
The state space speed is evaluate in all cells of the `grid`. Cells with small speed
(when compared to the "max" speed) resultin in this cell being subdivided more.
To avoid problems with spikes in the speed, the `q`-th quantile of the velocities
is used as the "max" speed (use `q = 1` for true maximum).
The subdivisions in the resulting grid are clamped to at most value `maxlevel`.

This approach is designed for _continuous time_ systems in which different areas of
the state space flow may have significantly different velocity. In case of
originally coarse grids, this may lead [`BasinMapRecurrences`](@ref)
being stuck in some state space regions with
a small motion speed and false identification of attractors.
"""
function subdivision_based_grid(ds::DynamicalSystem, grid; maxlevel = 4, q = 0.99)
    lvl_array = make_lvl_array(ds, grid, maxlevel, q)
    return SubdivisionBasedGrid(grid, lvl_array)
end

function SubdivisionBasedGrid(grid::NTuple{D, <:AbstractRange}, lvl_array::Array{Int, D}) where {D}
    size(lvl_array) == length.(grid) .- 1 ||
        error("`lvl_array` must have size `length.(grid) .- 1` (one entry per cell).")
    unique_lvls = unique(lvl_array)
    any(<(0), unique_lvls) && error("Level array cannot contain negative values!")

    grid_steps = Dict{Int, Vector{Int}}()
    for i in unique_lvls
        grid_steps[i] = [(length(axis) - 1) * 2^i for axis in grid]   # cells per dim
    end
    grid_maxima = SVector{D, Float64}(maximum.(grid))
    grid_minima = SVector{D, Float64}(minimum.(grid))

    function scale_axis(axis, multiplier)
        ncells = (length(axis) - 1) * 2^multiplier
        return range(first(axis), last(axis), length = ncells + 1)
    end
    multiplier = maximum(keys(grid_steps))
    max_grid = Tuple([scale_axis(axis, multiplier) for axis in grid])
    return SubdivisionBasedGrid(grid_steps, grid_minima, grid_maxima, lvl_array, grid, max_grid)
end

function make_lvl_array(ds::DynamicalSystem, grid, maxlevel, q)
    isdiscretetime(ds) && error("Dynamical system must be continuous time.")
    ncells = length.(grid) .- 1 # remember grid is the cell edges
    indices = CartesianIndices(ncells)
    f, p = dynamic_rule(ds), current_parameters(ds)
    udummy = copy(current_state(ds))
    velocities = zeros(ncells)
    for ind in indices
        u0 = Attractors.generate_ic_on_grid(grid, ind)
        velocity = if !isinplace(ds)
            f(u0, p, 0.0)
        else
            f(udummy, u0, p, 0.0)
            udummy
        end
        speed = norm(velocity)
        velocities[ind] = isnan(speed) ? Inf : speed
    end
    maxvel = quantile(filter(x -> x != Inf, velocities), q)
    ratios = maxvel ./ velocities
    return [round(Int, log2(clamp(x, 1, 2^maxlevel))) for x in ratios]
end


#####################################################################################
# Mapping a state to a cartesian index according to the grid
#####################################################################################
function basin_cell_index(u, g::RegularGrid)
    D = length(g.grid_minima)
    @inbounds for i in 1:D
        if !(g.grid_minima[i] ≤ u[i] < g.grid_maxima[i])
            return CartesianIndex{D}(-1)
        end
    end
    ind = @. floor(Int, (u - g.grid_minima) / g.grid_steps) + 1
    return CartesianIndex{D}(ind...)
end

function basin_cell_index(u, g::IrregularGrid)
    D = length(g.grid)
    for (axis, coord) in zip(g.grid, u)
        if !(first(axis) ≤ coord < last(axis))
            return CartesianIndex{D}(-1)
        end
    end
    inds = map((x, y) -> searchsortedlast(x, y), g.grid, u)
    return CartesianIndex{D}(inds...)
end

function basin_cell_index(u, g::SubdivisionBasedGrid)
    D = length(g.grid)
    s = SVector{D, Float64}(step.(g.grid))
    i0 = basin_cell_index(u, RegularGrid(s, g.grid_minima, g.grid_maxima, g.grid))
    i0 == CartesianIndex{D}(-1) && return i0
    lvl = g.lvl_array[i0]
    max_level = maximum(keys(g.grid_steps))
    ind = @. floor(Int, (u - g.grid_minima) / (s / 2^lvl)) * 2^(max_level - lvl) + 1
    return CartesianIndex{D}(ind...)
end

#####################################################################################
# Pretty printing
#####################################################################################
Base.show(io::IO, g::RegularGrid) = Base.show(io, g.grid)
Base.show(io::IO, g::IrregularGrid) = Base.show(io, g.grid)
function Base.show(io::IO, g::SubdivisionBasedGrid)
    return println(io, "SubdivisionBasedGrid with $(maximum(g.lvl_array)) subdivisions")
end
