using Attractors, CairoMakie
using PredefinedDynamicalSystems

# Set up dynamical system: (magnetic pendulum)
ds = PredefinedDynamicalSystems.magnetic_pendulum(d=0.2, α=0.2, ω=0.8, N=3)
psys = ProjectedDynamicalSystem(ds, [1, 2], [0.0, 0.0])
xg = yg = range(-4, 4; length = 21)
grid = (xg, yg)
mapper = AttractorsViaRecurrences(psys, grid; Δt = 0.2)
grid_nfo = mapper.bsn_nfo.grid_nfo

# %%

using Attractors: RegularGrid, IrregularGrid, SubdivisionBasedGrid

fig = Figure()
ax = Axis(fig[1,1])

# function initialize_grid!(ax, grid::RegularGrid; kwargs...)
#     # Plot grid axes
#     for i in 1:2
#         g = grid.grid[i]
#         f! = i == 1 ? vlines! : hlines!
#         f!(ax, g; linewidth = 1, color = :black, kwargs...)
#         lims! = i == 1 ? xlims! : ylims!
#         lims!(ax, grid.grid_minima[i], grid.grid_maxima[i])
#     end
#     # Populate grid with observable rectangles given the cartesian indices
#     # of the grid. These rectanges can be then accessed later and alter
#     # their color!
#     idxs = CartesianIndices(length.(grid.grid))

# end

# Populate the grid with poly! rectangle plots. However! The rectangles
# correspond to the same "cells" of the grid. Additionally, all
# rectangles are colored with an _observable_, that can be accessed
# later using the `basin_cell_index` function. The observable
# holds the face color of the rectangle!
function initialize_cells2!(ax, grid; kwargs...)
    # These are all possible outputs of the `basin_cell_index` function
    idxs = all_cartesian_idxs(grid)
    color_obs = Matrix{Any}(undef, size(idxs)...)
    # We now need to reverse-engineer
    for i in idxs
        rect = cell_index_to_rect(i, grid)
        color = Observable(rand(Makie.RGBA))
        color_obs[i] = color
        poly!(ax, rect; color = color, strokecolor = :black, strokewidth = 0.5)
    end
    # Set the axis limits better
    mini, maxi = Attractors.minmax_grid_extent(grid)
    xlims!(ax, mini[1], maxi[1])
    ylims!(ax, mini[2], maxi[2])
    return color_obs
end
all_cartesian_idxs(grid::RegularGrid) = CartesianIndices(length.(grid.grid))

# Given a cartesian index, the output of `basin_cell_index`, create
# a `Rect` object that corresponds to that grid cell!
function cell_index_to_rect(n::CartesianIndex, grid::RegularGrid)
    x = grid.grid[1][n[1]]
    y = grid.grid[2][n[2]]
    rect = Rect(x, y, grid.grid_steps[1], grid.grid_steps[2])
    return rect
end


color_obs = initialize_cells2!(ax, grid_nfo)

fig

for i in 1:5
    obs = color_obs[i, 2]
    obs[] = Makie.RGB(0, 0, 1)
end

fig

# Works!!!