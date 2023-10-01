using Attractors, GLMakie
using PredefinedDynamicalSystems

# Set up dynamical system: (magnetic pendulum)
# ds = PredefinedDynamicalSystems.magnetic_pendulum(d=0.2, α=0.2, ω=0.8, N=3)
# ds = ProjectedDynamicalSystem(ds, [1, 2], [0.0, 0.0])
# set_state!(ds, [2.5, 3.1])
# xg = yg = range(-4, 4; length = 21)

# Set up dynamical system: tri-stable predator pray
function predator_prey_rule(u, p, t)
    r, c, μ, ν, α, β, χ, δ = p
    N, P = u
    common = α*N*P/(β+N)
    dN = r*N*(1 - (c/r)*N)*((N-μ)/(N+ν)) - common
    dP = χ*common - δ*P
    return SVector(dN, dP)
end

u0 = SVector(8.0, 0.01)
r = 1.6
# r, c, μ, ν, α, β, χ, δ = p
p = [r, 0.19, 0.03, 0.003, 800, 1.5, 0.004, 2.2]

Δt = 0.2
diffeq = (adaptive = false, dt = Δt)
ds = CoupledODEs(predator_prey_rule, u0, p; diffeq)

density = 41
xg = range(0, 15; length = density)
yg = range(0, 0.02; length = density)

grid = (xg, yg)
mapper = AttractorsViaRecurrences(ds, grid;
    Δt, mx_chk_hit_bas = 10, mx_chk_fnd_att = 10, mx_chk_loc_att = 10, sparse = false,
)
grid_nfo = mapper.bsn_nfo.grid_nfo

using Attractors: RegularGrid, IrregularGrid, SubdivisionBasedGrid

fig = Figure()
ax = Axis(fig[1,1])

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
        color = Observable(Makie.RGBA(0,0,0,0))
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
    dx = grid.grid_steps[1]
    dy = grid.grid_steps[2]
    rect = Rect(x - dx/2, y - dy/2, dx, dy)
    return rect
end

color_obs = initialize_cells2!(ax, grid_nfo)

# plot the trajectory
state2marker = Dict(
    :att_search => :circle,
    :att_found => :dtriangle,
    :att_hit => :rect,
    :lost => :star5,
    :bas_hit => :xcross,
)

COLORS = ["#7143E0","#0A9A84","#AF9327","#791457","#6C768C",]
state2color = Dict(
    :att_search => COLORS[1],
    :att_found => COLORS[2],
    :lost => COLORS[3],
    :bas_hit => COLORS[4],
)

# Iteration and labelling
ds = mapper.ds
bsn_nfo = mapper.bsn_nfo
u0 = current_state(ds)

traj = Observable(SVector{2, Float64}[u0])
point = Observable([u0])
marker = Observable(:circle)
lines!(ax, traj; color = :black, linewidth = 1)
scatter!(ax, point; color = (:black, 0.5), markersize = 20, marker, strokewidth = 1.0, strokecolor = :black)

display(fig)

# %% loop
# The following code is similar to the source code of `recurrences_map_to_label!`

cell_label = 0

while cell_label == 0
    step!(ds, bsn_nfo.Δt)
    u = current_state(ds)

    # update FSM
    n = Attractors.basin_cell_index(u, bsn_nfo.grid_nfo)
    cell_label = Attractors.finite_state_machine!(bsn_nfo, n, u; mapper.kwargs...)

    state = bsn_nfo.state

    if cell_label ≠ 0 # FSM terminated.
        break
    end

    # update visuals:
    point[] = [u]
    push!(traj[], u)
    notify(traj)
    marker[] = state2marker[state]

    # update color of cell
    if state != :att_hit # note that here we do not update, the cell has been colored
        newcolr = state2color[state]
        color_obs[n][] = Makie.RGBA(Makie.RGB(to_color(newcolr)), 0.75)
    end
end

display(fig)