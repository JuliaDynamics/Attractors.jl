# [Animation illustrating `AttractorsViaRecurrences`](@id recurrences_animation)

The following Julia script inputs a 2D continuous time dynamical system and animates its time evolution while illustrating how [`AttractorsViaRecurrences`](@ref) works.

```julia
using Attractors, CairoMakie
using PredefinedDynamicalSystems
using OrdinaryDiffEqVerner

# Set up dynamical system: bi-stable predator pray
function predator_prey_rule(u, p, t)
    r, c, μ, ν, α, β, χ, δ = p
    N, P = u
    common = α*N*P/(β+N)
    dN = r*N*(1 - (c/r)*N)*((N-μ)/(N+ν)) - common
    dP = χ*common - δ*P
    return SVector(dN, dP)
end

u0 = SVector(8.0, 0.01)
r = 2.0
# r, c, μ, ν, α, β, χ, δ = p
p = [r, 0.19, 0.03, 0.003, 800, 1.5, 0.004, 2.2]

diffeq = (alg = Rodas5P(), abstol = 1e-9, rtol = 1e-9)
ds = CoupledODEs(predator_prey_rule, u0, p; diffeq)

u0s = [ # animation will start from these initial conditions
    [10, 0.012],
    [15, 0.02],
    [12, 0.01],
    [13, 0.015],
    [5, 0.02],
]

density = 31
xg = range(-0.1, 20; length = density)
yg = range(-0.001, 0.03; length = density)
Δt = 0.1
grid = (xg, yg)
mapper = AttractorsViaRecurrences(ds, grid;
    Δt, consecutive_attractor_steps = 10, consecutive_basin_steps = 10, sparse = false,
    consecutive_recurrences = 100, attractor_locate_steps = 100,
)

##########################################################################

function animate_attractors_via_recurrences(
        mapper::AttractorsViaRecurrences, u0s;
        colors = ["#FFFFFF", "#7143E0","#0A9A84","#AF9327","#791457", "#6C768C", "#4287f5",],
        filename = "recurrence_algorithm.mp4",
    )

    grid_nfo = mapper.bsn_nfo.BoA.grid

    fig = Figure()
    ax = Axis(fig[1,1])

    # Populate the grid with poly! rectangle plots. However! The rectangles
    # correspond to the same "cells" of the grid. Additionally, all
    # rectangles are colored with an _observable_, that can be accessed
    # later using the `basin_cell_index` function. The observable
    # holds the face color of the rectangle!

    # Only 6 colors; need 3 for base, and extra 2 for each attractor.
    # will choose initial conditions that are only in the first 2 attractors
    COLORS = map(c -> Makie.RGBA(Makie.RGB(to_color(c)), 0.9), colors)

    function initialize_cells2!(ax, grid; kwargs...)
        # These are all possible outputs of the `basin_cell_index` function
        idxs = all_cartesian_idxs(grid)
        color_obs = Matrix{Any}(undef, size(idxs)...)
        # We now need to reverse-engineer
        for i in idxs
            rect = cell_index_to_rect(i, grid)
            color = Observable(COLORS[1])
            color_obs[i] = color
            poly!(ax, rect; color = color, strokecolor = :black, strokewidth = 0.5)
        end
        # Set the axis limits better
        mini, maxi = Attractors.minmax_grid_extent(grid)
        xlims!(ax, mini[1], maxi[1])
        ylims!(ax, mini[2], maxi[2])
        return color_obs
    end

    all_cartesian_idxs(grid::Attractors.RegularGrid) = CartesianIndices(length.(grid.grid))

    # Given a cartesian index, the output of `basin_cell_index`, create
    # a `Rect` object that corresponds to that grid cell!
    function cell_index_to_rect(n::CartesianIndex, grid::Attractors.RegularGrid)
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

    # This function gives correct color to search, recurrence, and
    # the individual attractors. Ignores the lost state.
    function update_current_cell_color!(cellcolor, bsn_nfo)
        # We only alter the cell color at specific situations
        state = bsn_nfo.state
        if state == :att_search
            if cellcolor[] == COLORS[1] # empty
                cellcolor[] = COLORS[2] # visited
            elseif cellcolor[] == COLORS[2] # visited
                cellcolor[] = COLORS[3] # recurrence
            end
        elseif state == :att_found
            attidx = (bsn_nfo.current_att_label ÷ 2)
            attlabel = (attidx - 1)*2 + 1
            cellcolor[] = COLORS[3+attlabel]
        end
        return
    end

    # Iteration and labelling
    ds = mapper.ds
    bsn_nfo = mapper.bsn_nfo
    u0 = current_state(ds)

    traj = Observable(SVector{2, Float64}[u0])
    point = Observable([u0])

    marker = Observable(:circle)
    lines!(ax, traj; color = :black, linewidth = 1)
    scatter!(ax, point; color = (:black, 0.5), markersize = 20, marker, strokewidth = 1.0, strokecolor = :black)

    stateobs = Observable(:att_search)
    consecutiveobs = Observable(0)
    labeltext = @lift("state: $($(stateobs))\nconsecutive: $($(consecutiveobs))")

    Label(fig[0, 1][1,1], labeltext; justification = :left, halign = :left, tellwidth = false)
    # add text  with options
    kwargstext = prod("$(p[1])=$(p[2])\n" for p in mapper.kwargs)
    Label(fig[0, 1][1, 2], kwargstext; justification = :right, halign = :right, tellwidth = false)


    # make legend
    entries = [PolyElement(color = c) for c in COLORS[2:end]]
    labels = ["visited", "recurrence", "attr. 1", "basin 1", "attr. 2", "basin 2"]

    Legend(fig[:, 2][1, 1], entries, labels)


    # %% loop
    # The following code is similar to the source code of `recurrences_map_to_label!`

    cell_label = 0
    record(fig, filename) do io
        for u0 in u0s
            reinit!(ds, copy(u0))
            traj[] = [copy(u0)]

            while cell_label == 0
                step!(ds, bsn_nfo.Δt)
                u = current_state(ds)

                # update FSM
                n = Attractors.basin_cell_index(u, bsn_nfo.BoA.grid)
                cell_label = Attractors.finite_state_machine!(bsn_nfo, n, u; mapper.kwargs...)

                state = bsn_nfo.state

                if cell_label ≠ 0 # FSM terminated; we assume no lost/divergence in the system
                    stateobs[] = :terminated

                    # color-code initial condition if we converged to attractor
                    # or to basin (even or odd cell label)
                    u0n = Attractors.basin_cell_index(u0, bsn_nfo.BoA.grid)

                    basidx = (cell_label - 1)
                    color_obs[u0n][] = COLORS[3+basidx]

                    # Clean up: all "visited" cells become white again
                    visited_idxs = findall(v -> (v[] == COLORS[2] || v[] == COLORS[3]), color_obs)
                    for n in visited_idxs
                        color_obs[n][] = COLORS[1] # empty
                    end
                    # clean up trajectory line
                    traj[] = []

                    for i in 1:15; recordframe!(io); end
                    cell_label = 0
                    break
                end

                # update visuals:
                point[] = [u]
                push!(traj[], u)
                notify(traj)
                marker[] = state2marker[state]
                stateobs[] = state
                consecutiveobs[] = bsn_nfo.consecutive_match

                update_current_cell_color!(color_obs[n], bsn_nfo)

                recordframe!(io)
            end
        end
    end
end

animate_attractors_via_recurrences(mapper, u0s)
```

```@raw html
<video width="75%" height="auto" controls autoplay loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/attractors/recurrence_algorithm.mp4?raw=true" type="video/mp4">
</video>
```
