#####################################################################################
# Implementation of the Finite State Machine (low level code)
#####################################################################################
"""
    recurrences_map_to_label!(bsn_nfo::BasinsInfo, ds, u0; kwargs...) -> ic_label

Return the label of the attractor that the initial condition `u0` converges to,
or `-1` if it does not convergence anywhere (e.g., divergence to infinity or exceeding
`mx_chk_safety`).

Notice the numbering system `cell_label` is as in `finite_state_machine!`.
Even numbers are attractors, odd numbers are basins.
"""
function recurrences_map_to_label!(bsn_nfo::BasinsInfo, ds::DynamicalSystem, u0;
        mx_chk_safety = Int(1e6), Ttr = 0, kwargs...
    )
    # This routine identifies the attractor using the previously defined basin.

    # reinitialize everything
    reinit!(ds, u0)
    Ttr > 0 && step!(ds, Ttr)
    reset_basins_counters!(bsn_nfo)
    cell_label = 0
    bsn_nfo.safety_counter = 0

    while cell_label == 0
        step!(ds, bsn_nfo.Δt)

        # This clause here is added because sometimes the algorithm will never hault
        # for e.g., an ill conditioned grid where two or more attractors intersect
        # within the same grid cell. In such a case, when starting on the second attractor
        # the trajectory will forever reset between locating a new attractor and recurring
        # on the previously found one...
        bsn_nfo.safety_counter += 1
        if bsn_nfo.safety_counter ≥ mx_chk_safety
            # TODO: Set up some debugging framework here via environment variable
            # @warn """
            # `AttractorsViaRecurrences` algorithm exceeded safety count without haulting.
            # It may be that the grid is not fine enough and attractors intersect in the
            # same cell, or `mx_chk_safety` is not high enough for a very fine grid.
            # Here are some info on current status:\n
            # state: $(current_state(ds)),\n
            # parameters: $(current_parameters(ds)).
            # """
            cleanup_visited_cells!(bsn_nfo)
            return -1
        end

        if !successful_step(ds)
            cleanup_visited_cells!(bsn_nfo)
            return -1
        end

        new_y = current_state(ds)
        # The internal function `_possibly_reduced_state` exists solely to
        # accommodate the special case of a Poincare map with the grid defined
        # directly on the hyperplane, `plane::Tuple{Int, <: Real}`.
        y = _possibly_reduced_state(new_y, ds, bsn_nfo)
        n = basin_cell_index(y, bsn_nfo.grid_nfo)
        cell_label = finite_state_machine!(bsn_nfo, n, y; kwargs...)
    end
    return cell_label
end

_possibly_reduced_state(y, ds, grid) = y
function _possibly_reduced_state(y, ds::PoincareMap, bsn_nfo)
    grid = bsn_nfo.grid
    if ds.planecrossing.plane isa Tuple && length(grid) == dimension(ds)-1
        return y[ds.diffidxs]
    else
        return y
    end
end



"""
Main procedure which performs one step of the finite state machine.
Directly implements the algorithm of Datseris & Wagemakers 2021,
see the flowchart (Figure 2).
I.e., it obtains current state, updates state with respect to the cell
code the trajectory it is currently on, performs the check of counters,
and sets the f.s.m. to the new state, or increments counter if same state.

The basins and attractors are coded in the array with odd numbers for the basins and
even numbers for the attractors. The attractor `2n` has the corresponding basin `2n+1`.
This codification is changed when the basins and attractors are returned to the user.
Diverging trajectories and the trajectories staying outside the grid are coded with -1.

The label `1` (initial value) outlined in the paper is `0` here instead.
The function returns `0` unless the FSM has terminated its operation.
"""
function finite_state_machine!(
        bsn_nfo::BasinsInfo, n::CartesianIndex, u;
        mx_chk_att = 2, mx_chk_hit_bas = 10, mx_chk_fnd_att = 1000, mx_chk_loc_att = mx_chk_fnd_att÷10,
        horizon_limit = 1e6, mx_chk_lost = 20, store_once_per_cell = true,
        show_progress = true, # show_progress can be used when finding new attractor.
    )

    # if n[1] == -1 means we are outside the grid,
    # otherwise, we retrieve the label stored at the grid
    # (which by default is 0 unless we have visited the cell before)
    ic_label = n[1] == -1 ? -1 : bsn_nfo.basins[n]

    update_finite_state_machine!(bsn_nfo, ic_label)

    # This state means that we have visited a cell that contains a recorded attractor
    if bsn_nfo.state == :att_hit
        if ic_label == bsn_nfo.prev_label
             bsn_nfo.consecutive_match += 1
        end
        if bsn_nfo.consecutive_match ≥ mx_chk_att
            # We've hit an existing attractor `mx_chk_att` times in a row
            # so we assign the initial condition directly to the attractor
            hit_att = ic_label + 1
            cleanup_visited_cells!(bsn_nfo)
            reset_basins_counters!(bsn_nfo)
            return hit_att
        end
        bsn_nfo.prev_label = ic_label
        return 0
    end

    # This state is "searching for an attractor". It is the initial state.
    if bsn_nfo.state == :att_search
        if ic_label == 0
            # unlabeled box, label it with current odd label and reset counter
            bsn_nfo.basins[n] = bsn_nfo.visited_cell_label
            # also keep track of visited cells. This makes it easier to clean
            # up the basin array later!
            push!(bsn_nfo.visited_cells, n) # keep track of visited cells
            bsn_nfo.consecutive_match = 1
        elseif ic_label == bsn_nfo.visited_cell_label
            # hit a previously visited box with the current label, possible attractor?
            bsn_nfo.consecutive_match += 1
        end

        # If we accummulated enough recurrences, we claim that we
        # have found an attractor, and we switch to `:att_found`.
        if bsn_nfo.consecutive_match >= mx_chk_fnd_att
            bsn_nfo.basins[n] = bsn_nfo.current_att_label
            store_attractor!(bsn_nfo, u)
            bsn_nfo.state = :att_found
            bsn_nfo.consecutive_match = 1
        end
        bsn_nfo.prev_label = ic_label
        return 0
    end

    # This state can only be reached from `:att_search`. It means we have
    # enough recurrences to claim we have found an attractor.
    # We then locate the attractor by recording enough cells.
    if bsn_nfo.state == :att_found
        if ic_label == bsn_nfo.current_att_label
            # Visited a cell already labelled as new attractor; check `store_once_per_cell`
            store_once_per_cell || store_attractor!(bsn_nfo, u)
        elseif ic_label == 0 || ic_label == bsn_nfo.visited_cell_label
            # Visited a cells that was not labelled as the new attractor
            # label it and store it as part of the attractor
            bsn_nfo.basins[n] = bsn_nfo.current_att_label
            store_attractor!(bsn_nfo, u)
        elseif even(ic_label) && ic_label != bsn_nfo.current_att_label
            # Visited a cell labelled as an *existing* attractor! We have
            # attractors intersection in the grid! The algorithm can't handle this,
            # so we throw an error.
            error("""
            During the phase of locating a new attractor, found via sufficient recurrences,
            we encountered a cell of a previously-found attractor. This means that two
            attractors intersect in the grid, or that the precision with which we find
            and store attractors is not fine enough. Either decrease the grid spacing,
            or increase `mx_chk_fnd_att` (or both).

            Index of cell that this occured at: $(n).
            """)
        end
        # in the `:att_found` phase, the consecutive match is always increasing
        bsn_nfo.consecutive_match += 1
        if bsn_nfo.consecutive_match ≥ mx_chk_loc_att
            # We have recorded the attractor with sufficient accuracy.
            # We now set the empty counters for the new attractor.
            cleanup_visited_cells!(bsn_nfo)
            bsn_nfo.visited_cell_label += 2
            bsn_nfo.current_att_label += 2
            reset_basins_counters!(bsn_nfo)
            # We return the label corresponding to the *basin* of the attractor
            return ic_label + 1
        end
        return 0
    end

    # We've hit a cell labelled as a basin of an existing attractor.
    # Note, this clause can only trigger if the basin array is NOT sparse.
    if bsn_nfo.state == :bas_hit
        if bsn_nfo.prev_label == ic_label
            bsn_nfo.consecutive_match += 1
        else
            bsn_nfo.consecutive_match = 1
        end
        if  bsn_nfo.consecutive_match > mx_chk_hit_bas
            cleanup_visited_cells!(bsn_nfo)
            reset_basins_counters!(bsn_nfo)
            return ic_label
        end
        bsn_nfo.prev_label = ic_label
        return 0
    end

    # This state occurs when the dynamical system state is outside the grid
    if bsn_nfo.state == :lost
        bsn_nfo.consecutive_lost += 1
        if bsn_nfo.consecutive_lost > mx_chk_lost || norm(u) > horizon_limit
            cleanup_visited_cells!(bsn_nfo)
            reset_basins_counters!(bsn_nfo)
            # problematic IC: diverges or wanders outside the defined grid
            return -1
        end
        bsn_nfo.prev_label = ic_label
        return 0
    end
end

function store_attractor!(bsn_nfo::BasinsInfo{D, G, Δ, T}, u) where {D, G, Δ, T}
    # bsn_nfo.current_att_label is the number of the attractor multiplied by two
    attractor_id = bsn_nfo.current_att_label ÷ 2
    V = SVector{D, T}
    if haskey(bsn_nfo.attractors, attractor_id)
        push!(bsn_nfo.attractors[attractor_id], V(u))
    else
        # initialize container for new attractor
        bsn_nfo.attractors[attractor_id] = StateSpaceSet([V(u)])
    end
end

# Notice that seting a basin index to 0 _deletes the index_ if the
# array is a `SparseArray`, see the source code file!
function cleanup_visited_cells!(bsn_nfo::BasinsInfo)
    old_label = bsn_nfo.visited_cell_label
    while !isempty(bsn_nfo.visited_cells)
        ind = pop!(bsn_nfo.visited_cells)
        if bsn_nfo.basins[ind] == old_label
            bsn_nfo.basins[ind] = 0 # 0 is the unvisited label / empty label
        end
    end
end


function reset_basins_counters!(bsn_nfo::BasinsInfo)
    bsn_nfo.consecutive_match = 0
    bsn_nfo.consecutive_lost = 0
    bsn_nfo.prev_label = 0
    bsn_nfo.state = :att_search
end

function update_finite_state_machine!(bsn_nfo, ic_label)
    current_state = bsn_nfo.state
    if current_state == :att_found
        # this is a terminal state, once reached you don't get out
        return
    end

    next_state = :undef
    if ic_label == 0 || ic_label == bsn_nfo.visited_cell_label
        # unlabeled box or previously visited box with the current label
        next_state = :att_search
    elseif iseven(ic_label)
        # hit an attractor box
        next_state = :att_hit
    elseif ic_label == -1
        # out of the grid we do not reset the counter of other state
        # since the trajectory can follow an attractor that spans outside the grid
        bsn_nfo.state = :lost
        return
    elseif isodd(ic_label)
        # hit an basin box
        next_state = :bas_hit
    end

    if next_state != current_state
        # reset counter except in lost state (the counter freezes in this case)
        if current_state == :lost
            bsn_nfo.consecutive_lost = 1
        else
            bsn_nfo.consecutive_match = 1
        end
    end
    bsn_nfo.state = next_state
    return
end

