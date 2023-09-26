include("sparse_arrays.jl")
#####################################################################################
# Type definition and documentation
#####################################################################################
"""
    AttractorsViaRecurrences(ds::DynamicalSystem, grid::Tuple; kwargs...)

Map initial conditions of `ds` to attractors by identifying attractors on the fly based on
recurrences in the state space, as outlined by Datseris & Wagemakers [Datseris2022](@cite).

`grid` is a tuple of ranges partitioning the state space so that a finite state
machine can operate on top of it. For example
`grid = (xg, yg)` where `xg = yg = range(-5, 5; length = 100)` for a two-dimensional
system with regular grid or one can also define an irregular grid, for example 
`xg = vcat(range(-5, -2; length = 50), range(-2, 5; length = 50)),
yg = range(-5, 5; length = 100)`. The grid has to be the same dimensionality as
the state space, use a [`ProjectedDynamicalSystem`](@ref) if you 
want to search for attractors in a lower dimensional subspace.
## Keyword arguments

* `sparse = true`: control the interval representation of the state space grid. If true,
   uses a sparse array, whose memory usage is in general more efficient than a regular
   array obtained with `sparse=false`. In practice, the sparse representation should
   always be preferred when searching for [`basins_fractions`](@ref). Only for very low
   dimensional systems and for computing the full [`basins_of_attraction`](@ref) the
   non-sparse version should be used.

### Time evolution configuration

* `Ttr = 0`: Skip a transient before the recurrence routine begins.
* `Δt`: Approximate integration time step (second argument of the `step!` function).
  The keyword `Dt` can also be used instead if `Δ` (`\\Delta`) is not accessible.
  It is `1` for discrete time systems.
  For continuous systems, an automatic value is calculated using
  [`automatic_Δt_basins`](@ref). For very fine grids, this can become very small,
  much smaller than the typical integrator internal step size in case of adaptive
  integrators. In such cases, use `force_non_adaptive = true`.
* `force_non_adaptive = false`: Only used if the input dynamical system is `CoupledODEs`.
  If `true` the additional keywords `adaptive = false, dt = Δt` are given as `diffeq`
  to the `CoupledODEs`. This means that adaptive integration is turned off and `Δt` is
  used as the ODE integrator timestep. This is useful in (1) very fine grids, and (2)
  if some of the attractors are limit cycles. We have noticed that in this case the
  integrator timestep becomes commensurate with the limit cycle period, leading to
  incorrectly counting the limit cycle as more than one attractor.


### Finite state machine configuration

* `mx_chk_att = 2`: Μaximum checks of consecutives hits of an existing attractor cell
  before declaring convergence to that existing attractor.
* `mx_chk_hit_bas = 10`: Maximum check of consecutive visits of the same basin of
  attraction before declaring convergence to an existing attractor.
* `mx_chk_fnd_att = 100`: Maximum check of consecutive visits to a previously visited
  unlabeled cell before declaring we have found a new attractor.
* `mx_chk_loc_att = 100`: Maximum check of consecutive visits to cells marked as a new
  attractor, during the attractor identification phase, before declaring we that we have
  identified the new attractor with sufficient cells.
* `store_once_per_cell = true`: Control if multiple points in state space that belong to
  the same cell are stored or not in the attractor, after an attractor is found.
  If `true`, each visited cell will only store a point once, which is desirable for fixed
  points and limit cycles. If `false`, at least `mx_chk_loc_att` points are
  stored per attractor, leading to more densely stored attractors,
  which may be desirable for instance in chaotic attractors.
* `mx_chk_lost = 20`: Maximum check of iterations outside the defined grid before we
  declare the orbit lost outside and hence assign it label `-1`.
* `horizon_limit = 1e6`: If the norm of the integrator state reaches this
  limit we declare that the orbit diverged to infinity.
* `mx_chk_safety = Int(1e6)`: A safety counter that is always increasing for
  each initial condition. Once exceeded, the algorithm assigns `-1` and throws a warning.
  This clause exists to stop the algorithm never haulting for innappropriate grids,
  where a found attractor may intersect in the same cell with a new attractor the orbit
  traces (which leads to infinite resetting of all counters).


## Description

An initial condition given to an instance of `AttractorsViaRecurrences` is iterated
based on the integrator corresponding to `ds`. A recurrence in the state space means
that the trajectory has converged to an attractor. This is the basis for finding attractors.

A finite state machine (FSM) follows the
trajectory in the state space, and constantly maps it to the given `grid`. The FSM
decides when an initial condition has successfully converged into an attractor. An array,
internally called "basins", stores the state of the FSM on the grid, according to the
indexing system described in [Datseris2022](@cite). As the system is integrated more and more,
the information of the "basins" becomes richer and richer with more identified attractors
or with grid cells that belong to basins of already found attractors.
Notice that only in the special method
`basins_of_attraction(mapper::AttractorsViaRecurrences)` the information of the
attraction or exit basins is utilized. In other functions like `basins_fractions`
only the attractor locations are utilized, as the basins themselves are not stored.

The iteration of a given initial condition continues until one of the following happens:
-  The trajectory hits `mx_chk_fnd_att` times in a row grid cells previously visited:
   it is considered that an attractor is found and is labelled with a new ID. Then,
   iteration continues a bit more until we have identified the attractor with sufficient
   accuracy, i.e., until `mx_chk_loc_att` cells with the new ID have been visited.
-  The trajectory hits an already identified attractor `mx_chk_att` consecutive times:
   the initial condition is numbered with the attractor's ID.
-  The trajectory hits a known basin `mx_chk_hit_bas` times in a row: the initial condition
   belongs to that basin and is numbered accordingly. Notice that basins are stored and
   used only when `sparse = false`.
-  The trajectory spends `mx_chk_lost` steps outside the defined grid or the norm
   of the integrator state becomes > than `horizon_limit`: the initial
   condition's label is set to `-1`.
-  If none of the above happens, the initial condition is labelled `-1` after
   and `mx_chk_safety` integrator steps.
"""
struct AttractorsViaRecurrences{DS<:DynamicalSystem, B, G, K} <: AttractorMapper
    ds::DS
    bsn_nfo::B
    grid::G
    kwargs::K
end

function AttractorsViaRecurrences(ds::DynamicalSystem, grid;
        Dt = nothing, Δt = Dt, sparse = true, force_non_adaptive = false,  kwargs...
    )


    if grid isa Tuple  # regular or irregular
        D = length(grid)
        if all(t -> t isa AbstractRange, grid) # regular
            grid_steps = SVector{D,Float64}(step.(grid))
            grid_maxima = SVector{D,Float64}(maximum.(grid))
            grid_minima = SVector{D,Float64}(minimum.(grid))
            
            finalgrid = RegularGrid(grid_steps, grid_minima, grid_maxima, grid)
            
        elseif any(t -> t isa AbstractVector, grid) && all(axis -> issorted(axis), grid) # irregular
            finalgrid = IrregularGrid(grid)
        else
            error("Incorrect grid specification")
        end
    elseif grid isa SubdivisionBasedGrid
        
        finalgrid = grid
    else
        error("Incorrect grid specification")
    end
    
    bsn_nfo = initialize_basin_info(ds, finalgrid, Δt, sparse)
    if ds isa CoupledODEs && force_non_adaptive
        newdiffeq = (ds.diffeq..., adaptive = false, dt = bsn_nfo.Δt)
        ds = CoupledODEs(ds, newdiffeq)
    end
    return AttractorsViaRecurrences(ds, bsn_nfo, finalgrid, kwargs)
end

function (mapper::AttractorsViaRecurrences)(u0; show_progress = true)
    # Call low level code of `basins_of_attraction` function. Notice that in this
    # call signature the interal basins info array of the mapper is NOT updated.
    lab = get_label_ic!(mapper.bsn_nfo, mapper.ds, u0; show_progress, mapper.kwargs...)
    # Transform to integers indexing from odd-even indexing
    return iseven(lab) ? (lab ÷ 2) : (lab - 1) ÷ 2
end

function Base.show(io::IO, mapper::AttractorsViaRecurrences)
    ps = generic_mapper_print(io, mapper)
    println(io, rpad(" system: ", ps), nameof(typeof(mapper.ds)))
    println(io, rpad(" grid: ", ps), mapper.grid)
    println(io, rpad(" attractors: ", ps), mapper.bsn_nfo.attractors)
    return
end

extract_attractors(m::AttractorsViaRecurrences) = m.bsn_nfo.attractors

"""
    basins_of_attraction(mapper::AttractorsViaRecurrences; show_progress = true)

This is a special method of `basins_of_attraction` that using recurrences does
_exactly_ what is described in the paper by Datseris & Wagemakers [Datseris2022](@ref).
By enforcing that the internal grid of `mapper` is the same as the grid of initial
conditions to map to attractors, the method can further utilize found exit and attraction
basins, making the computation faster as the grid is processed more and more.
"""
function basins_of_attraction(mapper::AttractorsViaRecurrences; show_progress = true)
    basins = mapper.bsn_nfo.basins
    if basins isa SparseArray;
        throw(ArgumentError("""
            Sparse version of AttractorsViaRecurrences is incompatible with
            `basins_of_attraction(mapper)`."""
        ))
    end
    
    if (mapper.bsn_nfo.grid_nfo isa SubdivisionBasedGrid)
        grid = mapper.grid.max_grid
    else
        grid = mapper.grid.grid
    end
    
    I = CartesianIndices(basins)
    progress = ProgressMeter.Progress(
        length(basins); desc = "Basins of attraction: ", dt = 1.0
    )

    # TODO: Here we can have a slightly more efficient iteration by
    # iterating over `I` in different ways. In this way it always starts from the edge of
    # the grid, which is the least likely location for attractors. We need to
    # iterate I either randomly or from its center.
    for (k, ind) in enumerate(I)
        if basins[ind] == 0
            show_progress && ProgressMeter.update!(progress, k)
            y0 = generate_ic_on_grid(grid, ind)
            basins[ind] = get_label_ic!(
                mapper.bsn_nfo, mapper.ds, y0; show_progress, mapper.kwargs...
            )
        end
    end

    # remove attractors and rescale from 1 to max number of attractors
    ind = iseven.(basins)
    basins[ind] .+= 1
    basins .= (basins .- 1) .÷ 2
    return basins, mapper.bsn_nfo.attractors
end
 
#####################################################################################
# Definition of `BasinInfo` and initialization
#####################################################################################

abstract type Grid end

Base.@kwdef struct RegularGrid{D} <:Grid
    grid_steps::SVector{D, Float64}
    grid_minima::SVector{D, Float64}
    grid_maxima::SVector{D, Float64}
    grid::NTuple{D,AbstractRange}
end

struct IrregularGrid{D} <:Grid
    grid::NTuple{D,Vector{Float64}}
end

Base.@kwdef struct SubdivisionBasedGrid{D} <:Grid
    grid_steps::Dict{Int, Vector{Int}}
    grid_minima::SVector{D, Float64}
    grid_maxima::SVector{D, Float64}
    lvl_array::Union{Matrix{Int}, Array{Int, D}}
    grid::NTuple{D,AbstractRange}
    max_grid::NTuple{D,AbstractRange}
end

    

mutable struct BasinsInfo{D, G <:Grid, Δ, T, Q, A <: AbstractArray{Int, D}}
    basins::A # sparse or dense
    grid_nfo::G
    Δt::Δ
    state::Symbol
    current_att_label::Int
    visited_cell::Int
    consecutive_match::Int
    consecutive_lost::Int
    prev_label::Int
    safety_counter::Int
    attractors::Dict{Int, StateSpaceSet{D, T}}
    visited_list::Q
end



function initialize_basin_info(ds::DynamicalSystem, grid_nfo, Δtt, sparse)
    
    grid = grid_nfo.grid

    Δt = if isnothing(Δtt)
        if grid_nfo isa IrregularGrid
            throw(error("Provide the Dt for irregular grid"))
        end
        isdiscretetime(ds) ? 1 : automatic_Δt_basins(ds, grid)
    else
        Δtt
    end

    D = dimension(ds)
    T = eltype(current_state(ds))
    G = length(grid)
    if D ≠ G && (ds isa PoincareMap && G ∉ (D, D-1))
        error("Grid and dynamical system do not have the same dimension!")
    end
   
    D = length(grid)

    multiplier = 0
    if grid_nfo isa SubdivisionBasedGrid
        multiplier = maximum(keys(grid_nfo.grid_steps))
    end
    basins_array = if sparse
        SparseArray{Int}(undef, (map(length, grid ).*(2^multiplier)))
    else
        zeros(Int, (map(length, grid).*(2^multiplier))...)

    end
    bsn_nfo = BasinsInfo(
        basins_array,
        grid_nfo,
        Δt,
        :att_search,
        2,4,0,1,0,0,
        Dict{Int, StateSpaceSet{G, T}}(),
        Vector{CartesianIndex{G}}(),
    )
    
    reset_basins_counters!(bsn_nfo)
    return bsn_nfo
end



using LinearAlgebra: norm

"""
    automatic_Δt_basins(ds::DynamicalSystem, grid; N = 5000) → Δt

Calculate an optimal `Δt` value for [`basins_of_attraction`](@ref).
This is done by evaluating the dynamic rule `f` (vector field) at `N` randomly chosen
points of the grid. The average `f` is then compared with the diagonal length of a grid
cell and their ratio provides `Δt`.

Notice that `Δt` should not be too small which happens typically if the grid resolution
is high. It is okay for [`basins_of_attraction`](@ref) if the trajectory skips a few cells.
But if `Δt` is too small the default values for all other keywords such
as `mx_chk_hit_bas` need to be increased drastically.

Also, `Δt` that is smaller than the internal step size of the integrator will cause
a performance drop.
"""
function automatic_Δt_basins(ds, grid; N = 5000)
    isdiscretetime(ds) && return 1
    if ds isa ProjectedDynamicalSystem
        # TODO:
        error("Automatic Δt finding is not implemented for `ProjectedDynamicalSystem`.")
    end
    steps = step.(grid)
    s = sqrt(sum(x^2 for x in steps)) # diagonal length of a cell
    indices = CartesianIndices(length.(grid))
    random_points = [generate_ic_on_grid(grid, ind) for ind in rand(indices, N)]
    dudt = 0.0
    udummy = copy(current_state(ds))
    f, p = dynamic_rule(ds), current_parameters(ds)
    for point in random_points
        deriv = if !isinplace(ds)
            f(point, p, 0.0)
        else
            f(udummy, point, p, 0.0)
            udummy
        end
        dudt += norm(deriv)
    end
    Δt = 10*s*N/dudt
    return Δt
end



"""
    function subdivision_based_grid(ds::DynamicalSystem, grid; maxlevel = 4) 

Construct a grid structure `SubdivisionBasedGrid`  that can be directly passed 
as the `grid` argument to the `mapper` function.

The approach is designed for continuous time systems in which different areas of 
the state space flow may have significantly different velocity. In case of 
originally coarse grids it may lead to algorithm being stuck in some regions with
a small motion speed and false identification of attractors. To prevent this from
happening we provide an algorithm expansion to dynamically evaluate different regions speed
of motion to handle areas of the grid which should be more coarse or dense than others.

To achieve this, function make use of `make_irregular_array` which automatically constructs 
an array of discretization levels indices for grid originally specified by user with 
`grid` argument. Upon construction function automatically stores necessary parameters
to further adopt mapping of initial conditions to specific grid density levels.
"""
function subdivision_based_grid(ds::DynamicalSystem, grid; maxlevel = 4)
    lvl_array = make_irregular_array(ds, grid, maxlevel)
    unique_lvls = Set(lvl_array)
    grid_steps = Dict{Int, Vector{Int}}()
    for i in unique_lvls
        grid_steps[i] = [length(axis)*2^i for axis in grid]
    end
    D = length(grid)
    grid_maxima = SVector{D,Float64}(maximum.(grid))
    grid_minima = SVector{D,Float64}(minimum.(grid))

    function scale_axis(axis, multiplier)
        new_length = length(axis) * (2^multiplier)
        return range(first(axis), last(axis), length=new_length)
    end 
    multiplier = maximum(keys(grid_steps))
    scaled_axis = [scale_axis(axis, multiplier) for axis in grid]
    max_grid = Tuple(scaled_axis)

    return SubdivisionBasedGrid(grid_steps, grid_minima, grid_maxima, lvl_array, grid, max_grid)
end


function make_irregular_array(ds, grid, maxlevel = 4)
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
# Implementation of the Finite State Machine (low level code)
#####################################################################################
"""
    get_label_ic!(bsn_nfo::BasinsInfo, ds, u0; kwargs...) -> ic_label

Return the label of the attractor that the initial condition `u0` converges to,
or `-1` if it does not convergence anywhere (e.g., divergence to infinity or exceeding
`mx_chk_safety`).

Notice the numbering system `cell_label` is as in `_identify_basin_of_cell!`
so before the label processing done in e.g., `basins_of_attraction`.
"""
function get_label_ic!(bsn_nfo::BasinsInfo, ds::DynamicalSystem, u0;
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
            relabel_visited_cell!(bsn_nfo, bsn_nfo.visited_cell, 0)
            return -1
        end

        step!(ds, bsn_nfo.Δt)
        if !successful_step(ds)
            relabel_visited_cell!(bsn_nfo, bsn_nfo.visited_cell, 0)
            return -1
        end

        new_y = current_state(ds)
        # The internal function `_possibly_reduced_state` exists solely to
        # accommodate the special case of a Poincare map with the grid defined
        # directly on the hyperplane, `plane::Tuple{Int, <: Real}`.
        if (bsn_nfo.grid_nfo isa RegularGrid) || (bsn_nfo.grid_nfo isa SubdivisionBasedGrid)
            grid_min = bsn_nfo.grid_nfo.grid_minima
        else bsn_nfo.grid_nfo isa IrregularGrid
            grid_min = minimum.(bsn_nfo.grid_nfo.grid)
        end
        y = _possibly_reduced_state(new_y, ds, grid_min)
        n = basin_cell_index(y, bsn_nfo.grid_nfo)
        cell_label = _identify_basin_of_cell!(bsn_nfo, n, y; kwargs...)
        
    end
    return cell_label
end

_possibly_reduced_state(y, ds, grid) = y
function _possibly_reduced_state(y, ds::PoincareMap, grid)
    if ds.planecrossing.plane isa Tuple && length(grid) == dimension(ds)-1
        return y[ds.diffidxs]
    else
        return y
    end
end



"""
Main procedure. Directly implements the algorithm of Datseris & Wagemakers 2021,
see the flowchart (Figure 2).

The basins and attractors are coded in the array with odd numbers for the basins and
even numbers for the attractors. The attractor `2n` has the corresponding basin `2n+1`.
This codification is changed when the basins and attractors are returned to the user.
Diverging trajectories and the trajectories staying outside the grid are coded with -1.

The label `1` (initial value) outlined in the paper is `0` here instead.
"""
function _identify_basin_of_cell!(
        bsn_nfo::BasinsInfo, n::CartesianIndex, u;
        mx_chk_att = 2, mx_chk_hit_bas = 10, mx_chk_fnd_att = 100, mx_chk_loc_att = 100,
        horizon_limit = 1e6, mx_chk_lost = 20, store_once_per_cell = true,
        show_progress = true, # show_progress only used when finding new attractor.
    )

    #if n[1] == -1 means we are outside the grid
    ic_label = n[1] == -1 ? -1 : bsn_nfo.basins[n]

    check_next_state!(bsn_nfo, ic_label)

    if bsn_nfo.state == :att_hit
        if ic_label == bsn_nfo.prev_label
             bsn_nfo.consecutive_match += 1
        end
        if bsn_nfo.consecutive_match ≥ mx_chk_att
            # Wait if we hit the attractor a mx_chk_att times in a row just
            # to check if it is not a nearby trajectory
            hit_att = ic_label + 1
            relabel_visited_cell!(bsn_nfo, bsn_nfo.visited_cell, 0)
            reset_basins_counters!(bsn_nfo)
            return hit_att
        end
        bsn_nfo.prev_label = ic_label
        return 0
    end

    if bsn_nfo.state == :att_search
        if ic_label == 0
            # unlabeled box, label it with current odd label and reset counter
            bsn_nfo.basins[n] = bsn_nfo.visited_cell
            push!(bsn_nfo.visited_list, n) # keep track of visited cells
            bsn_nfo.consecutive_match = 1
        elseif ic_label == bsn_nfo.visited_cell
            # hit a previously visited box with the current label, possible attractor?
            bsn_nfo.consecutive_match += 1
        end

        if bsn_nfo.consecutive_match >= mx_chk_fnd_att
            bsn_nfo.basins[n] = bsn_nfo.current_att_label
            store_attractor!(bsn_nfo, u, show_progress)
            bsn_nfo.state = :att_found
            bsn_nfo.consecutive_match = 1
        end
        bsn_nfo.prev_label = ic_label
        return 0
    end

    if bsn_nfo.state == :att_found
        if ic_label == 0 || ic_label == bsn_nfo.visited_cell
            # Maybe chaotic attractor, periodic or long recursion.
            # label this box as part of an attractor
            bsn_nfo.basins[n] = bsn_nfo.current_att_label
            bsn_nfo.consecutive_match = 1
            store_attractor!(bsn_nfo, u, show_progress)
        elseif iseven(ic_label) && (bsn_nfo.consecutive_match <  mx_chk_loc_att)
            # We make sure we hit the attractor another mx_chk_loc_att consecutive times
            # just to be sure that we have the complete attractor
            bsn_nfo.consecutive_match += 1
            store_once_per_cell || store_attractor!(bsn_nfo, u, show_progress)
        elseif iseven(ic_label) && bsn_nfo.consecutive_match >= mx_chk_loc_att
            # We have checked the presence of an attractor: tidy up everything
            # and get a new cell
            relabel_visited_cell!(bsn_nfo, bsn_nfo.visited_cell, 0)
            # pick the next label for labeling the basin.
            bsn_nfo.visited_cell += 2
            bsn_nfo.current_att_label += 2
            reset_basins_counters!(bsn_nfo)
            return ic_label + 1;
        end
        return 0
    end

    if bsn_nfo.state == :bas_hit
        # hit a labeled basin point of the wrong basin, happens all the time,
        # we check if it happens mx_chk_hit_bas times in a row or if it happens
        # N times along the trajectory whether to decide if it is another basin.
        if bsn_nfo.prev_label == ic_label
            bsn_nfo.consecutive_match += 1
        else
            bsn_nfo.consecutive_match = 1
        end
        if  bsn_nfo.consecutive_match > mx_chk_hit_bas
            relabel_visited_cell!(bsn_nfo, bsn_nfo.visited_cell, 0)
            reset_basins_counters!(bsn_nfo)
            return ic_label
        end
        bsn_nfo.prev_label = ic_label
        return 0
    end

    if bsn_nfo.state == :lost
        bsn_nfo.consecutive_lost += 1  
        if bsn_nfo.consecutive_lost > mx_chk_lost || norm(u) > horizon_limit
            relabel_visited_cell!(bsn_nfo, bsn_nfo.visited_cell, 0)
            reset_basins_counters!(bsn_nfo)
            # problematic IC: diverges or wanders outside the defined grid
            return -1
        end
        bsn_nfo.prev_label = ic_label
        return 0
    end
end

function store_attractor!(bsn_nfo::BasinsInfo{D, G, Δ, T},
    u, show_progress = true) where {D, G, Δ, T}
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

function relabel_visited_cell!(bsn_nfo::BasinsInfo, old_label, new_label)
    while !isempty(bsn_nfo.visited_list)
        ind = pop!(bsn_nfo.visited_list)
        if bsn_nfo.basins[ind] == old_label
            bsn_nfo.basins[ind] = new_label
        end
    end
end




function basin_cell_index(y_grid_state, grid_nfo::RegularGrid)
    iswithingrid = true
    @inbounds for i in eachindex(grid_nfo.grid_minima)

        if !(grid_nfo.grid_minima[i] ≤ y_grid_state[i] ≤ grid_nfo.grid_maxima[i])
            iswithingrid = false
            break
        end
    end
    if iswithingrid
        # Snap point to grid
        ind = @. round(Int, (y_grid_state - grid_nfo.grid_minima)/grid_nfo.grid_steps) + 1
        return CartesianIndex(ind...)
    else
        return CartesianIndex(-1)
    end
end



function basin_cell_index(y_grid_state, grid_nfo::SubdivisionBasedGrid)


    D = length(grid_nfo.grid)
    initial_index = basin_cell_index(y_grid_state, RegularGrid(SVector{D,Float64}(step.(grid_nfo.grid)), grid_nfo.grid_minima, grid_nfo.grid_maxima, grid_nfo.grid))
    if initial_index == CartesianIndex(-1)
        return CartesianIndex(-1)
    end
    cell_area = grid_nfo.lvl_array[initial_index]
    grid_maxima = grid_nfo.grid_maxima
    grid_minima = grid_nfo.grid_minima
    grid_steps = grid_nfo.grid_steps
    max_level = maximum(keys(grid_steps))
    grid_step = (grid_maxima - grid_minima .+1)./ grid_steps[cell_area]
    iswithingrid = true
    @inbounds for i in eachindex(grid_minima)

        if !(grid_minima[i] ≤ y_grid_state[i] ≤ grid_maxima[i])
            iswithingrid = false
            break
        end
    end
    if iswithingrid
        ind = @. round(Int,(y_grid_state - grid_minima)/grid_step, RoundDown) * (2^(max_level-cell_area)) + 1
        return CartesianIndex(ind...) 
    else
        return CartesianIndex(-1)
    end  
end


function basin_cell_index(y_grid_state, grid_nfo::IrregularGrid)
    
    for axis in grid_nfo.grid
        if !issorted(axis)
            throw(error("Please provide sorted vectors for each axis"))
        end
    end

    for (axis, coord) in zip(grid_nfo.grid, y_grid_state)
        if coord < first(axis) || coord > last(axis)
            return CartesianIndex(-1)
        end
    end

    cell_indices = map((x, y) -> searchsortedlast(x, y), grid_nfo.grid, y_grid_state)
    
    return CartesianIndex(cell_indices...)
        
end


function reset_basins_counters!(bsn_nfo::BasinsInfo)
    bsn_nfo.consecutive_match = 0
    bsn_nfo.consecutive_lost = 0
    bsn_nfo.prev_label = 0
    bsn_nfo.state = :att_search
end

function check_next_state!(bsn_nfo, ic_label)
    next_state = :undef
    current_state = bsn_nfo.state
    if current_state == :att_found
        # this is a terminal state, once reached you don't get out
        return
    end

    if ic_label == 0 || ic_label == bsn_nfo.visited_cell
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
end

