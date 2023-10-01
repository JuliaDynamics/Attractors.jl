#####################################################################################
# Type definition and documentation
#####################################################################################
"""
    AttractorsViaRecurrences(ds::DynamicalSystem, grid::Tuple; kwargs...)

Map initial conditions of `ds` to attractors by identifying attractors on the fly based on
recurrences in the state space, as outlined by Datseris & Wagemakers [Datseris2022](@cite).

`grid` is instructions for partitioning the state space into finite-sized cells
so that a finite state machine can operate on top of it. Possibilities are:

1. A tuple of sorted `AbstractRange`s for a regular grid.
  Example is `grid = (xg, yg)` where `xg = yg = range(-5, 5; length = 100)`
  for a two-dimensional system
2. A tuple of sorted `AbstractVector`s for an irregular grid, for example
  `grid = (xg, yg)` with `xg = vcat(range(-5, -2; length = 50), range(-2, 5; length = 50)),
  yg = range(-5, 5; length = 100)`.
3. An instance of the special grid type
[`SubdividedBasedGrid`](@ref), which can be created either manually or by using [`subdivision_based_grid`](@ref).
  This automatically analyzes and adapts grid discretization
  levels in accordance with state space flow speed in different regions.

The grid has to be the same dimensionality as
the state space, use a [`ProjectedDynamicalSystem`](@ref) if you
want to search for attractors in a lower dimensional subspace.

## Keyword arguments

* `sparse = true`: control the storage type of the state space grid. If true,
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
  This is ignored if `sparse = true`, as basins are not stored internally.
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
        Dt = nothing, Δt = Dt, sparse = true, force_non_adaptive = false, kwargs...
    )

    if grid isa Tuple  # regular or irregular
        if all(t -> t isa AbstractRange, grid) && all(axis -> issorted(axis), grid) # regular
            finalgrid = RegularGrid(grid)
        elseif any(t -> t isa AbstractVector, grid) && all(axis -> issorted(axis), grid) # irregular
            finalgrid = IrregularGrid(grid)
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
    # Call low level code. Notice that in this
    # call signature the interal basins info array of the mapper is NOT updated
    # with the basins of attraction info. Only with the attractors info.
    lab = recurrences_map_to_label!(mapper.bsn_nfo, mapper.ds, u0; show_progress, mapper.kwargs...)
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
            basins[ind] = recurrences_map_to_label!(
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

mutable struct BasinsInfo{D, G<:Grid, Δ, T, A <: AbstractArray{Int, D}}
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
    visited_list::Vector{CartesianIndex{D}}
end

function initialize_basin_info(ds::DynamicalSystem, grid_nfo, Δtt, sparse)
    Δt = if isnothing(Δtt)
        isdiscretetime(ds) ? 1 : automatic_Δt_basins(ds, grid_nfo)
    else
        Δtt
    end

    grid = grid_nfo.grid # this is always a Tuple irrespectively of grid type
    D = dimension(ds)
    T = eltype(current_state(ds))
    G = length(grid)
    if D ≠ G && (ds isa PoincareMap && G ∉ (D, D-1))
        error("Grid and dynamical system do not have the same dimension!")
    end

    if grid_nfo isa SubdivisionBasedGrid
        multiplier = maximum(keys(grid_nfo.grid_steps))
    else
        multiplier = 0
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
points within the bounding box of the grid.
The average `f` is then compared with the average diagonal length of a grid
cell and their ratio provides `Δt`.

Notice that `Δt` should not be too small which happens typically if the grid resolution
is high. It is okay if the trajectory skips a few cells.
Also, `Δt` that is smaller than the internal step size of the integrator will cause
a performance drop.
"""
function automatic_Δt_basins(ds, grid_nfo::Grid; N = 5000)
    isdiscretetime(ds) && return 1
    if ds isa ProjectedDynamicalSystem
        error("Automatic Δt finding is not implemented for `ProjectedDynamicalSystem`.")
    end

    # Create a random sampler with min-max the grid range
    mins, maxs = minmax_grid_extent(grid_nfo)
    sampler, = statespace_sampler(HRectangle([mins...], [maxs...]))
    # Sample velocity at random points
    f, p, t0 = dynamic_rule(ds), current_parameters(ds), initial_time(ds)
    udummy = copy(current_state(ds))
    dudt = 0.0
    for _ in 1:N
        point = sampler()
        deriv = if !isinplace(ds)
            f(point, p, t0)
        else
            f(udummy, point, p, t0)
            udummy
        end
        dudt += norm(deriv)
    end

    s = mean_cell_diagonal(grid_nfo)
    Δt = 10*s*N/dudt
    return Δt
end

include("sparse_arrays.jl")
include("finite_state_machine.jl")
include("grids.jl")