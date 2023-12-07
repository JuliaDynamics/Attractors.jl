export edgetracking, bisect_to_edge, EdgeTrackingResults

"""
    EdgeTrackingResults(edge, track1, track2, time, bisect_idx)
Data type that stores output of the [`edgetracking`](@ref) algorithm.

## Fields
* `edge::StateSpaceSet`: the pseudo-trajectory representing the edge track
* `track1::StateSpaceSet`: the parallel pseudo-trajectory tracking the edge within basin 1
* `track2::StateSpaceSet`: the parallel pseudo-trajectory tracking the edge within basin 2
* `time::Vector`: time points of the above `StateSpaceSet`s
* `bisect_idx::Vector`: indices of `time` at which a re-bisection was done
"""
struct EdgeTrackingResults{D, T}
    edge::StateSpaceSet{D, T}
    track1::StateSpaceSet{D, T}
    track2::StateSpaceSet{D, T}
    time::Vector{Float64}
    bisect_idx::Vector{Int}
end

"""
    edgetracking(ds::DynamicalSystem, attractors::Dict; kwargs...)
Track along a basin boundary in a dynamical system `ds` with two or more attractors
in order to find an *edge state*. Results are returned in the form of
`EdgeTrackingResults`, which contains the pseudo-trajectory `edge` representing the track on
the basin boundary, along with additional output (see below).

The system's `attractors` are specified as a Dict of `StateSpaceSet`s, as in
[`AttractorsViaProximity`](@ref) or the output of `extract_attractors`. By default, the
algorithm is initialized from the first and second attractor in `attractors`. Alternatively,
the initial states can be set via keyword arguments `u1`, `u2` (see below). Note that the
two initial states must belong to different basins of attraction.

## Keyword arguments
* `bisect_thresh = 1e-7`: distance threshold for bisection
* `diverge_thresh = 1e-6`: distance threshold for parallel integration
* `u1`: first initial state (defaults to first point in first entry of `attractors`)
* `u2`: second initial state (defaults to first point in second entry of `attractors`)
* `maxiter = 100`: maximum number of iterations before the algorithm stops
* `abstol = 0.0`: distance threshold for convergence of the updated edge state
* `tmax = Inf`: maximum integration time of parallel trajectories until re-bisection 
* `Δt = 0.01`: time step passed to [`step!`](@ref) when evolving the two trajectories
* `ϵ_mapper = nothing`: `ϵ` parameter in [`AttractorsViaProximity`](@ref)
* `show_progress = true`: if true, shows progress bar and information while running
* `kwargs...`: additional keyword arguments to be passed to [`AttractorsViaProximity`](@ref)

## Description
The edge tracking algorithm is a numerical method to find
an *edge state* or (possibly chaotic) saddle on the boundary between two basins of
attraction. Introduced by [Battelino1988](@cite) and further described by
[Skufca2006](@cite), the
algorithm has been applied to, e.g., the laminar-turbulent boundary in plane Couette
flow [Schneider2008](@cite), Wada basins [Wagemakers2020](@cite), as well as Melancholia
states in conceptual [Mehling2023](@cite) and intermediate-complexity [Lucarini2017](@cite) 
climate models. 
Relying only on forward integration of the system, it works even in
high-dimensional systems with complicated fractal basin boundary structures.

The algorithm consists of two main steps: bisection and tracking. First, it iteratively 
bisects along a straight line in state space between the intial states `u1` and `u2` to find
the separating basin boundary. The bisection stops when the two updated states are less than
`bisect_thresh` (Euclidean distance in state space) apart from each other.
Next, a `ParallelDynamicalSystem` is initialized
from these two updated states and integrated forward until the two trajectories diverge
from each other by more than `diverge_thresh` (Euclidean distance). The two final states of
the parallel integration are then used as new states `u1` and `u2` for a new bisection, and 
so on, until a stopping criterion is fulfilled. 

Two stopping criteria are implemented via the keyword arguments `maxiter` and `abstol`.
Either the algorithm stops when the number of iterations reaches `maxiter`, or when the
state space position of the updated edge point changes by less than `abstol` (in
Euclidean distance) compared to the previous iteration. Convergence below `abstol` happens
after sufficient iterations if the edge state is a saddle point. However, the edge state
may also be an unstable limit cycle or a chaotic saddle. In these cases, the algorithm will
never actually converge to a point but (after a burn-in) continue populating the set
consituting the edge state by tracking along it.

A central idea behind this algorithm is that basin boundaries are typically the stable
manifolds of unstable sets, namely edge states or saddles. The flow along the basin boundary 
will thus lead to these sets, and the iterative bisection neutralizes the unstable
direction of the flow away from the basin boundary. If the system possesses multiple edge 
states, the algorithm will find one of them depending on where the initial bisection locates
the boundary.

## Output

Returns a data type `EdgeTrackingResults` containing five fields:
* `edge::StateSpaceSet`: the pseudo-trajectory representing the tracked edge segment
* `track1::StateSpaceSet`: the pseudo-trajectory tracking the edge within basin 1
* `track2::StateSpaceSet`: the pseudo-trajectory tracking the edge within basin 2
* `time::Vector`: time points of the above `StateSpaceSet`s
* `bisect_idx::Vector`: indices of `time` at which a re-bisection was done

Note that `edge` is simply the spatial mean (in state space) of `track1` and `track2`.
"""
function edgetracking(ds::DynamicalSystem, attractors::Dict;
    bisect_thresh=1e-6,
    diverge_thresh=1e-5,
    u1=collect(values(attractors))[1][1],
    u2=collect(values(attractors))[2][1],
    maxiter=100,
    abstol=0.0,
    T_transient=0.0,
    Δt=0.01,
    tmax=Inf,
    ϵ_mapper=nothing,
    show_progress=true,
    kwargs...)
    
    pds = ParallelDynamicalSystem(ds, [u1, u2])
    mapper = AttractorsViaProximity(ds, attractors, ϵ_mapper; kwargs...)
    
    edgetracking(pds, mapper;
        bisect_thresh, diverge_thresh, maxiter, abstol, T_transient, Δt, tmax, show_progress)
end;

"""
    edgetracking(pds::ParallelDynamicalSystem, mapper::AttractorMapper; kwargs...)
Low-level function for running the edge tracking algorithm, see [`edgetracking`](@ref)
for a description, keyword arguments and output type.

`pds` is a `ParallelDynamicalSystem` with two states. The `mapper` must be an
`AttractorMapper` of subtype `AttractorsViaProximity` or `AttractorsViaRecurrences`.
"""
function edgetracking(pds::ParallelDynamicalSystem, mapper::AttractorMapper;
    bisect_thresh=1e-6,
    diverge_thresh=1e-5,
    maxiter=100,
    abstol=0.0,
    T_transient=0.0,
    Δt=0.01,
    tmax=Inf,
    show_progress=true)
    
    if bisect_thresh >= diverge_thresh
        error("diverge_thresh must be larger than bisect_thresh.")
    end

    
    # initial bisection
    u1, u2 = bisect_to_edge(pds, mapper; bisect_thresh)
    edgestate = (u1 + u2)/2
    
    track1, track2 = [u1], [u2]
    time, bisect_idx = Float64[], Int[]
    progress = ProgressMeter.Progress(maxiter; desc = "Running edge tracking algorithm",
        enabled = show_progress)
    
    # edge track iteration loop
    displacement, counter, T = Inf, 1, 0.0
    while (displacement > abstol) && (maxiter > counter)
        t = 0
        set_state!(pds, u1, 1)
        set_state!(pds, u2, 2)
        distance = diffnorm(pds)
        # forward integration loop
        while (distance < diverge_thresh) && (t < tmax)
            step!(pds, Δt)
            distance = diffnorm(pds)
            t += Δt
            T += Δt
            if T >= T_transient
                push!(track1, current_state(pds, 1))
                push!(track2, current_state(pds, 2))
                push!(time, T)
            end
        end
        # re-bisect
        u1, u2 = bisect_to_edge(pds, mapper; bisect_thresh)
        T += Δt
        if T >= T_transient
            push!(track1, current_state(pds, 1))
            push!(track2, current_state(pds, 2))
            push!(time, T)
            push!(bisect_idx, length(time))
        end
        displacement = diffnorm(edgestate, (u1 + u2)/2)
        edgestate = (u1 + u2)/2
        counter += 1
        
        ProgressMeter.next!(progress;
            showvalues = [(:Iteration, counter), (:"Edge point", edgestate)])
        (counter == maxiter) && @warn("Reached maximum number of $(maxiter) iterations.")
    end

    (counter < maxiter) && println("Edge-tracking converged after $(counter) iterations.")
    
    track1 = StateSpaceSet(reduce(hcat, track1)')
    track2 = StateSpaceSet(reduce(hcat, track2)')

    return EdgeTrackingResults(
        StateSpaceSet((vec(track1) .+ vec(track2))./2),
        track1,
        track2,
        time,
        bisect_idx)
end;

"""
    bisect_to_edge(pds::ParallelDynamicalSystem, mapper::AttractorMapper; kwargs...) -> u1, u2
Finds the basin boundary between two states `u1, u2 = current_states(pds)` by bisecting
along a straight line in phase space. The states `u1` and `u2` must belong to different
basins.

Returns two new states `u1, u2` located on either side of the basin boundary that lie less
than `bisect_thresh` (Euclidean distance in state space) apart from each other.

## Keyword arguments
* `bisect_thresh = 1e-7`: The maximum (Euclidean) distance between the two returned states.

## Description
`pds` is a `ParallelDynamicalSystem` with two states. The `mapper` must be an
`AttractorMapper` of subtype `AttractorsViaProximity` or `AttractorsViaRecurrences`.

!!! info
    If the straight line between `u1` and `u2` intersects the basin boundary multiple
    times, the method will find one of these intersection points. If more than two attractors
    exist, one of the two returned states may belong to a different basin than the initial
    conditions `u1` and `u2`. A warning is raised if the bisection involves a third basin.
"""
function bisect_to_edge(pds::ParallelDynamicalSystem, mapper::AttractorMapper;
    bisect_thresh=1e-6)
    u1, u2 = current_states(pds)
    idx1, idx2 = mapper(u1), mapper(u2)
    
    if (idx1 == idx2)
        if idx1 == -1
            error("AttractorMapper returned label -1 (could not match the initial condition with any attractor).
            Try changing the settings of AttractorsViaProximity or increasing bisect_thresh, diverge_thresh.")
        else
            error("Both initial conditions belong to the same basin of attraction.
                Attractor label: $(idx1)
                u1 = $(u1)
                u2 = $(u2)")
        end
    end
    
    distance = diffnorm(u1, u2)
    while distance > bisect_thresh
        u_new = (u1 + u2)/2
        idx_new = mapper(u_new)
        # slightly shift u_new if it lands too close to the boundary
        retry_counter = 1
        while (idx_new == -1) && retry_counter < 3 # ToDO: make kwarg
            @warn "Shifting new point $(retry_counter)"
            u_new += bisect_thresh*(u1 - u2)
            idx_new = mapper(u_new)
            retry_counter += 1
        end
        # update u1 or u2
        if idx_new == idx1
            u1 = u_new
        else
            if idx_new != idx2
                if idx_new == -1
                    error("AttractorMapper returned label -1 (could not match the initial condition with any attractor.)
                    Try changing the settings of AttractorsViaProximity or increasing bisect_thresh, diverge_thresh.")
                else
                    @warn "New bisection point belongs to a third basin of attraction."
                end
            end
            u2 = u_new
        end
        distance = diffnorm(u1, u2)
    end
    return u1, u2
end;

function diffnorm(u1, u2)
    d = zero(eltype(u1))
    @inbounds for i in eachindex(u1)
        d += (u1[i] - u2[i])^2
    end
    sqrt(d)
end;

diffnorm(pds::ParallelDynamicalSystem) = diffnorm(current_state(pds, 1),
    current_state(pds, 2))