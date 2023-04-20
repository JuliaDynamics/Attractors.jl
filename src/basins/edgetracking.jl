export edgetracking, bisect_to_edge

"""
    edgetracking(ds::DynamicalSystem, u1, u2, attractors::Dict; kwargs...)
Runs the edge tracking algorithm^[1,2,3] for a multistable DynamicalSystem `ds`.

The algorithm is initialized from two states `u1`, `u2` in state space that must belong to
different basins of attraction. The system's `attractors` are specified as a dictionary of
`StateSpaceSet`s.

## Description
The edge tracking algorithm is a simple numerical method to find the *edge state* or
(possibly chaotic) saddle on the boundary between two basins of attraction. It is first
described by Skufca et al. (2006)[^1], and has been applied to e.g. the laminar-turbulent
boundary in plane Couette flow[^2] as well as Melancholia states in a climate model
separating the snowball and warm climate states[^3]. Relying only on forward integration of
the system, it works even in high-dimensional systems with complicated basin boundary
structures.

The algorithm consists of two main steps: bisection and tracking. First, it iteratively 
bisects along a
straight line in state space between the intial states `u1` and `u2` to find the basin
boundary. The bisection stops when the two updated states are less than `eps1` (in
Euclidean distance) apart from each other. Next, a `ParallelDynamicalSystem` is initialized
from these two updated states and integrated forward until the two trajectories diverge
from each other by more than `eps2` (Euclidean distance). The two final states of the
parallel integration are then used as new states `u1` and `u2` for a new bisection, and 
so on, until a stopping criterion is fulfilled.

## Keyword arguments
* `eps1=1e-9`: bisection distance threshold
* `eps2=1e-8`: trajectory divergence distance threshold
* `maxiter=100`: maximum number of iterations before the algorithm stops
* `abstol=0.0`: convergence threshold for returned edge state (distance in state space)
* `tmax = Inf`: maximum integration time of parallel trajectories until re-bisection 
* `dt=0.01`: integration time step
* `ϵ_mapper=0.1`: `ϵ` parameter in [`AttractorsViaProximity`](@ref)
* `verbose=false`: if true, prints info while running
* `output_level`: what data to return (see below)
* `kwargs...`: additional keyword arguments to be passed to `AttractorsViaProximity`

## References
[^1]: [Skufca et al. (2006)](https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.96.174101?casa_token=RUn26KnFdNEAAAAA%3AoXsTlmEWVMkEYbOtR-j2PH2vYOOPOy1a2R_37ncnf4gsiHp6GR66M-IBpzXocLoMQC_oHhk8MIFRa_8)
[^2]: [Schneider et al. (2008)](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.78.037301?casa_token=mLQmTv_cBGUAAAAA%3AVKnQs290sq1MNm-5k8hW7nJeLtVX54I7l-SEGol_HUPSCwziPi-EGDE8ucrDiVMIXZGUbzzam8benFw)
[^3]: [Luarini and Bodai (2017)](https://iopscience.iop.org/article/10.1088/1361-6544/aa6b11)

## Output 
Output can be controlled via the `output_level` argument.
* `output_level = 0`: returns only the final state on the edge
* `output_level = 1`: returns the track along the edge as a vector of states
* `output_level = 2`: returns `[track, bisect1, bisect2]`, where `track` is the track along  
  the edge and `bisect1`, `bisect2` are vectors containing the bisection points on either  
  side of the basin boundary

!!! warning
    May behave erroneously when run with `solver = SimpleATsit5()`, which is the default
    solver for `AttractorsViaProximity`. The recommended solver is `Vern9()`.
"""
function edgetracking(ds::DynamicalSystem, u1, u2, attractors::Dict;
    eps1 = 1e-9,
    eps2 = 1e-8,
    maxiter = 100,
    abstol = 0.0,
    ϵ_mapper = 0.1,
    #diffeq = (;alg = Vern9()),
    dt = 0.01,
    tmax = Inf,
    output_level=2,
    verbose = false,
    kwargs...)
    
    pds = ParallelDynamicalSystem(ds, [u1, u2])
    mapper = AttractorsViaProximity(ds, attractors, ϵ_mapper; kwargs...)
    
    edgetracking(pds, mapper;
        eps1=eps1, eps2=eps2, maxiter=maxiter, abstol=abstol, dt=dt, tmax=tmax,
        output_level=output_level, verbose=verbose)
end;

"""
    edgetracking(pds::ParallelDynamicalSystem, mapper::AttractorMapper; kwargs...)
Low-level function for the [`edgetracking`](@ref) algorithm.
"""
function edgetracking(pds::ParallelDynamicalSystem, mapper::AttractorMapper;
    eps1=1e-9,
    eps2=1e-8,
    maxiter=100,
    abstol=0.0,
    dt=0.01,
    tmax=Inf,
    output_level=2,
    verbose=false)
    
    verbose && println("=== Starting edge tracking algorithm ===")
    
    u1, u2 = bisect_to_edge(pds, mapper; eps=eps1)
    edgestate = (u1 + u2)/2
    
    if output_level>0
        bisect_pts1 = [u1]
        bisect_pts2 = [u2]
        edge = [edgestate]
    end
    verbose && println("... Iteration 1: Edge at $(edgestate)")

    correction = abstol + 1
    counter = 1
    # edge track iteration loop
    while (correction > abstol) & (maxiter > counter)
        reinit!(pds, [u1,u2])
        state = edgestate
        distance = dist(current_state(pds, 1) - current_state(pds, 2))
        T = 0
        # forward integration loop
        while (distance < eps2) && (T < tmax)
            step!(pds, dt)
            distance = dist(current_state(pds, 1) - current_state(pds, 2))
            T += dt
            if output_level>1
                push!(edge, (current_state(pds, 1) + current_state(pds, 2))/2)
            end
        end
        u1, u2 = bisect_to_edge(pds, mapper; eps=eps1)
        edgestate = (u1 + u2)/2
        correction = dist(edgestate - state)
        counter += 1

        if output_level>0
            push!(bisect_pts1, u1)
            push!(bisect_pts2, u2)
            push!(edge, edgestate)
        end
        
        verbose && println("... Iteration $(counter): Edge at $(edgestate)")
        (counter == maxiter) && @warn("Reached maximum number of $(maxiter) iterations;
        did not converge.")
    end

    (counter < maxiter) && println("Edge-tracking converged after $(counter) iterations.")

    (output_level == 0) && return edgestate
    (output_level == 1) && return edge
    (output_level == 2) && return [edge, bisect_pts1, bisect_pts2]
end;

"""
    bisect_to_edge(pds::ParallelDynamicalSystem, mapper::AttractorMapper; eps=1e-9)
Finds the basin boundary between two states `u1, u2 = current_states(pds)` by bisecting along a
straight line in phase space. The states `u1` and `u2` must belong to different basins.
Returns two new states located on either side of the basin boundary at a maximum 
(Euclidean) distance of `eps` between each other.

`pds` is a `ParallelDynamicalSystem` with two states. The `mapper` must be an `AttractorMapper`
of subtype `AttractorsViaProximity` or `AttractorsViaRecurrences`.

Note: If the straight line between `u1` and `u2` intersects the basin boundary multiple
times, the method will find one of these intersection points. If more than two attractors
exist, one of the two returned states may belong to a different basin than the initial
conditions `u1` and `u2`. A warning is raised if the bisection involves a third basin.

# Keyword arguments
* `eps = 1e-9`: The maximum (Euclidean) distance between the two returned states.
"""
function bisect_to_edge(pds::ParallelDynamicalSystem, mapper::AttractorMapper; eps=1e-9)
    u1, u2 = current_states(pds)
    idx1, idx2 = mapper(u1), mapper(u2)
    
    if idx1 == idx2
        error("Both initial conditions belong to the same basin of attraction.")
    end
    
    distance = dist(u1-u2)
    while distance > eps
        u_new = (u1 + u2)/2
        idx_new = mapper(u_new)
        # slightly shift u_new if it lands directly on the boundary
        if idx_new == -1
            u_new += eps*(u_1 - u_2)
        end
        # update u1 or u2
        if idx_new == idx1
            u1 = u_new
        else 
            u2 = u_new
            if idx_new != idx2
                @warn "New bisection point belongs to a third basin of attraction."
            end
        end    
        distance = dist(u1-u2)
    end
    [u1, u2]
end;

dist(x) = sqrt(sum(abs, x.^2));