"""
    AttractorsViaProximity(ds::DynamicalSystem, attractors::Dict; kwargs...)

Map initial conditions to attractors based on whether the trajectory reaches `ε`-distance
close to any of the user-provided `attractors`, which have to be in a form of a dictionary
mapping attractor labels to `StateSpaceSet`s containing the attractors.

## Keywords
* `ε = nothing`: Distance below which a trajectory has converged to an attractor, see below.
  Type `\\varepsilon<TAB>` to input `ε`.
* `Ttr = 0`: Transient time to first evolve the system for before checking for proximity.
* `Δt = 1`: Step time given to `step!`.
* `stop_at_Δt = false`: Third argument given to `step!`.
* `horizon_limit = 1e3`: If the maximum distance of the trajectory from any of the given
  attractors exceeds this limit, it is assumed
  that the trajectory diverged (gets labelled as `-1`).
* `consecutive_lost_steps = 10000`: If the `ds` has been stepped this many times without
  coming `ε`-near to any attractor,  it is assumed
  that the trajectory diverged (gets labelled as `-1`).
* `distance = StrictlyMinimumDistance()`: Distance function for evaluating the distance
  between the trajectory end-point and the given attractors. Can be anything given to
  [`set_distance`](@ref).

## Description

The system gets stepped, and at each step the distance of the current state to all
attractors is computed via `set_distance` using the `distance` keyword.
If any of these distances is `< ε`, then the label of the nearest
attractor is returned.

`attractors` do not have to be "true" attractors. Any arbitrary sets
in the state space can be provided.

If an `ε::Real` is not provided by the user, a value is computed
automatically as 1/10th of the minimum distance between all `attractors`.
This operation can be expensive for large `StateSpaceSet`s.
If `length(attractors) == 1`, then `ε` becomes 1/10 of the diagonal of the box containing
the attractor. If `length(attractors) == 1` and the attractor is a single point,
an error is thrown.

The [`convergence_time`](@ref) is `Inf` if an initial condition has not converged.
As such, the convergence time is always a float type even for discrete time systems.
"""
struct AttractorsViaProximity{DS<:DynamicalSystem, AK, SSS<:AbstractStateSpaceSet, N, K, M, SS<:AbstractStateSpaceSet, T} <: AttractorMapper
    ds::DS
    attractors::Dict{AK, SSS}
    ε::Float64
    Δt::N
    Ttr::N
    consecutive_lost_steps::Int
    horizon_limit::Float64
    search_trees::K
    dist::Vector{Float64}
    idx::Vector{Int}
    maxdist::Float64
    distance::M
    cset::SS # current state of dynamical system as `StateSpaceSet`.
    stop_at_Δt::Bool
    latest_convergence_time::Base.RefValue{T}
end

AttractorsViaProximity(ds::DynamicalSystem, attractors::Dict, ε; kw...) =
AttractorsViaProximity(ds, attractors; ε = ε, kw...)

function AttractorsViaProximity(ds::DynamicalSystem, attractors::Dict;
        Δt=1, Ttr=0, consecutive_lost_steps=10000, horizon_limit=1e3, verbose = false,
        distance = StrictlyMinimumDistance(), stop_at_Δt = false, ε = nothing,
    )
    if !(valtype(attractors) <: AbstractStateSpaceSet)
        error("The input attractors must be a dictionary with values of `StateSpaceSet`s.")
    end
    if dimension(ds) ≠ dimension(first(attractors)[2])
        error("Dimension of the dynamical system and candidate attractors must match.")
    end

    if distance isa Union{Hausdorff, StrictlyMinimumDistance}
        # We pre-initialize `KDTree`s for performance
        search_trees = Dict(k => KDTree(vec(att), distance.metric) for (k, att) in attractors)
    else
        search_trees = Dict(k => nothing for (k, att) in attractors)
    end

    if isnothing(ε)
        ε = _deduce_ε_from_attractors(attractors, search_trees, verbose)
    else
        ε isa Real || error("ε must be a Real number")
    end

    mapper = AttractorsViaProximity(
        ds, attractors,
        ε, Δt, eltype(Δt)(Ttr), consecutive_lost_steps, horizon_limit,
        search_trees, [Inf], [0], 0.0, distance, StateSpaceSet([current_state(ds)]),
        stop_at_Δt, Ref(float(current_time(ds))), # we make it float so that it can handle `Inf`.
    )

    return mapper
end

reset_mapper!(::AttractorsViaProximity) = nothing

function _deduce_ε_from_attractors(attractors, search_trees, verbose = false)
    if length(attractors) != 1
        verbose && @info("Computing minimum distance between attractors to deduce `ε`...")
        # Minimum distance between attractors
        # notice that we do not use `StateSpaceSet_distance` because
        # we have more than two StateSpaceSets and want the absolute minimum distance
        # between one of them.
        dist, idx = [Inf], [0]
        minε = Inf
        for (k, A) in attractors
            for (m, tree) in search_trees
                k == m && continue
                for p in A # iterate over all points of attractor
                    Neighborhood.NearestNeighbors.knn_point!(
                        tree, p, false, dist, idx, Neighborhood.NearestNeighbors.always_false
                    )
                    dist[1] < minε && (minε = dist[1])
                end
            end
        end
        verbose && @info("Minimum distance between attractors computed: $(minε)")
        ε = minε/10
    else
        attractor = first(attractors)[2] # get the single attractor
        mini, maxi = minmaxima(attractor)
        ε = sqrt(sum(abs, maxi .- mini))/10
        if ε == 0
            throw(ArgumentError("""
            Computed `ε = 0` in automatic estimation, probably because there is
            a single attractor that also is a single point. Please provide `ε` manually.
            """))
        end
    end
    return ε
end


# TODO: Implement `show_progress`
function (mapper::AttractorsViaProximity)(u0; show_progress = false)
    ds = referenced_dynamical_system(mapper)
    reinit!(ds, u0)
    t0 = current_time(ds)
    maxdist = zero(eltype(current_state(ds)))
    mapper.latest_convergence_time[] = Inf # default return value
    mapper.Ttr > 0 && step!(ds, mapper.Ttr)
    lost_count = 0
    while lost_count < mapper.consecutive_lost_steps
        step!(ds, mapper.Δt, mapper.stop_at_Δt)
        successful_step(ds) || return -1 # first check for Inf or NaN
        lost_count += 1
        u = current_state(ds)
        # then update the stored set
        mapper.cset[1] = u
        # then compute all distances
        for (k, tree) in mapper.search_trees # this is a `Dict`
            A = mapper.attractors[k]
            # we use internal method from StateSpaceSets.jl
            d = set_distance(mapper.cset, A, mapper.distance; tree2 = tree)
            if d < mapper.ε
                mapper.latest_convergence_time[] = current_time(ds) - t0
                return k
            elseif maxdist < d
                maxdist = d
                # exit if the distance is too large
                maxdist > mapper.horizon_limit && return -1
            end
        end
    end
    return -1
end

function Base.show(io::IO, mapper::AttractorsViaProximity)
    ps = generic_mapper_print(io, mapper)
    println(io, rpad(" ε: ", ps), mapper.ε)
    println(io, rpad(" Δt: ", ps), mapper.Δt)
    println(io, rpad(" Ttr: ", ps), mapper.Ttr)
    attstrings = split(sprint(show, MIME"text/plain"(), mapper.attractors), '\n')
    println(io, rpad(" attractors: ", ps), attstrings[1])
    for j in 2:length(attstrings)
        println(io, rpad(" ", ps), attstrings[j])
    end
    return
end

extract_attractors(mapper::AttractorsViaProximity) = mapper.attractors

convergence_time(mapper::AttractorsViaProximity) = mapper.latest_convergence_time[]
