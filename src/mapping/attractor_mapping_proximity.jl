"""
    AttractorsViaProximity(ds::DynamicalSystem, attractors::Dict [, ε]; kwargs...)

Map initial conditions to attractors based on whether the trajectory reaches `ε`-distance
close to any of the user-provided `attractors`. They have to be in a form of a dictionary
mapping attractor labels to `StateSpaceSet`s containing the attractors.

The system gets stepped, and at each step the minimum distance to all
attractors is computed. If any of these distances is `< ε`, then the label of the nearest
attractor is returned.

If an `ε::Real` is not provided by the user, a value is computed
automatically as half of the minimum distance between all attractors.
This operation can be expensive for large `StateSpaceSet`s.
If `length(attractors) == 1`, then `ε` becomes 1/10 of the diagonal of the box containing
the attractor. If `length(attractors) == 1` and the attractor is a single point,
an error is thrown.

## Keywords

* `Ttr = 100`: Transient time to first evolve the system for before checking for proximity.
* `Δt = 1`: Step time given to `step!`.
* `horizon_limit = 1e3`: If the maximum distance of the trajectory from any of the given
  attractors exceeds this limit, it is assumed
  that the trajectory diverged (gets labelled as `-1`).
* `consecutive_lost_steps = 1000`: If the integrator has been stepped this many times without
  coming `ε`-near to any attractor,  it is assumed
  that the trajectory diverged (gets labelled as `-1`).
"""
struct AttractorsViaProximity{DS<:DynamicalSystem, AK, D, T, N, K} <: AttractorMapper
    ds::DS
    attractors::Dict{AK, StateSpaceSet{D, T}}
    ε::Float64
    Δt::N
    Ttr::N
    consecutive_lost_steps::Int
    horizon_limit::Float64
    search_trees::K
    dist::Vector{Float64}
    idx::Vector{Int}
    maxdist::Float64
end
function AttractorsViaProximity(ds::DynamicalSystem, attractors::Dict, ε = nothing;
        Δt=1, Ttr=100, consecutive_lost_steps=1000, horizon_limit=1e3, verbose = false
    )
    dimension(ds) == dimension(first(attractors)[2]) ||
            error("Dimension of the dynamical system and candidate attractors must match")
    search_trees = Dict(k => KDTree(att.data, Euclidean()) for (k, att) in attractors)

    if isnothing(ε)
        ε = _deduce_ε_from_attractors(attractors, search_trees, verbose)
    else
        ε isa Real || error("ε must be a Real number")
    end

    mapper = AttractorsViaProximity(
        ds, attractors,
        ε, Δt, eltype(Δt)(Ttr), consecutive_lost_steps, horizon_limit,
        search_trees, [Inf], [0], 0.0,
    )

    return mapper
end

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
        ε = minε/2
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
    reinit!(mapper.ds, u0)
    maxdist = 0.0
    mapper.Ttr > 0 && step!(mapper.ds, mapper.Ttr)
    lost_count = 0
    while lost_count < mapper.consecutive_lost_steps
        step!(mapper.ds, mapper.Δt)
        lost_count += 1
        u = current_state(mapper.ds)
        # first check for Inf or NaN
        any(x -> (isnan(x) || isinf(x)), u) && return -1
        for (k, tree) in mapper.search_trees # this is a `Dict`
            Neighborhood.NearestNeighbors.knn_point!(
                tree, u, false, mapper.dist, mapper.idx, Neighborhood.NearestNeighbors.always_false
            )
            if mapper.dist[1] < mapper.ε
                return k
            elseif maxdist < mapper.dist[1]
                maxdist = mapper.dist[1]
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

convergence_time(mapper::AttractorsViaProximity) = current_time(mapper.ds) - initial_time(mapper.ds)
