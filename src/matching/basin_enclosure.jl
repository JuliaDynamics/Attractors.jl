"""
    MatchByBasinEnclosure(; kw...) <: IDMatcher

A matcher that matches attractors by whether they are enclosed in
the basin of a new attractor or not.

## Keyword arguments

- `ε = nothing`: distance threshold given to [`AttractorsViaProximity`](@ref).
  If `nothing`, it is estimated as half the minimum distance of centroids
  (in contrast to the default more accurate estimation in [`AttractorsViaProximity`](@ref)).
- `Δt = 1, consecutive_lost_steps = 1000`: also given to [`AttractorsViaProximity`](@ref).
  Note that attractors that did not converge to one of the current attractors
  within this number of steps
  do not get assigned ID -1 as in [`AttractorsViaProximity`](@ref). Rather, they
  get assigned the next available free ID.
- `distance = Centroid()`: metric to estimate distances between state space sets
  in case there are co-flowing attractors, see below.
- `seeding = A -> A[end]`: how to select a point from the attractor to see if
  it is enclosed in the basin of a new attractor.

## Description

An attractor `A₋` is a set in a state space that occupies a particular region
(or, a single point, if it is a fixed point).
This region is always within the basin of said attractor.
When the parameter of the dynamical system is incremented,
the attractors `A₊` in the new parameter have basins that may have changed in shape and size.

The previous attractor `A₋` is "matched" (i.e., has its ID changed)
to a new attractor `A₊` attractor if `A₋` is located inside the basin of attraction of `A₊`.
To see if `A₋` is in the basin of `A₊`, we first pick a point from `A₊` using the `seeding`
keyword argument. By default this is the last point on the attractor, but it could be anything
else, including the centroid of the attractor (`mean(A)`).
This point is given as an initial condition to an [`AttractorsViaProximity`](@ref) mapper
that maps initial conditions to the `₊` attractors when they are `ε`-close to them.

There can be the situation where multiple `₋` attractors get matched to the same `₊`
attractor, which we call "coflowing attractors". In this scenario matching is prioritized
for the `₋` attractor that is closest to the `₊` in terms of state space set distance,
which is estimated with the `distance` keyword, which can be anything
[`setsofsets_distances`](@ref) accepts. The closest `₋` coflowing attractor
gets assigned the same ID as the `₊` one, while the rest get different unique IDs.

Basin enclosure is a concept similar to "basin instability" in [Ritchie2023](@cite).
"""
@kwdef struct BasinEnclosure{E, D} <: IDMatcher
    ε::E = nothing
    distance::D = Centroid()
    seeding::S = A -> A[end]
    Δt::Float64 = 1
    consecutive_lost_steps::Int = 1000
end

function _match_attractors(
        current_attractors, prev_attractors, matcher::MatchByBasinEnclosure,
        mapper, p, pprev
    )
    ds = referenced_dynamical_system(mapper)
    if matcher.ε === nothing
        e = ε_from_centroids(attractors)
    else
        e = matcher.ε
    end
    proximity = AttractorsViaProximity(ds, current_attractors, e;
        horizon_limit = Inf, Ttr = 0, consecutive_lost_steps = matcher.consecutive_lost_steps
    )
    rmap = Dict(k => proximity(matcher.seeding(A)) for (k, A) in prev_attractors)
    # we now process the replacement map `rmap` for co-flowing or diverged attractors.
    next_id = next_free_id(current_attractors, prev_attractors)
    # take care of diverged attractors
    for (old_ID, new_ID) in rmap
        if new_ID < 0 # diverged attractors get -1 ID.
            rmap[old_ID] = next_id
            next_id += 1
        end
    end
    # next, take care of co-flowing attractors
    if unique(values(rmap)) != length(rmap) # a value is repeated
        # Do coflowing and assign `next_id` to the least distant coflowing
        error("Logic for co-flowing attractors is not implemented yet.")
        # First, find all old_IDs that are mapped to the same new_ID
        # For these estimate the distance of corresponding attractors.
        # All attractors beyond the first get assigned new ID.
    end
    return rmap
end

function ε_from_centroids(attractors)
    distances = setsofsets_distances(attractors, attractors, Centroid())
    alldists = sort!(vcat([collect(values(d)) for (k,d) in distances]...))
    filter!(!iszero, alldists)
    return minimum(alldists)/2
end
