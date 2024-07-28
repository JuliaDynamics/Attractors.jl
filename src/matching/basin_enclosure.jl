export MatchByBasinEnclosure

"""
    MatchByBasinEnclosure(; kw...) <: IDMatcher

A matcher that matches attractors by whether they are enclosed in
the basin of a new attractor or not.

## Keyword arguments

- `ε = nothing`: distance threshold given to [`AttractorsViaProximity`](@ref).
  If `nothing`, it is estimated as a quarter of the minimum distance of centroids
  (in contrast to the default more accurate estimation in [`AttractorsViaProximity`](@ref)).
- `Δt = 1, consecutive_lost_steps = 1000`: also given to [`AttractorsViaProximity`](@ref).
  We have not yet decided what should happen to attractors that did not converge to one of the current attractors
  within this number of steps. At the moment they
  get assigned the next available free ID but this may change in future releases.
- `distance = Centroid()`: metric to estimate distances between state space sets
  in case there are co-flowing attractors, see below.
- `seeding = A -> A[end]`: how to select a point from the attractor to see if
  it is enclosed in the basin of a new attractor.

## Description

An attractor `A₋` is a set in a state space that occupies a particular region
(or, a single point, if it is a fixed point).
This region is always within the basin of attraction of said attractor.
When the parameter of the dynamical system is incremented,
the attractors `A₊` in the new parameter have basins that may have changed in shape and size.

The new attractor `A₊` is "matched" (i.e., has its ID changed)
to the old attractor `A₋` attractor if `A₋` is located inside the basin of attraction of `A₊`.
To see if `A₋` is in the basin of `A₊`, we first pick a point from `A₋` using the `seeding`
keyword argument. By default this is the last point on the attractor, but it could be anything
else, including the centroid of the attractor (`mean(A)`).
This point is given as an initial condition to an [`AttractorsViaProximity`](@ref) mapper
that maps initial conditions to the `₊` attractors when
the trajectories from the initial conditions are `ε`-close to the `₊` attractors.

There can be the situation where multiple `₋` attractors converge to the same `₊`
attractor, which we call "coflowing attractors". In this scenario matching is prioritized
for the `₋` attractor that is closest to the `₊` in terms of state space set distance,
which is estimated with the `distance` keyword, which can be anything
[`MatchBySSSetDistance`](@ref) accepts. The closest `₊` attractor gets the
ID of the `₋` closest attractor that converge to it.

Basin enclosure is a concept similar to "basin instability" in [Ritchie2023](@cite).
"""
@kwdef struct MatchByBasinEnclosure{E, D, S, T} <: IDMatcher
    ε::E = nothing
    distance::D = Centroid()
    seeding::S = A -> A[end]
    Δt::T = 1
    consecutive_lost_steps::Int = 1000
end

function matching_map(
        current_attractors, prev_attractors, matcher::MatchByBasinEnclosure;
        ds, p, pprev = nothing, next_id = next_free_id(current_attractors, prev_attractors)
    )
    if matcher.ε === nothing
        e = ε_from_centroids(current_attractors)
    else
        e = matcher.ε
    end
    set_parameters!(ds, p)
    proximity = AttractorsViaProximity(ds, current_attractors, e;
        horizon_limit = Inf, Ttr = 0, consecutive_lost_steps = matcher.consecutive_lost_steps
    )
    # we start building the "flow" map mapping previous attractors
    # to where they flowed to in current attractors
    # (notice that `proximity(u)` returns IDs of current attractors)
    flow = Dict(k => proximity(matcher.seeding(A)) for (k, A) in prev_attractors)
    # of course, the matching map is the inverse of `flow`
    rmap = Dict{Int, Int}()
    # but we need to take care of diverging and co-flowing attractors.
    # Let's start with diverging ones
    for (old_ID, new_ID) in flow
        if new_ID < 0 # diverged attractors get -1 ID
            # but here we assign them just to the next available integer
            rmap[new_ID] = next_id
            next_id += 1
            delete!(flow, old_ID) # and remove from flow map
        end
    end
    # next up are the co-flowing attractors
    grouped_flows = _grouped_flows(flow)
    # notice the keys of `grouped_flows` are new IDs, same as with `rmap`.
    for (new_ID, old_flowed_to_same) in grouped_flows
        if length(old_flowed_to_same) == 0
            continue # none of the old IDs converged to the current `new_ID`
        elseif length(old_flowed_to_same) == 1
            rmap[new_ID] = only(old_flowed_to_same)
        else # need to resolve coflowing using distances
            a₊ = Dict(new_ID => current_attractors[new_ID])
            a₋ = Dict(old_ID => prev_attractors[old_ID] for old_ID in old_flowed_to_same)
            ssmatcher = MatchBySSSetDistance(; distance = matcher.distance)
            @show length(a₊), length(a₋)
            matched_rmap = matching_map(a₊, a₋, ssmatcher)
            # this matcher has only one entry, so we use it to match
            # (we don't care what happens to the rest of the old_IDs, as the `rmap`
            # only cares about what adjustments need to happen to the new_IDs)
            new_ID, old_ID = only(matched_rmap)# our main `rmap`
            rmap[new_ID] = old_ID
        end
    end
    return rmap
end

function ε_from_centroids(attractors::AbstractDict)
    if length(attractors) == 1 # `attractors` has only 1 attractor
        attractor = first(attractors)[2] # get the single attractor
        mini, maxi = minmaxima(attractor)
        ε = sqrt(sum(abs, maxi .- mini))/10
        if ε == 0
            throw(ArgumentError("""
            Computed `ε = 0` in automatic estimation for `AttractorsViaFeaturizing`, probably because there is
            only a single attractor that also is a single point. Please provide `ε` manually.
            """))
        end
        return ε
    end
    # otherwise compute cross-distances
    distances = setsofsets_distances(attractors, attractors, Centroid())
    alldists = sort!(vcat([collect(values(d)) for (k,d) in distances]...))
    filter!(!iszero, alldists)
    return minimum(alldists)/4
end

# group flows so that all old IDs that go to same new ID are in one vector
function _grouped_flows(flows) # separated into
    grouped = Dict{Int, Vector{Int}}()
    oldids = collect(keys(flows))
    for k in values(flows)
        grouped[k] = findall(isequal(k), oldids)
    end
    return grouped
end