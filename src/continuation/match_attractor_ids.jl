# Notice this file uses heavily `dict_utils.jl`!
export match_statespacesets!, match_basins_ids!, replacement_map, match_continuation!

###########################################################################################
# Matching attractors and key swapping business
###########################################################################################
"""
    match_statespacesets!(a₊::AbstractDict, a₋; distance = Centroid(), threshold = Inf)

Given dictionaries `a₊, a₋` mapping IDs to `StateSpaceSet` instances,
match the IDs in dictionary `a₊` so that its sets that are the closest to
those in dictionary `a₋` get assigned the same key as in `a₋`.
Typically the +,- mean after and before some change of parameter of a system.

Return the replacement map, a dictionary mapping old keys of `a₊` to
the new ones that they were mapped to. You can obtain this map, without modifying
the dictionaries, by calling the [`replacement_map`](@ref) function directly.

## Keyword arguments

- `distance = Centroid()`: given to [`setsofsets_distances`](@ref).
- `threshold = Inf`: attractors with distance larger than the `threshold` are guaranteed
  to not be mapped to each other.

## Description

When finding attractors and their fractions in Attractors.jl,
different attractors get assigned different IDs. However
which attractor gets which ID is somewhat arbitrary. Finding the attractors of the
same system for slightly different parameters could label "similar" attractors (at
the different parameters) with different IDs.
`match_statespacesets!` tries to "match" them by modifying the IDs,
i.e., the keys of the given dictionaries. Do note however that there is nothing
in this function that is limited to attractors in the formal mathematical sense.
Any dictionary with `StateSpaceSet` values is a valid input and these sets
may represent attractors, trajectories, group of features, or anything else.

The matching happens according to the output of the [`setsofsets_distances`](@ref)
function with the keyword `distance`. `distance` can be whatever that function accepts,
i.e., one of `Centroid, Hausdorff, StrictlyMinimumDistance` or any arbitrary user-
provided function that given two sets it returns a positive number (their distance).
State space sets are then matched according to this distance.
First, all possible pairs (old, new, distance) are sorted according to their distance.
The pair with smallest distance is matched. Sets in matched pairs are removed from the
matching pool to ensure a unique mapping. Then, the next pair with least
remaining distance is matched, and the process repeats until all pairs are exhausted.

Additionally, you can provide a `threshold` value. If the distance between two attractors
is larger than this `threshold`, then it is guaranteed that the attractors will get assigned
different key in the dictionary `a₊` (which is the next available integer).
"""
function match_statespacesets!(a₊::AbstractDict, a₋; kwargs...)
    rmap = replacement_map(a₊, a₋; kwargs...)
    swap_dict_keys!(a₊, rmap)
    return rmap
end

# Note that `next_id` is an internal argument not exposed to Public API.
# This is used when we ignore previously existing attractors in
# `rematch_continuation!`
"""
    replacement_map(a₊, a₋; distance = Centroid(), threshold = Inf) → rmap
Return a dictionary mapping keys in `a₊` to new keys in `a₋`,
as explained in [`match_statespacesets!`](@ref).
"""
function replacement_map(a₊::Dict, a₋::Dict;
        distance = Centroid(), threshold = Inf, next_id = nothing
    )
    distances = setsofsets_distances(a₊, a₋, distance)
    keys₊, keys₋ = keys.((a₊, a₋))
    nextid = isnothing(next_id) ? max(maximum(keys₊), maximum(keys₋)) + 1 : next_id
    replacement_map(keys₊, keys₋, distances::Dict, threshold, nextid)
end

function replacement_map(keys₊, keys₋, distances::Dict, threshold, next_id = max(maximum(keys₊), maximum(keys₋)) + 1)
    # Transform distances to sortable collection. Sorting by distance
    # ensures we prioritize the closest matches
    sorted_keys_with_distances = Tuple{Int, Int, Float64}[]
    for i in keys(distances)
        for j in keys(distances[i])
            push!(sorted_keys_with_distances, (i, j, distances[i][j]))
        end
    end
    sort!(sorted_keys_with_distances; by = x -> x[3])

    # Iterate through distances, find match with least distance, and "remove" (skip)
    # all remaining same indices a'la Eratosthenis sieve
    # In the same loop we match keys according to distance of values,
    # but also ensure that keys that have too high of a value distance are guaranteeed
    # to have different keys, and ensure that there is unique mapping happening!
    rmap = Dict{eltype(keys₊), eltype(keys₋)}()
    done_keys₊ = eltype(keys₊)[] # stores keys of a₊ already processed
    used_keys₋ = eltype(keys₋)[] # stores keys of a₋ already used
    for (oldkey, newkey, dist) in sorted_keys_with_distances
        # used keys can't be re-used to match again
        (oldkey ∈ done_keys₊ || newkey ∈ used_keys₋) && continue
        if dist < threshold
            push!(used_keys₋, newkey)
        else
            # The distance exceeds threshold, so we will assign a new key
            # (notice that this assumes the sorting by distance we did above,
            # otherwise it wouldn't work!)
            newkey = next_id
            next_id += 1
        end
        push!(done_keys₊, oldkey)
        rmap[oldkey] = newkey
    end

    # if not all keys were processed, we map them to the next available integers
    if length(done_keys₊) ≠ length(keys₊)
        unprocessed = setdiff(collect(keys₊), done_keys₊)
        for oldkey in unprocessed
            rmap[oldkey] = next_id
            next_id += 1
        end
    end
    return rmap
end



###########################################################################################
# Matching with basins and possibly overlaps
###########################################################################################
"""
    match_basins_ids!(b₊::AbstractArray, b₋; threshold = Inf)
Similar to [`match_statespacesets!`](@ref) but operate on basin arrays instead
(the arrays typically returned by [`basins_of_attraction`](@ref)).

This method matches IDs of attractors whose basins of attraction before and after `b₋,b₊`
have the most overlap (in pixels). This overlap is normalized in 0-1 (with 1 meaning
100% overlap of pixels). The `threshold` in this case is compared to the inverse
of the overlap (so, for `threshold = 2` attractors that have less than 50% overlap get
different IDs guaranteed).
"""
function match_basins_ids!(b₊::AbstractArray, b₋; threshold = Inf)
    ids₊, ids₋ = unique(b₊), unique(b₋)
    distances = _similarity_from_overlaps(b₊, ids₊, b₋, ids₋)
    rmap = replacement_map(ids₊, ids₋, distances, threshold)
    replace!(b₊, rmap...)
    return rmap
end

function _similarity_from_overlaps(b₊, ids₊, b₋, ids₋)
    size(b₋) == size(b₊) || error("Sizes of b₊ and  b₋ must match")
    distances = Dict{eltype(ids₊), Dict{eltype(ids₋), Float64}}()
    for i in ids₊
        Bi = findall(isequal(i), b₊)
        d = valtype(distances)()
        # Compute normalized overlaps of each basin with each other basis
        for j in ids₋
            Bj = findall(isequal(j), b₋)
            overlap = length(Bi ∩ Bj)/length(Bj)
            d[j] = 1 / overlap # distance is inverse overlap
        end
        distances[i] = d
    end
    return distances
end


###########################################################################################
# Continuation matching (which can be used after continuation)
###########################################################################################
"""
    match_continuation!(fractions_curves::Vector{<:Dict}, attractors_info::Vector{<:Dict}; kwargs...)

Loop over all entries in the given arguments (which are typically the direct outputs of
[`continuation`](@ref) with [`RecurrencesFindAndMatch`](@ref)), and match the
attractor IDs in both the attractors container and the basins fractions container.
This means that we loop over each entry of the vectors (besides the first),
and in each entry we attempt to match the dictionary keys to the keys of the
previous dictionary using [`match_statespacesets`](@ref).

The keywords `distance, threshold` are propagated to [`match_statespacesets`](@ref).
However, there is a unique keyword for `match_continuation!`: `use_vanished::Bool`.
If `true`, then attractors that existed before but have vanished are kept in "memory"
when it comes to matching: the new attractors are compared to the latest istance
of all attractors that have ever existed, and get match to their closest ones
as per [`match_statespacesets!`](@ref).
If `false`, vanished attractors are ignored. Note that in this case new attractors
that cannot be matched to any previous attractors will get an appropriately
incremented ID. E.g., if we started with three attractors, and attractor 3 vanished,
and at some later parameter value we again have three attractors, the new third
attractor will _not_ have ID 3, but 4 (i.e., the next available ID).

By default `use_vanished = !isinf(threshold)` and since the default value for
`threshold` is `Inf`, `use_vanished` is `false`.


    rematch_continuation!(attractors_info::Vector{<:Dict}; kwargs...)

This is a convenience method that only uses and modifies the state space set dictionary
container without the need for a basins fractions container.
"""
function match_continuation!(attractors_info; kwargs...)
    fractions_curves = [Dict(k => nothing for k in keys(d)) for d in attractors_info]
    match_continuation!(fractions_curves, attractors_info; kwargs...)
end
function match_continuation!(
        fractions_curves::Vector{<:Dict}, attractors_info::Vector{<:Dict};
        threshold = Inf, use_vanished = !isinf(threshold), kwargs...
    )
    if !use_vanished
        _rematch_ignored!(fractions_curves, attractors_info; threshold, kwargs...)
    else
        _rematch_with_past!(fractions_curves, attractors_info; threshold, kwargs...)
    end
    # This part normalizes so that keys increment by +1
    rmap = retract_keys_to_consecutive(fractions_curves)
    for (da, df) in zip(attractors_info, fractions_curves)
        swap_dict_keys!(da, rmap)
        swap_dict_keys!(df, rmap)
    end
end

function _rematch_ignored!(fractions_curves, attractors_info; kwargs...)
    next_id = 1
    for i in 1:length(attractors_info)-1
        a₊, a₋ = attractors_info[i+1], attractors_info[i]
        # Here we always compute a next id. In this way, if an attractor dissapears
        # and re-appears, it will get a different (incremented) id as it should!
        next_id_a = max(maximum(keys(a₊)), maximum(keys(a₋))) + 1
        next_id = max(next_id+1, next_id_a)
        rmap = match_statespacesets!(a₊, a₋; next_id, kwargs...)
        swap_dict_keys!(fractions_curves[i+1], rmap)
    end
end

function _rematch_with_past!(fractions_curves, attractors_info; kwargs...)
    # this dictionary stores all instances of previous attractors and is updated
    # at every step. It is then given to the matching function as if it was
    # the current attractors
    latest_ghosts = copy(attractors_info[1])
    for i in 1:length(attractors_info)-1
        a₊, a₋ = attractors_info[i+1], attractors_info[i]
        # update ghosts
        for (k, A) in a₋
            latest_ghosts[k] = A
        end
        rmap = match_statespacesets!(a₊, latest_ghosts; kwargs...)
        swap_dict_keys!(fractions_curves[i+1], rmap)
    end
end
