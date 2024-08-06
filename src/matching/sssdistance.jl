export MatchBySSSetDistance

"""
    MatchBySSSetDistance(; distance = Centroid(), threshold = Inf, use_vanished = false)

A matcher type that matches IDs by the distance of their corresponding state space sets.

## Keyword arguments

- `distance = Centroid()`: distance to match by, given to [`setsofsets_distances`](@ref).
- `threshold = Inf`: sets with distance larger than the `threshold` are guaranteed
  to not be mapped to each other.
- `use_vanished = !isinf(threshold)`: value of the keyword `use_vanished` when
  used in [`match_sequentially!`](@ref).

## Description

In this matcher the values compared are [`StateSpaceSet`](@ref)s
which in most cases represent attractors in the state space, but may also
represent any other set such as a group of features.

Here is how this matcher works:
(recall in this conversation that sets/attractors are stored in dictionaries,
mapping keys/IDs to the sets, and we want to match keys in the "new" dictionary (`a₊`)
to those in the "old" dictionary (`a₋`)).

The distance between all possible pairs of sets between the "old" sets and "new" sets
is computed as a formal distance between sets.
This is controlled by the `distance` option, itself given to the lower-level
[`setsofsets_distances`](@ref) function, so `distance` can be whatever that function accepts.
That is, one of [`Centroid`](@ref), [`Hausdorff`](@ref), [`StrictlyMinimumDistance`](@ref),
or any arbitrary user-provided function `f` that given two sets `f(A, B)` it returns a
positive number (their distance).

Sets (in particular, their corresponding IDs) are then matched according to this distance.
First, all possible ID pairs (old, new) are sorted according to the distance of their corresponding sets.
The pair with smallest distance is matched. IDs in matched pairs are removed from the
matching pool to ensure a unique mapping. Then, the next pair with least
remaining distance is matched, and the process repeats until all pairs are exhausted.

Additionally, you can provide a `threshold` value. If the distance between two sets
is larger than this `threshold`, then it is guaranteed that the two sets will get assigned
different ID in the replacement map, and hence, the set in `a₊` gets the next available
integer as its ID.
"""
@kwdef struct MatchBySSSetDistance{M, T<:Real} <: IDMatcher
    distance::M = Centroid()
    threshold::T = Inf
    use_vanished::Bool = !isinf(threshold)
end

_use_vanished(m::MatchBySSSetDistance) = m.use_vanished

function matching_map(a₊::AbstractDict, a₋::AbstractDict, matcher::MatchBySSSetDistance;
        kw...
    )
    distances = setsofsets_distances(a₊, a₋, matcher.distance)
    keys₊, keys₋ = sort.(collect.(keys.((a₊, a₋))))
    _matching_map_distances(keys₊, keys₋, distances::Dict, matcher.threshold; kw...)
end

function _matching_map_distances(keys₊, keys₋, distances::Dict, threshold;
        next_id = next_free_id(keys₊, keys₋), kw... # rest of keywords aren't used.
    )
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
    # but also ensure that keys that have too high of a value distance are guaranteed
    # to have different keys, and ensure that there is unique mapping happening!
    rmap = Dict{eltype(keys₊), eltype(keys₋)}()
    done_keys₊ = eltype(keys₊)[] # stores keys of a₊ already processed
    used_keys₋ = eltype(keys₋)[] # stores keys of a₋ already used
    for (oldkey, newkey, dist) in sorted_keys_with_distances
        # used keys can't be re-used to match again
        (oldkey ∈ done_keys₊ || newkey ∈ used_keys₋) && continue
        if dist < threshold
            # If the distance is small enough, we store the "new" key
            # in the replacement map (see after `if`), and we also
            # ensure that this new key is marked as "used"
            push!(used_keys₋, newkey)
        else
            # If the distance is too large, we assign a new key
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


