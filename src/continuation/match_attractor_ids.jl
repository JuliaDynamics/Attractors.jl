# Notice this file uses heavily `dict_utils.jl`!
export match_attractor_ids!, match_basins_ids!, replacement_map, rematch!

###########################################################################################
# Matching attractors and key swapping business
###########################################################################################
"""
    match_attractor_ids!(a₊::AbstractDict, a₋; distance = Centroid(), threshold = Inf)

Given dictionaries `a₊, a₋` mapping IDs to attractors (`StateSpaceSet` instances),
match attractor IDs in dictionary `a₊` so that its attractors that are the closest to
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

When finding attractors and their fractions in DynamicalSystems.jl,
different attractors get assigned different IDs. However
which attractor gets which ID is somewhat arbitrary. Finding the attractors of the
same system for slightly different parameters could label "similar" attractors (at
the different parameters) with different IDs.
`match_attractors_ids!` tries to "match" them by modifying the attractor IDs,
i.e., the keys of the given dictionaries.

The matching happens according to the output of the [`setsofsets_distances`](@ref)
function with the keyword `distance`. `distance` can be whatever that function accepts,
i.e., one of `Centroid, Hausdorff, StrictlyMinimumDistance` or any arbitrary user-
provided function that given two sets it returns a positive number (their distance).
Attractors are then matched according to this distance, with unique mapping.
The closest attractors (before and after) are mapped to each
other, and are removed from the matching pool, and then the next pair with least
remaining distance is matched, and so on.

Additionally, you can provide a `threshold` value. If the distance between two attractors
is larger than this `threshold`, then it is guaranteed that the attractors will get assigned
different key in the dictionary `a₊`.
"""
function match_attractor_ids!(a₊::AbstractDict, a₋; kwargs...)
    rmap = replacement_map(a₊, a₋; kwargs...)
    swap_dict_keys!(a₊, rmap)
    return rmap
end

"""
    match_attractor_ids!(attractors_info::Vector{<:Dict}; kwargs...)

Convenience method that progresses through the dictionaries in `attractors_info` in sequence
and matches them using same keywords as the above method.
"""
function match_attractor_ids!(as::Vector{<:Dict}; kwargs...)
    for i in 1:length(as)-1
        a₊, a₋ = as[i+1], as[i]
        match_attractor_ids!(a₊, a₋; kwargs...)
    end
end

"""
    replacement_map(a₊, a₋; distance = Centroid(), threshold = Inf) → rmap
Return a dictionary mapping keys in `a₊` to new keys in `a₋`,
as explained in [`match_attractor_ids!`](@ref).
"""
function replacement_map(a₊::Dict, a₋::Dict; distance = Centroid(), threshold = Inf)
    distances = setsofsets_distances(a₊, a₋, distance)
    keys₊, keys₋ = keys.((a₊, a₋))
    replacement_map(keys₊, keys₋, distances::Dict, threshold)
end

function replacement_map(keys₊, keys₋, distances::Dict, threshold)
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
    next_id = max(maximum(keys₊), maximum(keys₋)) + 1
    done_keys₊ = eltype(keys₊)[] # stores keys of a₊ already processed
    used_keys₋ = eltype(keys₋)[] # stores keys of a₋ already used
    for (oldkey, newkey, dist) in sorted_keys_with_distances
        (oldkey ∈ done_keys₊ || newkey ∈ used_keys₋) && continue
        if  dist < threshold
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
Similar to [`match_attractor_ids!`](@ref) but operate on basin arrays instead
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
# Rematch! (which can be used after continuation)
###########################################################################################
"""
    rematch!(fractions_curves, attractors_info; kwargs...)

Given the outputs of [`continuation`](@ref) witn [`RecurrencesFindAndMatch`](@ref),
perform the matching step of the process again with the (possibly different) keywords
that [`match_attractor_ids!`](@ref) accepts. This "re-matching" is possible because in
[`continuation`](@ref) finding the attractors and their basins is a completely independent
step from matching them with their IDs in the previous parameter value.
"""
function rematch!(fractions_curves, attractors_info; kwargs...)
    for i in 1:length(attractors_info)-1
        a₊, a₋ = attractors_info[i+1], attractors_info[i]
        rmap = match_attractor_ids!(a₊, a₋; kwargs...)
        swap_dict_keys!(fractions_curves[i+1], rmap)
    end
    rmap = retract_keys_to_consecutive(fractions_curves)
    for (da, df) in zip(attractors_info, fractions_curves)
        swap_dict_keys!(da, rmap)
        swap_dict_keys!(df, rmap)
    end
end
