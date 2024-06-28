export SSSetMatcher, replacement_map, replacement_map!, match_sequentially!

"""
    SSSetMatcher

Supertype of all "matchers" that match IDs between state space sets (typically attractors).

Matchers implement an extendable interface. A new mapper type only needs to
extend the function [`replacement_map`](@ref).
This function is used by the higher level function [`match_sequentially!`](@ref)
that is called at the end of [`continuation`](@ref) by most continuation algorithms.

Currently available matchers are:

- [`MatchBySSDistance`](@ref)
- [`MatchByBasinEnclosure`](@ref)

As you user you typically only care about given an instance of `SSSetMatcher`
to a continuation algorithm such as [`AttractorsContinueAndMatch`](@ref),
and you don't have to worry about the matching functions themselves.
"""
abstract type SSSetMatcher end

"""
    replacement_map(a₊, a₋, matcher; i = nothing) → rmap

Given dictionaries `a₊, a₋` mapping IDs to values,
return a _replacement map_: a dictionary mapping the IDs (keys) in dictionary `a₊`
to IDs (keys) in dictionary `a₋`, so that
so that values in `a₊` that are the "closest" to values in `a₋` get assigned the
same key as in `a₋`. In this way keys of `a₊` are "matched" to keys of `a₋`.
Use [`swap_dict_keys`](@ref) to apply `rmap` to `a₊`.

How matching happens, i.e., how "closeness" is defined, depends on the algorithm `matcher`.

The values contained in `a₊, a₋` can be anything supported by `matcher`.
Within Attractors.jl they are typically `StateSpaceSet`s representing attractors.

Typically the +,- mean after and before some change of parameter of a dynamical system.
For some matchers, such as [`MatchByBasinEnclosure`](@ref), the value of the parameter
is important. `i` in this case can be given as the index of the parameter range
corresponding to the result `a₋`, assuming `a₊, a₋` are two subsequent results
of a [`continuation`](@ref) output.
"""
function replacement_map(a₊::AbstractDict, a₋, matcher::SSSetMatcher; kw...)
    # For developers: a private keyword `next_id` is also given to `replacement_map`
    # that is utilized in the `match_sequentially!` function.
    throw(ArgumentError("Not implemented for $(typeof(matcher))"))
end

"""
    replacement_map!(a₊, a₋, matcher; i = nothing) → rmap

Convenience function that first calls [`replacement_map`](@ref) and then
replaces the IDs in `a₊` with this `rmap`.
"""
function replacement_map!(a₊::AbstractDict, a₋, matcher::SSSetMatcher; kw...)
    rmap = replacement_map(a₊, a₋, matcher; kw...)
    swap_dict_keys!(a₊, rmap)
    return rmap
end

"""
    match_sequentially!(dicts::Vector{Dict{Int, Any}}, matcher::SSSetMatcher; kw...)

Match the `dicts`, a vector of dictionaries mapping IDs (integers) to values,
according to the given `matcher` by sequentially applying the
[`replacement_map`](@ref) function to all elements of `dicts` besides the first one.

In the context of Attractors.jl `dicts` are typically dictionaries mapping
IDs to attractors (`StateSpaceSet`s), however the function is generic and would
work for any values that `matcher` works with.

Return `rmaps`, which is a vector of dictionaries.
`rmaps[i]` contains the [`replacement_map`](@ref) for `attractors[i+1]`,
i.e., the pairs of `old => new` IDs.

## Keyword arguments

- `retract_keys::Bool = true`: If `true` at the end the function will "retract" keys (i.e., make the
  integers smaller integers) so that all unique IDs
  are the 1-incremented positive integers. E.g., if the IDs where 1, 6, 8, they will become
  1, 2, 3. The special ID -1 is unaffected by this.
"""
function match_sequentially!(
        attractors::AbstractVector{<:Dict}, matcher::SSSetMatcher; retract_keys = true
    )
    # this generic implementation works for any matcher!!!
    # the matcher also provides the `use_vanished` keyword if it makes sense!
    rmaps = Dict{Int,Int}[]
    if !use_vanished(matcher) # matchers implement this!
        rmaps = _rematch_ignored!(attractors, matcher)
    else
        rmaps = _rematch_with_past!(attractors, matcher)
    end
    if retract_keys
        retracted = retract_keys_to_consecutive(attractors) # already matched input
        for (rmap, attrs) in zip(rmaps, attractors)
            swap_dict_keys!(attrs, retracted)
            # for `rmap` the situation is more tricky, because we have to change the
            # value of the _values_ of the dictionary, not the keys!
            for (k, v) in rmap
                if v ∈ keys(retracted)
                    # so we make that the replacement map points to the
                    # retracted key instead of whatever it pointed to originally,
                    # if this key exists in the retracted mapping
                    rmap[k] = retracted[v]
                end
            end
        end
        # `attractors` have 1 more element than `rmaps`
        swap_dict_keys!(attractors[end], retracted)
    end
    return rmaps
end

"""
    match_sequentially!(continuation_quantity::Vector{Dict}, rmaps::Vector{Dict})

Do the same as in `match_sequentially!` above, now given the vector of replacement maps,
and for any arbitrary quantity that has been tracked in the continuation.
`continuation_quantity` can for example be `fractions_curves` from [`continuation`](@ref).
"""
function match_sequentially!(continuation_quantity::AbstractVector{<:Dict}, rmaps::Vector{Dict{Int, Int}})
    if length(rmaps) ≠ length(continuation_quantity) - 1
        throw(ArgumentError("the replacement maps should be 1 less than the continuation quantities"))
    end
    for (i, rmap) in enumerate(rmaps)
        quantity = continuation_quantity[i+1]
        swap_dict_keys!(quantity, rmap)
    end
    return rmaps
end

# Concrete implementation of `match_sequentially!`:
function _rematch_ignored!(attractors_info, matcher)
    next_id = 1
    rmaps = Dict{keytype(attractors_info[1]), keytype(attractors_info[1])}[]
    for i in 1:length(attractors_info)-1
        a₊, a₋ = attractors_info[i+1], attractors_info[i]
        # If there are no attractors, skip the matching
        (isempty(a₊) || isempty(a₋)) && continue
        # Here we always compute a next id. In this way, if an attractor disappears
        # and reappears, it will get a different (incremented) ID as it should!
        next_id_a = max(maximum(keys(a₊)), maximum(keys(a₋)))
        next_id = max(next_id, next_id_a) + 1
        rmap = replacement_map!(a₊, a₋, matcher; next_id, i)
        push!(rmaps, rmap)
    end
    return rmaps
end

function _rematch_with_past!(attractors_info, matcher)
    # this dictionary stores all instances of previous attractors and is updated
    # at every step. It is then given to the matching function as if it was
    # the current attractors
    latest_ghosts = deepcopy(attractors_info[1])
    rmaps = Dict{keytype(attractors_info[1]), keytype(attractors_info[1])}[]
    for i in 1:length(attractors_info)-1
        a₊, a₋ = attractors_info[i+1], attractors_info[i]
        # update ghosts
        for (k, A) in a₋
            latest_ghosts[k] = A
        end
        rmap = replacement_map!(a₊, latest_ghosts, matcher; i)
        push!(rmaps, rmap)
    end
    return rmaps
end

include("match_basins.jl")
include("sssdistance.jl")
include("basin_enclosure.jl")