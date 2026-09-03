# Matching is part of an extendable interface. Currently, all matchers that can be
# used in global continuation have a type parameter flag for using the vanished
# approach to matching, or the next-id approach. If in the future more approaches
# become possible, we should make a dedicated Enum type and make this the type parameter.

# For now, the only parts exposed are these functions:
export matching_map, matching_map!, match_sequentially!, IDMatcher, DontMatch
# all of which take as input the treshold and distance of the
# `MatchBySSSetDistance` matcher.

"""
    IDMatcher

Supertype of all "matchers" that match can IDs labelling attractors.
Currently available matchers:

- [`MatchBySSSetDistance`](@ref)
- [`MatchByBasinEnclosure`](@ref)
- [`MatchByBasinOverlap`](@ref)
- [`DontMatch`](@ref)

Matchers implement an extendable interface based on the function [`matching_map`](@ref).
This function is used by the higher level function [`match_sequentially!`](@ref),
which can be called after any call to a global continuation to match attractors
differently, if the matching used originally during the continuation was not optimal.
"""
abstract type IDMatcher end

"""
    DontMatch <: IDMatcher
    DontMatch()

Do no matching, i.e. return an empty dictionary when called in [`matching_map`](@ref).
"""
struct DontMatch <: IDMatcher end

"""
    matching_map(
        a₊::Dict, a₋::Dict, matcher::IDMatcher;
        ds::DynamicalSystem, p, pprev, next_id
    ) → rmap

Given dictionaries `a₊, a₋` mapping IDs to values,
return a _matching map_: a dictionary mapping the IDs (keys) in dictionary `a₊`
to IDs (keys) in dictionary `a₋`, so that
so that values in `a₊` that are the "closest" to values in `a₋` get assigned the
same key as in `a₋`. In this way keys of `a₊` are "matched" to keys of `a₋`.
Use [`swap_dict_keys!`](@ref) to apply `rmap` to `a₊`
or to other dictionaries with same keys as `a₊`.

How matching happens, i.e., how "closeness" is defined, depends on the algorithm `matcher`.

The values contained in `a₊, a₋` can be anything supported by `matcher`.
Within Attractors.jl they are typically `StateSpaceSet`s representing attractors.
Typically the +,- mean after and before some change of parameter of a dynamical system.

`matching_map` always returns an empty dictionary of either `a₊, a₋` is empty.

## Keyword arguments

- `ds`: the dynamical system that generated `a₊, a₋`.
- `p, pprev`: the parameters corresponding to `a₊, a₋`. Both need to be iterables mapping
  parameter index to parameter value (such as `Dict, Vector{Pair}`, etc., so whatever
  can be given as input to `DynamicalSystems.set_parameters!`).
- `next_id = next_free_id(a₊, a₋)`: the ID to give to values of  `a₊` that cannot be
  matched to `a₋` and hence must obtain a new unique ID.

Some matchers like [`MatchBySSSetDistance`](@ref) do not utilize `ds, p, pprev` in any way
while other matchers like [`MatchByBasinEnclosure`](@ref) do, and those require
expliticly giving values to `ds, p, pprev` as their default values
is just `nothing`.
"""
function matching_map(a₊::AbstractDict, a₋::AbstractDict, matcher::DontMatch; kw...)
    # For developers: a private keyword `next_id` can be given to `matching_map`
    # that is utilized in the `match_sequentially!` function.
    return Dict{keytype(a₊), keytype(a₋)}()
end

"""
    matching_map!(a₊, a₋, matcher; kw...) → rmap

Convenience function that first calls [`matching_map`](@ref) and then
replaces the IDs in `a₊` with this `rmap`.
"""
function matching_map!(a₊::AbstractDict, a₋, matcher::IDMatcher; kw...)
    rmap = matching_map(a₊, a₋, matcher; kw...)
    swap_dict_keys!(a₊, rmap)
    return rmap
end

"""
    match_sequentially!(dicts::Vector{Dict{Int, Any}}, matcher::IDMatcher; kw...)

Match the `dicts`, a vector of dictionaries mapping IDs (integers) to values,
according to the given `matcher` by sequentially applying the
[`matching_map`](@ref) function to all elements of `dicts` besides the first one.

In the context of Attractors.jl `dicts` are typically dictionaries mapping
IDs to attractors (`StateSpaceSet`s), however the function is generic and would
work for any values that `matcher` works with.

Return `rmaps`, which is a vector of dictionaries.
`rmaps[i]` contains the [`matching_map`](@ref) for `dicts[i+1]`,
i.e., the pairs of `old => new` IDs.

## Keyword arguments

- `pcurve = nothing`: the curve of parameters along which the continuation occured,
  from which to extract the `p, pprev` values given to [`matching_map`](@ref).
  See [`global_continuation`](@ref) if you are unsure what this means.
- `ds = nothing`: dynamical system propagated to [`matching_map`](@ref).
- `retract_keys::Bool = true`: If `true` at the end the function will "retract" keys (i.e., make the
  integers smaller integers) so that all unique IDs
  are the 1-incremented positive integers. E.g., if the IDs where 1, 6, 8, they will become
  1, 2, 3. The special ID -1 is unaffected by this.
"""
function match_sequentially!(
        attractors_cont::AbstractVector{<:Dict}, matcher::IDMatcher;
        ds = nothing, p = nothing, pprev = nothing, # parameter and ds keywords
        # TODO: Remove these keywords:
        retract_keys = _retract_keys(matcher), use_vanished = _use_vanished(matcher),
    )
    # this generic implementation works for any matcher!!!
    # the matcher also provides the `use_vanished` keyword if it makes sense!
    rmaps = Dict{keytype(attractors_cont[1]), keytype(attractors_cont[1])}[]
    tracker = init_matching_tracker(attractors_cont, matcher)
    for i in 1:(length(attractors_cont) - 1)
        a₊, a₋ = attractors_cont[i + 1], attractors_cont[i]
        p, pprev = pcurve[i + 1], pcurve[i]
        tracker = update_matching_tracker(tracker, matcher, a₊, a₋)
        rmap = tracked_matching_map!(a₊, a₋, matcher, tracker, ds, p, pprev)
        push!(rmaps, rmap)
    end
    # note this changes keys in both attractors and rmaps
    retract_keys && retract_keys!(attractors_cont, rmaps)
    return rmaps
end

# The following three helper functions allow the sequential mapping to have simple
# and elegant code that handles both sequential behaviours: vanished attractors or
# incremental IDs. Before we had two individual functions for these behaviours.
# The "tracking" variable depends on whether # we matching with ghosts or without,
# in which case we need unique next ids.
function init_matching_tracker(attractors::AbstractVector{<:Dict}, matcher::IDMatcher)
    if _use_vanished(matcher) # return ghost attractor container
        return latest_ghosts = deepcopy(attractors_cont[1])
    else
        return next_id = 1
    end
end

function update_matching_tracker(tracker, matcher, a₊, a₋)
    if _use_vanished(matcher) # then it's the `latest_ghosts`
        latest_ghosts = tracker
        for (k, A) in a₋
            latest_ghosts[k] = A
        end
        return latest_ghosts
    else # then it's the `next_id`
        next_id = tracker
        next_id_a = max(maximum(keys(a₊)), maximum(keys(a₋)))
        next_id = max(next_id, next_id_a) + 1
        return next_id
    end
end

function tracked_matching_map!(a₊, a₋, matcher, tracker, ds, p, pprev)
    if _use_vanished(matcher) # then the `tracker` is the `latest_ghosts`
        rmap = matching_map!(a₊, tracker, matcher; pprev, p, ds)
    else # then the `tracker` is the `next_id`
        rmap = matching_map!(a₊, a₋, matcher; next_id = tracker, ds, p, pprev)
    end
    return rmap
end

function retract_keys!(attractors, rmaps)
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
    return
end


"""
    match_sequentially!(continuation_quantity::Vector{Dict}, rmaps::Vector{Dict})

Do the same as in `match_sequentially!` above, now given the vector of matching maps,
and for any arbitrary quantity that has been tracked in the global continuation.
`continuation_quantity` can for example be `fractions_cont` from [`global_continuation`](@ref).
"""
function match_sequentially!(continuation_quantity::AbstractVector{<:Dict}, rmaps::Vector{Dict{Int, Int}})
    if length(rmaps) ≠ length(continuation_quantity) - 1
        throw(ArgumentError("the matching maps should be 1 less than the global_continuation quantities"))
    end
    for (i, rmap) in enumerate(rmaps)
        quantity = continuation_quantity[i + 1]
        swap_dict_keys!(quantity, rmap)
    end
    return rmaps
end

# The following two functions were prior concrete implementations of `match_sequentially!`,
# that are kept here for educational purposes!
# function _rematch_ignored!(
#         attractors_cont, matcher;
#         ds = nothing, pcurve = eachindex(attractors_cont),
#     )
#     next_id = 1
#     rmaps = Dict{keytype(attractors_cont[1]), keytype(attractors_cont[1])}[]
#     for i in 1:(length(attractors_cont) - 1)
#         a₊, a₋ = attractors_cont[i + 1], attractors_cont[i]
#         p, pprev = pcurve[i + 1], pcurve[i]
#         # If there attractors, update the max id
#         if !(isempty(a₊) || isempty(a₋))
#             # we always compute a next id. In this way, if an attractor disappears
#             # and reappears, it will get a different (incremented) ID as it should!
#             next_id_a = max(maximum(keys(a₊)), maximum(keys(a₋)))
#             next_id = max(next_id, next_id_a) + 1
#         end
#         # matching_map returns empty dict if the inputs are empty dicts
#         rmap = matching_map!(a₊, a₋, matcher; next_id, ds, p, pprev)
#         push!(rmaps, rmap)
#     end
#     return rmaps
# end

# # Another concrete implementation of `match_sequentially!`:
# function _rematch_with_past!(
#         attractors_cont, matcher;
#         ds = nothing, pcurve = eachindex(attractors_cont),
#     )
#     # this dictionary stores all instances of previous attractors and is updated
#     # at every step. It is then given to the matching function as if it was
#     # the current attractors
#     latest_ghosts = deepcopy(attractors_cont[1])
#     rmaps = Dict{keytype(attractors_cont[1]), keytype(attractors_cont[1])}[]
#     for i in 1:(length(attractors_cont) - 1)
#         a₊, a₋ = attractors_cont[i + 1], attractors_cont[i]
#         # first update ghosts
#         for (k, A) in a₋
#             latest_ghosts[k] = A
#         end
#         # and then match
#         p, pprev = pcurve[i + 1], pcurve[i]
#         rmap = matching_map!(a₊, latest_ghosts, matcher; pprev, p, ds)
#         push!(rmaps, rmap)
#     end
#     return rmaps
# end

_use_vanished(matcher) = false
_retract_keys(matcher) = true

include("basin_overlap.jl")
include("sssdistance.jl")
include("basin_enclosure.jl")
