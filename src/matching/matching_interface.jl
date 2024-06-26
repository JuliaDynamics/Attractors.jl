"""
    SSSetMatcher

Supertype of all "matchers" that match IDs between state space sets (typically attractors).
This is for example useful after performing a [`continuation`](@ref).

Matchers implement an extendable interface. A new mapper type only needs to
extend the function [`replacement_map`](@ref).

As a user, you typically want to use the higher
level function [`match_continuation!`](@ref) directly.

Currently available matchers are:

- [`MatchBySSDistance`](@ref)
- [`MatchByBasinEnclosure`](@ref)
"""
abstract type SSSetMatcher end

"""
    replacement_map(a₊, a₋, matcher; pi = nothing) → rmap

Given dictionaries `a₊, a₋` mapping IDs to `StateSpaceSet` instances,
create a _replacement map_: a dictionary mapping the IDs (keys) in dictionary `a₊`
to IDs (keys) in dictionary `a₋`, so that
so that sets in `a₊` that are the "closest" to sets in `a₋` get assigned the
same key as in `a₋`.
Use [`swap_dict_keys`](@ref) to apply `rmap` to `a₊`.

How matching happens, i.e., how "closeness" is defined, depends on the algorithm `matcher`.

Typically the +,- mean after and before some change of parameter of a system.
For some matchers, such as [`MatchByBasinEnclosure`](@ref), the value of the parameter
is important. `pi` in this case can be given as the index of the parameter
corresponding to the result `a₋`, assuming `a₊, a₋` are two subsequent results
of a [`continuation`](@ref) output.

Return the replacement map, a dictionary mapping old keys of `a₊` to
the new ones that they were mapped to. You can obtain this map, without modifying
the dictionaries, by calling the [`replacement_map`](@ref) function directly.
"""
function replacement_map(a₊::AbstractDict, a₋, matcher::SSSetMatcher; kw...)
    throw(ArgumentError("Not implemented for $(typeof(matcher))"))
end

"""
    match_statespacesets!(a₊, a₋, matcher; pi = nothing) → rmap

Convenience function that estimates the [`replacement_map`](@ref)
and then applies it to `a₊` using the [`swap_dict_keys!`](@ref) function.
"""
function match_statespacesets!(a₊::AbstractDict, a₋, matcher::SSSetMatcher; kw...)
    rmap = replacement_map(a₊, a₋, matcher; kw...)
    swap_dict_keys!(a₊, rmap)
    return rmap
end

# TODO: is it possible to make the "ghost" argument generic that applies
# to all matchers? if so it can be a keyword of match continuation.
# I need to see what it does exactly.
# The retract keys ca be universal and given exclusively to `match_continuation`.

"""
    match_continuation!(set_dicts::Vector{Dict}, matcher::SSSetMatcher) → rmaps

Match the `set_dicts`, a vector of dictionaries mapping IDs to
state space sets (typically the output of [`continuation`](@ref)), according
to the given `matcher` which is any subtype of [`SSSetMatcher`](@ref).
`attractors` is a vector of dictionaries. Each dictionary maps attractor ids
to state space sets representing the attractors.

Return `rmaps`, which is a vector of dictionaries.
`rmaps[i]` contains the [`replacement_map`](@ref) for `attractors[i+1]`,
i.e., the pairs of `old => new` IDs.
"""
function match_continuation!(attractors::AbstractVector{<:Dict}, matcher::SSSetMatcher)
    # this generic implementation works for any matcher which only depends on
    # current and previous attractor state
    rmaps = Dict{Int,Int}[]
    for i in 2:length(attractors)
        a₊, a₋ = attractors[i], attractors[i - 1]
        rmap = match_sssets!(a₊, a₋, matcher)
        push!(rmaps, rmap)
    end
    return rmaps
end

"""
    match_continuation!(continuation_quantity::Vector{Dict}, rmaps::Vector{Dict})

Do the same as in `match_continuation!` above, now given the vector of replacement maps,
and for any arbitrary quantity that has been tracked in the continuation.
`continuation_quantity` is most typically `fractions_curves` from [`continuation`](@ref).
"""
function match_continuation!(continuation_quantity::AbstractVector{<:Dict}, rmaps::Vector{Dict{Int, Int}})
    for (quantity, rmap) in zip(continuation_quantity, rmaps)
        swap_dict_keys!(quantity, rmap)
    end
    return continuation_quantity
end
