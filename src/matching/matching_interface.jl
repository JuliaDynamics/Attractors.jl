"""
    SSSetMatcher

Supertype of all "matchers" that match IDs between state space sets (typically attractors).
This is for example useful after performing a [`continuation`](@ref).

Matchers implement the core [`match_statespacesets!`](@ref) function,
however it is typically the case that you want to use the [`match_continuation!`](@ref)
function instead.

Currently available matchers are:

- [`MatchBySSDistance`](@ref)
- [`MatchByBasinEnclosure`](@ref)
"""
abstract type SSSetMatcher end

# TODO: Give [`match_statespacesets!`] the index `i`  that by default is
# `nothing` and some matchers use it.
# Some matchers may reference the parameter vector and dynamical system if
# they need to evolve it.`
# This means that the low level function to actually extend is `match_statespacesets!`
# and not the continuation one!

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
