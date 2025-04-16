# This function is deprecated, however it will NEVER be removed.
# In the future, this call signature will remain here but its source code
# will simply be replaced by an error message.
function basins_of_attraction(grid::Tuple, ds::DynamicalSystem; kwargs...)
    @warn("""
    The function `basins_of_attraction(grid::Tuple, ds::DynamicalSystem; ...)` has
    been replaced by the more generic
    `basins_of_attraction(mapper::AttractorMapper, grid::Tuple)` which works for
    any instance of `AttractorMapper`. The `AttractorMapper` itself requires as
    input an kind of dynamical system the user wants, like a `StroboscopicMap` or
    `CoupledODEs` or `DeterministicIteratedMap` etc.

    For now, we do the following for you:
    ```
    mapper = AttractorsViaRecurrences(ds, grid; sparse = false)
    basins_of_attraction(mapper)
    ```
    and we are completely ignoring any keywords you provided (which could be about the
    differential equation solve, or the metaparameters of the recurrences algorithm).

    We strongly recommend that you study the documentation of Attractors.jl
    and update your code. The only reason we provide this backwards compatibility
    is because our first paper "Effortless estimation of basins of attraction"
    uses this function signature in the script in the paper (which we can't change anymore).
    """)
    mapper = AttractorsViaRecurrences(ds, grid; sparse = false)
    return basins_of_attraction(mapper)
end

function continuation(args...; kwargs...)
    @warn("""
    The function `Attractors.continuation` is deprecated in favor of
    `global_continuation`, in preparation for future developments where both
    local/linearized and global continuations will be possible within DynamicalSystems.jl.
    """)
    return global_continuation(args...; kwargs...)
end

export continuation

@deprecate AttractorsBasinsContinuation GlobalContinuationAlgorithm
@deprecate RecurrencesSeededContinuation RecurrencesFindAndMatch
@deprecate GroupAcrossParameterContinuation FeaturizeGroupAcrossParameter
@deprecate match_attractor_ids! match_statespacesets!
@deprecate GroupAcrossParameter FeaturizeGroupAcrossParameter
@deprecate rematch! match_continuation!
@deprecate replacement_map  matching_map

function match_statespacesets!(as::Vector{<:Dict}; kwargs...)
    error("This function was incorrect. Use `match_sequentially!` instead.")
end

export match_continuation!, match_statespacesets!, match_basins_ids!

function match_continuation!(args...; kwargs...)
    @warn "match_continuation! has been deprecated for `match_sequentially!`,
    which has a more generic name that better reflects its capabilities."
    return match_sequentially!(args...; kwargs...)
end

function match_statespacesets!(a_afte, a_befo; kwargs...)
    @warn "match_statespacesets! is deprecated. Use `matching_map!` with `MatchBySSSetDistance`."
    return matching_map!(a_afte, a_befo, MatchBySSSetDistance(kwargs...))
end

function match_basins_ids!(b₊::AbstractArray, b₋; threshold = Inf)
    @warn "`match_basins_ids!` is deprecated, use `matching_map` with `MatchByBasinOverlap`."
    matcher = MatchByBasinOverlap(threshold)
    return matching_map!(b₊, b₋, matcher)
end

@deprecate minimal_fatal_shock minimal_critical_shock
@deprecate MFSBlackBoxOptim MCSBlackBoxOptim
@deprecate MFSBruteForce MCSBruteForce