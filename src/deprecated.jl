# This function is deprecated, however it will NEVER be removed.
# In the future, this call signature will remain here but its source code
# will simply be replaced by an error message.
function basins_of_attraction(grid::Tuple, ds::DynamicalSystem; kwargs...)
    @warn("""
    The function `basins_of_attraction(grid::Tuple, ds::DynamicalSystem; ...)` has
    been replaced by the more generic
    `basins_of_attraction(mapper::AttractorMapper, grid::Tuple`) which works for
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

@deprecate RecurrencesSeededContinuation RecurrencesFindAndMatch
@deprecate GroupAcrossParameterContinuation FeaturizeGroupAcrossParameter
@deprecate match_attractor_ids! match_statespacesets!
@deprecate GroupAcrossParameter FeaturizeGroupAcrossParameter
@deprecate rematch! match_continuation!

function match_statespacesets!(as::Vector{<:Dict}; kwargs...)
    error("This function was incorrect. Use `match_continuation!` instead.")
end