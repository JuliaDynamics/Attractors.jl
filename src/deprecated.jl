# This function is deprecated, however it will NEVER be removed.
# In the future, this call signature will remain here but its source code
# will simply be replaced by an error message.
function basins_of_attraction(grid::Tuple, ds::DynamicalSystem; kwargs...)
    error("""
    The function `basins_of_attraction(grid::Tuple, ds::DynamicalSystem; ...)` has
    been replaced by the more generic
    `basins_of_attraction(mapper::AttractorMapper, grid::Tuple`) which works for
    any instance of `AttractorMapper`. The `AttractorMapper` itself requires as
    input an kind of dynamical system the user wants, like a `StroboscopicMap` or
    `CoupledODEs` or `DeterministicIteratedMap` etc.

    Please use that method instead.
    """)
end


@deprecate RecurrencesSeededContinuation RecurrencesFindAndMatch
@deprecate GroupAcrossParameterContinuation GroupAcrossParameter