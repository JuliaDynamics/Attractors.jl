# This function is will NEVER be removed. It refers to original article code.
function basins_of_attraction(grid::Tuple, ds::DynamicalSystem; kwargs...)
    error(
        """
        The function `basins_of_attraction(grid::Tuple, ds::DynamicalSystem; ...)` has
        been replaced by the more generic
        `basins_of_attraction(mapper::BasinMap, grid::Tuple)` which works for
        any instance of `BasinMap`. The `BasinMap` itself requires as
        input an kind of dynamical system the user wants, like a `StroboscopicMap` or
        `CoupledODEs` or `DeterministicIteratedMap` etc.

        For now, we do the following for you:
        ```
        mapper = BasinMapRecurrences(ds, grid; sparse = false)
        basins_of_attraction(mapper)
        ```
        and we are completely ignoring any keywords you provided (which could be about the
        differential equation solve, or the metaparameters of the recurrences algorithm).

        We strongly recommend that you study the documentation of Attractors.jl
        and update your code. The only reason we provide this backwards compatibility
        is because our first paper "Effortless estimation of basins of attraction"
        uses this function signature in the script in the paper (which we can't change anymore).
        """
    )
    mapper = BasinMapRecurrences(ds, grid; sparse = false)
    return basins_of_attraction(mapper)
end

@deprecate AttractorsViaRecurrences BasinMapRecurrences
@deprecate AttractorsViaFeaturizing BasinMapFeaturizeGroup
@deprecate AttractorsViaProximity BasinMapProximity
@deprecate stability_measures_along_continuation stability_quantifiers_along_continuation
@deprecate StabilityMeasuresAccumulator StabilityQuantifiersAccumulator
