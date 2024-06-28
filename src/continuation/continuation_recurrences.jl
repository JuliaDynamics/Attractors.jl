export RecurrencesFindAndMatch, RAFM

"""
    RecurrencesFindAndMatch <: AttractorsBasinsContinuation
    RecurrencesFindAndMatch(mapper::AttractorsViaRecurrences; kwargs...)

A method for [`continuation`](@ref) as in [Datseris2023](@cite) that is based on the
recurrences algorithm for finding attractors ([`AttractorsViaRecurrences`](@ref))
and then matching them according to their state space distance.

## Keyword arguments

- `distance = Centroid(), threshold = Inf`: passed to [`MatchBySSSetDistance`](@ref).
- `seeds_from_attractor`: A function that takes as an input an attractor and returns
  an iterator of initial conditions to be seeded from the attractor for the next
  parameter slice. By default, we sample only the first stored point on the attractor.

## Description

`RecurrencesFindAndMatch` is a wrapper type. It is has been generalized by
[`AttractorsContinueAndMatch`](@ref). It is still exported for backwards compatibility
and to have a clear reference to the original algorithm developed in [Datseris2023](@cite).

The source code of `RecurrencesFindAndMatch` is trival:
it takes the given mapper, it initializes a [`MatchBySSSetDistance`](@ref),
and along with `seeds_from_attractor` it makes the [`AttractorsContinueAndMatch`](@ref)
instance. This is the process described in [Datseris2023](@cite),
whereby attractors are found using the recurrences algorithm [`AttractorsViaRecurrences`](@ref)
and they are then matched by their distance in state space [`MatchBySSSetDistance`](@ref).
"""
function RecurrencesFindAndMatch(
        mapper::AttractorsViaRecurrences; distance = Centroid(), threshold = Inf,
        info_extraction = nothing, seeds_from_attractor = _default_seeding_process,
    )
    if info_extraction !== nothing
        @warn "`info_extraction` is ignored in `RecurrencesFindAndMatch`.
        You can extract info after the attractors have been found."
    end
    matcher = MatchBySSSetDistance(; distance, threshold)
    return AttractorsContinueAndMatch(mapper, matcher, seeds_from_attractor)
end

"Alias for [`RecurrencesFindAndMatch`](@ref)."
const RAFM = RecurrencesFindAndMatch
