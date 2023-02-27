export continuation, AttractorsBasinsContinuation

# In the end, it is better to have a continuation type that contains
# how to match, because there are other keywords that always go into the
# continuation... Like in the recurrences the keyword of seeds per attractor,
# or in the clustering some other stuff like the parameter scaling
abstract type AttractorsBasinsContinuation end

"""
    continuation(abc::AttractorsBasinsContinuation, prange, pidx, ics; kwargs...)

Find and continue attractors (or feature-based representations of attractors)
and the fractions of their basins of attraction across a parameter range.

The continuation type `abc` is a subtype of `AttractorsBasinsContinuation`
and contains an [`AttractorMapper`](@ref). The mapper contains information
on how to find the attractors and basins of a dynamical system. Additional
arguments and keyword arguments given when creating `abc` further tune the continuation
and how attractors are matched across different parameter values.

The basin fractions and the attractors (or some representation of them) are continued
across the parameter range `prange`, for the parameter of the system with index `pidx`.
`ics` is as in [`basins_fractions`](@ref), i.e., it is either a function generating
initial conditions or a set containing them.

`ics` is the initial conditions as in [`basins_fractions`](@ref).

Possible subtypes of `AttractorsBasinsContinuation` are:

- [`RecurrencesSeededContinuation`](@ref)
- [`GroupAcrossParameterContinuation`](@ref)

## Return

1. `fractions_curves :: Vector{Dict{Int, Float64}}`. The fractions of basins of attraction.
   `fractions_curves[i]` is a dictionary mapping attractor IDs to their basin fraction
   at the `i`-th parameter.
2. `attractors_info <: Vector{Dict{Int, <:Any}}`. Information about the attractors.
   `attractors_info[i]` is a dictionary mapping attractor ID to information about the
   attractor at the `i`-th parameter.
   The type of information stored depends on the chosen continuation type.

## Keyword arguments

- `show_progress = true`: display a progress bar of the computation.
* `samples_per_parameter = 100`: amount of initial conditions sampled at each parameter
  if `ics` is a function.
"""
function continuation end

include("match_attractor_ids.jl")
include("continuation_recurrences.jl")
include("continuation_grouping.jl")
include("aggregate_attractor_fractions.jl")