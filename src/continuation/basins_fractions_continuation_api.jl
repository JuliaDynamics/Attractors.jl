export basins_fractions_continuation

# In the end, it is better to have a continuation type that contains
# how to match, because there are other keywords that always go into the
# continuation... Like in the recurrences the keyword of seeds per attractor,
# or in the clustering some other stuff like the parameter scaling
abstract type AttractorsBasinsContinuation end

"""
    basins_fractions_continuation(continuation, prange, pidx, ics; kwargs...)

Find and continue attractors and the fractions of their basins of attraction
across a parameter range.
The given `continuation` contains a reference to a dynamical system,
as well as how to find its attractors. I.e., it contains an [`AttractorMapper`](@ref).
Given this `continuation`, the basin fractions (and the attractors for the
`RecurrencesContinuation` method) are continued across the parameter range `prange`,
for the parameter of the system with index `pidx`.
`ics` is as in [`basins_fractions`](@ref), i.e., it is either a function generating
initial conditions or a dataset containing them.

## Return

1. `fractions_curves <: Vector{Dict{Int, Float64}}`. The fractions of basins of attraction.
   `fractions_curves[i]` is a dictionary mapping attractor IDs to their basin fraction
   at the `i`-th parameter.
2. `attractors_info <: Vector{Dict{Int, <:Any}}`. Information about the attractors.
   `attractors_info[i]` is a dictionary mapping attractor ID to information about the
   attractor at the `i`-th parameter.
   The type of information stored depends on the chosen continuation method.

## Keyword arguments

* `samples_per_parameter = 100`: Amount of initial conditions sampled at each parameter.
* `show_progress = true`: Whether to show a progress bar for the process.

## Continuation methods

- [`RecurrencesContinuation`](@ref). For this sampler, `ics` is optional.
  If not given, one is created using the `grid` of [`AttractorsViaRecurrences`](@ref):
  ```
  sampler, = statespace_sampler(min_bounds = minimum.(grid), max_bounds = maximum.(grid))
  ```
- [`FeaturizingContinuation`](@ref).

"""
function basins_fractions_continuation end

include("match_attractor_ids.jl")
include("continuation_recurrences.jl")
include("continuation_featurizing.jl")
include("aggregate_attractor_fractions.jl")
