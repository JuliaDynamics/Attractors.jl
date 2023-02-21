export continuation

"""
   AttractorsBasinsContinuation

Supertype of continuation methods given to [`continuation`](@ref).
Possible subtypes are:

- [`RecurrencesContinuation`](@ref)
- [`FeaturizingContinuation`](@ref)

These are given to the main [`continuation`](@ref) function.
"""
abstract type AttractorsBasinsContinuation end

"""
   MatchingMethod

Supertype of types specifying how to match attractors/features in [`continuation`](@ref).
Concrete subtypes are given as options when instantiating an
[`AttractorsBasinsContinuation`](@ref) subtype.

Possible subtypes are:

- [`ParameterSliceCrossDistance`](@ref)
- [`ClusterOverAllParameters`](@ref)
"""
abstract type MatchingMethod end

"""
   continuation(abc::AttractorsBasinsContinuation, prange, pidx, ics; kwargs...)

Find and continue attractors and the fractions of their basins of attraction
across a parameter range.

The continuation type `abc` is a subtype of [`AttractorsBasinsContinuation`](@ref)
and contains an [`AttractorMapper`](@ref). The mapper contains information
on how to find the attractors and basins of a dynamical system. Additional
arguments and keyword arguments given when creating `abc` further tune the continuation.

In the `continuation` function, the basin fractions and the attractors (or some representation
of them in the case of featurizing) are continued across the parameter range `prange`,
for the parameter of the system with index `pidx`.
`ics` is as in [`basins_fractions`](@ref), i.e., it is either a function generating
initial conditions or a set containing them.

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
* `samples_per_parameter = 100`: amount of initial conditions sampled at each parameter.
- `cont_method = :grouping`: selects the method to perform the continuation. `:grouping`
  is meant to group the features accross the parameter range while `:matching` will match
  the clusters of attractors from one parameter slice to the next.
"""
function continuation end

include("match_attractor_ids.jl")
include("continuation_recurrences.jl")
include("continuation_featurizing.jl")
include("aggregate_attractor_fractions.jl")
