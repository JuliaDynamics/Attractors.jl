export global_continuation, GlobalContinuationAlgorithm

"""
    GlobalContinuationAlgorithm

Supertype of all algorithms used in [`global_continuation`](@ref).
Each algorithm typically references an [`AttractorMapper`](@ref),
as well as contains more information for how to continue/track/match attractors
across a parameter range.

See [`global_continuation`](@ref) for more.
"""
abstract type GlobalContinuationAlgorithm end

"""
    global_continuation(gca::GlobalContinuationAlgorithm, prange, pidx, ics; kwargs...)

Find and continue attractors (or representations of attractors)
and the fractions of their basins of attraction across a parameter range.
`global_continuation` is the central function of the framework for global stability analysis
illustrated in [Datseris2023](@cite).

The global_continuation type `gca` is a subtype of `GlobalContinuationAlgorithm`
and references an [`AttractorMapper`](@ref). The mapper contains information
on how to find the attractors and basins of a dynamical system. Additional
arguments that control how to continue/track/match attractors
are given when creating `gca`.

The basin fractions and the attractors (or some representation of them) are continued
across the parameter range `prange`, for the parameter of the system with index `pidx`
(any index valid in `DynamicalSystems.set_parameter!` can be used).

`ics` are the initial conditions to use when globally sampling the state space.
Like in [`basins_fractions`](@ref) it can be either a set vector of initial conditions,
or a 0-argument function that generates random initial conditions.

Possible subtypes of `GlobalContinuationAlgorithm` are:

- [`RecurrencesFindAndMatch`](@ref)
- [`FeaturizeGroupAcrossParameter`](@ref)

## Return

1. `fractions_cont::Vector{Dict{Int, Float64}}`. The fractions of basins of attraction.
   `fractions_cont[i]` is a dictionary mapping attractor IDs to their basin fraction
   at the `i`-th parameter.
2. `attractors_cont::Vector{Dict{Int, <:Any}}`. Information about the attractors.
   `attractors_cont[i]` is a dictionary mapping attractor ID to information about the
   attractor at the `i`-th parameter.
   The type of information stored depends on the chosen global_continuation type.

## Keyword arguments

- `show_progress = true`: display a progress bar of the computation.
- `samples_per_parameter = 100`: amount of initial conditions sampled at each parameter
  from `ics` if `ics` is a function instead of set initial conditions.
"""
function global_continuation end

include("continuation_afam.jl")
include("continuation_recurrences.jl")
include("continuation_grouping.jl")
include("aggregate_attractor_fractions.jl")