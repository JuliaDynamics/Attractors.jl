export global_continuation, GlobalContinuationAlgorithm, continuation_series, stability_quantifiers_along_continuation
include("sampler_api.jl")

"""
    GlobalContinuationAlgorithm

Supertype of all algorithms used in [`global_continuation`](@ref).
Each algorithm typically references a [`BasinMap`](@ref),
as well as contains more information for how to continue/track/match attractors
across a parameter range.

See [`global_continuation`](@ref) for more.
"""
abstract type GlobalContinuationAlgorithm end

"""
    global_continuation(gca::GlobalContinuationAlgorithm, pcurve, icsampler; kwargs...)

Find and continue attractors (or representations of attractors) and properties of their
basins across a parameter curve `pcurve` according to algorithm `gca` and
by sampling initial conditions using `icsampler`.

Possible subtypes of a `GlobalContinuationAlgorithm` are:

- [`AttractorSeedContinueMatch`](@ref)
- [`FeaturizeGroupAcrossParameter`](@ref)

`pcurve` is a `Vector` of dictionaries, each dictionary mapping parameter indices to values.
This defines an arbitrary curve in the parameter space of the dynamical system.

`icsampler` is a subtype of [`InitialConditionSampler`](@ref) and provides instructions
for how to sample initial conditions to explore the state space during the continuation.

Return:

1. `fractions_cont::Vector{Dict{Int, Float64}}`. The fractions of basins of attraction.
   `fractions_cont[i]` is a dictionary mapping basin IDs to their basin fraction
   at the `i`-th parameter combination.
   -  This output is different if you are using [`StabilityQuantifiersAccumulator`](@ref)
      in combination with [`AttractorSeedContinueMatch`](@ref). See the docstring
      of [`StabilityQuantifiersAccumulator`](@ref) for more details.
2. `attractors_cont::Vector{Dict{Int, SSSet}}`. The continued attractors.
   `attractors_cont[i]` is a dictionary mapping basin ID to its
   attractor set at the `i`-th parameter combination.

See the function [`continuation_series`](@ref) if you wish to transform the output(s)
to an alternative format. There is no difference between single or multi parameter
global continuation. Use [`hilbert_pcurve`](@ref) to cover multiparameter spaces.

## Keyword arguments

- `show_progress = true`: display a progress bar of the computation.

## Description

`global_continuation` is the central function for performing global continuation
as outlined in our article [Datseris2026](@cite).

The global continuation algorithm typically references a [`BasinMap`](@ref)
which is used to find basins and corresponding attractors of a dynamical system. Additional
arguments that control how to continue/track/match attractors across a parameter range
are given when creating `gca`.

The basin properties and the attractors (or some representation of them) are continued
across the parameter curve, whose elements are simply given to `DynamicalSystems.set_parameters!`
to update system parameters. This is fundamentally different to local (traditional) continuation,
see the online documentation or our article for details.
"""
function global_continuation end

"""
    continuation_series(continuation_info, fillval = NaN)

Transform a continuation quantity (a vector of dictionaries, each dictionary
mapping basin IDs to values of same type as `fillval`) to a dictionary of vectors
where the `k` dictionary entry is the series of the continuation quantity corresponding to
basin with ID `k`. `fillval` denotes the value to assign in the series
if the basin with ID `k` does not exist at this particular series index.
If the `continuation_info` is the attractors themselves, you likely want to
use as `fillval` some empty state space set such as `StateSpaceSet{D}()`,
or just the `missing` value (in which case the output will be of `Union` type).

This function is typically used on the output of [`global_continuation`](@ref).
"""
function continuation_series(continuation_info::AbstractVector{<:AbstractDict}, fillval = NaN, ukeys = unique_keys(continuation_info))
    V = valtype(eltype(continuation_info))
    T = promote_type(typeof(fillval), V)
    series = Dict(k => fill!(Vector{T}(undef, length(continuation_info)), fillval) for k in ukeys)
    for i in eachindex(continuation_info)
        for k in ukeys
            series[k][i] = get(continuation_info[i], k, fillval)
        end
    end
    return series
end

include("hilbert_pcurve.jl")
include("continuation_ascm_generic.jl")
include("continuation_recurrences.jl")
include("continuation_grouping.jl")
include("aggregate_continuation.jl")
include("continuation_stability_quantifiers.jl")

