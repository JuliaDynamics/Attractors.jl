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
by sampling initial conditions using `icsampler`

Possible subtypes of a `GlobalContinuationAlgorithm` are:

- [`AttractorSeedContinueMatch`](@ref)
- [`FeaturizeGroupAcrossParameter`](@ref)

`pcurve` is a vector of dictionaries, each dictionary mapping parameter indices to values.
This defines an arbitrary curve in the whole parameter space of the dynamical system.

`icsampler` is a subtype of [`InitialConditionSampler`](@ref) and provides instructions
for how to sample initial conditions to explore the state space during the continuation.

Return:

1. `fractions_cont::Vector{Dict{Int, Float64}}`. The fractions of basins of attraction.
   `fractions_cont[i]` is a dictionary mapping attractor IDs to their basin fraction
   at the `i`-th parameter combination.
   -  This output is different if you are using [`StabilityQuantifiersAccumulator`](@ref)
      in combination with [`AttractorSeedContinueMatch`](@ref). See the docstring
      of [`StabilityQuantifiersAccumulator`](@ref) for more details.
2. `attractors_cont::Vector{Dict{Int, <:Any}}`. The continued attractors.
   `attractors_cont[i]` is a dictionary mapping attractor ID to the
   attractor set at the `i`-th parameter combination.

See the function [`continuation_series`](@ref) if you wish to transform the output(s)
to an alternative format. There is no difference between single or multi parameter
global continuation. Use [`hilbert_pcurve`](@ref) to cover multiparameter spaces.

## Keyword arguments

- `show_progress = true`: display a progress bar of the computation.

## Description

`global_continuation` is the central function of the framework for global stability analysis
illustrated in [Datseris2023](@cite).

The global continuation algorithm typically references an [`BasinMap`](@ref)
which is used to find the attractors and basins of a dynamical system. Additional
arguments that control how to continue/track/match attractors across a parameter range
are given when creating `gca`.

The basin fractions and the attractors (or some representation of them) are continued
across the parameter range `prange`, for the parameter of the system with index `pidx`
(any index valid in `DynamicalSystems.set_parameter!` can be used). In contrast to
traditional continuation (see online Tutorial for a comparison), global continuation
can be performed over arbitrary user-defined curves in parameter space.
The second call signature with `pcurve` allows for this possibility. In this case
`pcurve` is a vector of iterables, where each iterable maps parameter indices
to parameter values. These iterables can be dictionaries, named tuples, `Vector{Pair}`,
anything that can be given in `set_parameters!`. We recommend to use dictionaries.
Regardless, the sequence of the iterables defines a curve in parameter space.
In fact, the version with `prange, pidx` simply defines
`pcurve = [Dict(pidx => p) for p in prange]` and calls the second method.
"""
function global_continuation(alg::GlobalContinuationAlgorithm, prange::AbstractVector, pidx, sampler; kw...)
    # everything is propagated to the curve setting
    pcurve = [Dict(pidx => p) for p in prange]
    return global_continuation(alg, pcurve, sampler; kw...)
end


"""
    continuation_series(continuation_info, fillval = NaN)

Transform a continuation quantity (a vector of dictionaries, each dictionary
mapping attractor IDs to values of same type as `fillval`) to a dictionary of vectors
where the `k` dictionary entry is the series of the continuation quantity corresponding to
attractor with ID `k`. `fillval` denotes the value to assign in the series
if the attractor with ID `k` does not exist at this particular series index.
If the `continuation_info` is the attractors themselves, you likely want to
use as `fillval` some empty state space set such as `StateSpaceSet{D}()`.
"""
function continuation_series(continuation_info::AbstractVector{<:AbstractDict}, fillval = NaN, ukeys = unique_keys(continuation_info))
    series = Dict(k => fill(fillval, length(continuation_info)) for k in ukeys)
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
