export global_continuation, GlobalContinuationAlgorithm, continuation_series
export GlobalContinuationOutput

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
    GlobalContinuationOutput

The type of output of [`global_continuation`](@ref). It contains the fields:

- `attractors::Vector{Dict{Int, SSSet}}`. The continued and matched attractors.
  `attractors[i]` is a dictionary mapping basin ID to its
  attractor set at the `i`-th parameter combination.
- `fractions::Vector{Dict{Int, Float64}}`. The fractions of basins of attraction.
  `fractions[i]` is a dictionary mapping basin IDs to their basin fraction
  at the `i`-th parameter combination.
- `quantifiers::Dict{String, Vector{Dict}}`. Quantifiers of the attractors
  or their basins that have been (potentially) continued. This entry is empty,
  unless continuation is done with [`StabilityQuantifiersAccumulator`](@ref).
  Then, it is a dictionary mapping strings (names of quantifiers) to vectors of dictionaries.
- `other::Dict{String, Vector}`. Any other information was generated during the continuation,
  that may relate to the sampling strategy, featurizing, or other.
  Each quantity's name (string) is mapped to a vector of the same length as `pcurve`.
- `pcurve::Vector{<:Dict}`. The parameter curve the continuation continued over.

The various algorithms used in global continuation inform in their documentation strings
about additional information added to this output.

See the function [`continuation_series`](@ref) if you wish to transform the output(s)
to an alternative format.
"""
struct GlobalContinuationOutput{SSS, F<:AbstractFloat, V<:Vector, A, P}
    attractors::Vector{Dict{Int, SSS}}
    fractions::Vector{Dict{Int, F}}
    quantifiers::Dict{String, V}
    other::Dict{String, A}
    pcurve::Vector{P}
end
# this type extents `iterate` in deprecated.jl file.

function Base.show(io::IO, gco::GlobalContinuationOutput)
    println(io, "GlobalContinuationOutput with fields:")
    println(io, " attractors")
    println(io, " fractions")
    println(io, " quantifiers")
    println(io, " other")
    println(io, " pcurve")
end


# Internal function for adding info; dispatches on `something`
function add_extra_continuation_info!(extras::Dict{String,Any}, something)
    return nothing
end

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

`icsampler` is a subtype of [`InitialConditionsSampler`](@ref) and provides instructions
for how to sample initial conditions to explore the state space during the continuation.

Return an instance of [`GlobalContinuationOutput`](@ref) that contains the continued
attractors, their basins fractions, and any other additional information.

There is no difference between single or multi parameter
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
across the parameter curve, whose elements are given to `DynamicalSystems.set_parameters!`
to update system parameters. Thus, global continuation operates on a prescribed parameter
curve.
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

