export InitialConditionsSampler
export RandomICsSampler, PrescribedICs, PerParameterICs

using Random: Xoshiro

"""
    InitialConditionsSampler

Data structure deciding how to sample initial conditions during
[`global_continuation`](@ref).
Concerete subtypes are:

- [`RandomICsSampler`](@ref)
- [`PrescribedICs`](@ref)
- [`PerParameterICs`](@ref)
- [`BayesianUpdateSampler`](@ref)

`InitialConditionsSampler` defines a currently experimental extendable interface
based on the internal functions
`generate_ics, update_sampler!, resampling_required, weighted_fractions`.

`length(sampler)` returns the number of initial conditions that will be generated
by the sampler.
"""
abstract type InitialConditionsSampler end
Base.length(s::InitialConditionsSampler) = s.N

"""
    generate_ics(sampler::InitialConditionsSampler, params)

Generate initial condititions from the given sampler, optionally utilizing
a container of parameters of a dynamical system.
"""
function generate_ics end

"""
    update_sampler!(sampler::InitialConditionsSampler, args...)

Todo, decide `args`.
"""
update_sampler!(sampler, args...) = nothing
resampling_required(sampler) = false

"""
    weighted_fractions(sampler::InitialConditionsSampler, counts) → Dict{Int, Float64}

Since the sampling is not uniform across the domain the fractions have to be computed
taking into account the weight of each box.

The default implementation is just `counts ./ sum(counts)`, which is right
whenever the sampler covers the region uniformly and in one round.
"""
function weighted_fractions(sampler, counts)
    n = sum(values(counts); init = 0)
    return Dict{Int, Float64}(k => c / n for (k, c) in counts)
end

"""
    RandomICsSampler(f::Function, N::Int) <: InitialConditionsSampler

Wrapper around a function `f`, to be called as
`f() -> u`. When called, it generates a random initial condition.
The sampler generates overall `N` initial conditions.

The following convenience signature is also provided

    RandomICsSampler(N::Int, args...; kw...)

which propagates `args, kw` to [`statespace_sampler`](@ref) and uses the generated
sampler as the function `f`.
"""
struct RandomICsSampler{F} <: InitialConditionsSampler
    f::F
    N::Int
end
RandomICsSampler(N::Int, args...; kw...) = RandomICsSampler(statespace_sampler(args...; kw...)[1], N)
generate_ics(p::RandomICsSampler, args...) = (p.f() for _ in 1:p.N)

"""
    PrescribedICs(u0s::AbstractVector) <: InitialConditionsSampler

Wrapper around a container of initial conditions that simply provides
`u0s` as the sampled initial conditions.
"""
struct PrescribedICs{V<:AbstractVector} <: InitialConditionsSampler
    ics::V
end
generate_ics(p::PrescribedICs, args...) = p.ics
Base.length(s::PrescribedICs) = length(s.ics)

"""
    PerParameterICs(f, N::Int) <: InitialConditionsSampler

Wrapper around a function `f`, to be called as
`f(parameters, N)`. When used in [`basins_fractions`](@ref),
it inputs the `current_parameters` of the dynamical system.
When used in [`global_continuation`](@ref) it inputs the current
element of `pcurve` (which is expected to be a dictionary).
The sampler generates overall `N` initial conditions.
"""
struct PerParameterICs{F} <: InitialConditionsSampler
    f::F
    N::Int
end
generate_ics(p::PerParameterICs, params, args...) = p.f(params, p.N)

include("bayesian.jl")