export InitialConditionSampler, generate_ics, update_sampler!
export RandomICSampler, PerParameterICs, PerParameterInitialConditions

"""
    InitialConditionSampler

Data structure deciding how to sample initial conditions during
[`global_continuation`](@ref). It defines an extendable interface
based on the internal functions [`generate_ics`](@ref) and [`inform_sampler!`](@ref).

Concerete subtypes are:

- [`RandomICSampler`]
- [`PrescribedICs`]
- [`PerParameterICs`](@ref)
"""
abstract type InitialConditionSampler end
Base.length(s::InitialConditionSampler) = s.N

"""
    generate_ics(sampler::InitialConditionsSampler, p::Dict)

Generate initial condititions from the given sampler, optionally utilizing
the parameters of the current continuation step.
"""
function generate_ics end

"""
    update_sampler!(sampler::InitialConditionsSampler, args...)

Todo, decide `args`.
"""
update_sampler!(sampler, args...) = nothing

"""
    RandomICSampler(f::Function, N::Int) <: InitialConditionSampler

Wrapper around a function `f`, to be called as
`f() -> u`. When called, it generates a random initial condition.
The sampler generates overall `N` initial conditions.

The following convenience signature is also provided

    RandomICSampler(N::Int, args...; kw...)

which propagates `args, kw` to [`statespace_sampler`](@ref) and uses the generated
sampler as the function `f`.
"""
struct RandomICSampler{F} <: InitialConditionSampler
    f::F
    N::Int
end
RandomICSampler(N::Int, args...; kw...) = RandomICSampler(statespace_sampler(args...; kw...)[1], N)
generate_ics(p::RandomICSampler, args...) = (p.f() for _ in 1:p.N)

"""
    PrescribedICs(u0s::AbstractVector) <: InitialConditionSampler

Wrapper around a container of initial conditions that simply provides
`u0s` as the sampled initial conditions.
"""
struct PrescribedICs{V<:AbstractVector} <: InitialConditionSampler
    ics::V
end
generate_ics(p::PrescribedICs, args...) = p.ics
Base.length(s::PrescribedICs) = length(s.ics)

"""
    PerParameterICs(f, N::Int) <: InitialConditionSampler

Wrapper around a function `f`, to be called as
`f(parameters, N)`.
It inputs the current parameter(s) of a [`global_continuation`](@ref)
(elements of `pcurve` which are always a dictionary),
and ouputs an iterable of `N` initial conditions.
The sampler generates overall `N` initial conditions.
"""
struct PerParameterInitialConditions{F} <: InitialConditionSampler
    f::F
    N::Int
end
generate_ics(p::PerParameterInitialConditions, params, args...) = p.f(params, p.N)

"""
Alex's paper.
"""
struct BayesianUpdateSampler{F} <: InitialConditionSampler
    f::F
    Ndense::Int
    Nsparse::Int
    γ::Float64
end
