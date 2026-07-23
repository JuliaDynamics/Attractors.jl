export InitialConditionSampler, generate_ics, update_sampler!
export PerParameterInitialConditions

"""
    InitialConditionSampler

Data structure deciding how to sample initial conditions during
[`global_continuation`](@ref). It defines an extendable interface
based on the internal functions [`generate_ics`](@ref) and [`inform_sampler!`](@ref).

Concerete subtypes are:

- [`RandomICSampler`]
- [`PrescribedICs`]
- [`PerParameterICs`](@ref)
- [`BayesianUpdateSampler`](@ref)


"""
abstract type InitialConditionSampler end

# what should be inputs here? parameters for sure, anything else?
function generate_ics end

function update_sampler!(sampler, args...)
    return nothing
end

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
struct RandomICSampler::F <: InitialConditionSampler
    f::F
    N::Int
end

RandomICSampler(N::Int, args...; kw...) = RandomICSampler(statespace_sampler(args...; kw...)[1], N)

"""
    PrescribedICs(u0s::AbstractVector) <: InitialConditionSampler

Wrapper around a container of initial conditions that simply provides
`u0s` as the sampled initial conditions.
"""
struct PrescribedICs{V<:AbstractVector} <: InitialConditionSampler
    ics::V
end

"""
    PerParameterICs(f, N::Int) <: InitialConditionSampler

Wrapper around a function `f`, to be called as
`f(parameters) -> u`.
It inputs the current parameter(s) of a [`global_continuation`](@ref)
(elements of `pcurve` which are always a dictionary),
and ouputs an initial condition. It generates `N` initial conditions
at a given parameter.
"""
struct PerParameterInitialConditions{F}
    generator::F
    N::Int
end

"""
    Alex's paper.
"""
struct BayesianUpdateSampler::F <: InitialConditionSampler
    f::F
    Ndense::Int
    Nsparse::Int
    γ::Float64
end
