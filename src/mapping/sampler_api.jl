export InitialConditionsSampler
export RandomICsSampler, PrescribedICs, PerParameterICs

"""
    InitialConditionsSampler

Data structure deciding how to sample initial conditions during
[`global_continuation`](@ref).
Concerete subtypes are:

- [`RandomICsSampler`](@ref)
- [`PrescribedICs`](@ref)
- [`PerParameterICs`](@ref)

`InitialConditionsSampler` defines a currently experimental extendable interface
based on the internal functions `generate_ics, update_sampler!, resampling_required`.

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
# function RandomICsSampler(N::Int, args...; kw...)
#     # we have to acknowledge that f() operates in-place here!
#     f, = statespace_sampler(args...; kw...)[1]
#     g = () -> copy(f())
#     return RandomICsSampler(g, N)
# end



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

"""
    BayesianUpdateSampler(dense_sampler, sparse_sampler, γ = 0.5) <: InitialConditionsSampler

From Alex's paper.
"""
struct BayesianUpdateSampler{D, S, R <: HRectangle} <: InitialConditionsSampler
    dense_n::Int
    sparse_n::Int
    γ::Float64
    β::Float64
    λ::Float64
    boxes::Vector{R}
    boxes_flags::Vector{Bool} # same size as `boxes`, true if resampling is required for this box
    resampling_necessary::Bool
    more_fields
end

# Sketch for boxes
# each box has either `dense_n` or `sparse_n`

# if you have a total of `N` i.c.
# and you know the `boxes_flags`, then you know exactly which ones of the total `N`
# are in each box.

# This way, you can generate a _single_ vector of initial conditions
# (although its size would change at each step of the continuation depending
# on which boxes get sparse or dense `n`)

function generate_ics(sampler::BayesianUpdateSampler)
    # This function utilizes `resampling_necessary`.
    all_ics = []
    # loop through boxes
    for (boxi, box) in enumerate(sampler.boxes)
        # for each box, compute η
        η = somehow
        if η < 0
            ics = random_ics_from_box(box, dense_n)
        else
            ics = random_ics_from_box(box, sparse_n)
        end
        append!(all_ics, ics)
    end
    return all_ics
end

function update_sampler!(sampler::BayesianUpdateSampler, labels, args...)
    # first, analyze labels per-box
    counter = 1
    for box_index in eachindex(sampler.boxes)
        n_box = boxes_flags[box_index] ? sampler.dense_n : sampler.sparse_n
        labels_for_box = view(labels, counter:(counter + n_box))
        # then make a decision on the η and box flag and whatever
        # update bayesian
        η = update_somehow(labels_for_box)
        # then update the box flag
        boxes_flags[box_index] = η < 0 ? true : false
        counter += n_box
    end
    # Processing stuff
    return nothing
end

resampling_required(sampler::BayesianUpdateSampler) = any(sampler.boxes_flags)