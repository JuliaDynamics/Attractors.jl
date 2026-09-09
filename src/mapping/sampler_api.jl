export InitialConditionsSampler
export RandomICsSampler, PrescribedICs, PerParameterICs
export BayesianUpdateSampler, sampler_history

using Random: Xoshiro
using SpecialFunctions: loggamma

"""
    InitialConditionsSampler

Data structure deciding how to sample initial conditions during
[`global_continuation`](@ref).
Concerete subtypes are:

- [`RandomICsSampler`](@ref)
- [`PrescribedICs`](@ref)
- [`PerParameterICs`](@ref)

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
    BayesianUpdateSampler(region, n_tiles::Int; sparse_n, kwargs...) <: InitialConditionsSampler

Sampler allocating initial conditions where the basins are actually changing.
`region` is tiled into `n_tiles^D` equally sized boxes, each carrying information
over the attractor labels found inside it. At every parameter of a
[`global_continuation`](@ref) each box is sampled sparsely and the resulting label
counts are tested against its prior with a log Bayes factor `η`. A box with negative
`η` means the data are better explained by no prior at all than by its history.
The basins in the box have changed and the sampler asks for a dense re-sample.

`region` is anything [`statespace_sampler`](@ref) accepts as a region: an `HRectangle`,
or a tuple of ranges/`(min, max)` pairs, one per dimension.

## Keyword arguments

- `sparse_n::Int`: initial conditions drawn per box during routine monitoring.
- `dense_n::Int = sparse_n^2`: initial conditions drawn per box when a box asks for a
  re-sample, and for every box at the first parameter of the continuation.
- `λ::Real = 0.7`: forgetting factor. The prior is decayed as `α ← λα` before each
  sparse update, so that evidence from far-away parameters is progressively discounted.
- `β::Real = 0.5`: Dirichlet base pseudo-count assigned to unseen labels.
- `global_reset::Bool = false`: a heuristic sitting on top of the per-box test.
  If an attractor appear or disapear from one parameter to another we can ask 
  the sampler to flag all boxes for a dense resampling. It will give a better 
  estimate of the basin entropy for example. Basins fraction are more robust and 
  the change will not be so drastic. 
- `seed = abs(rand(Int))`: seed for the per-box point generators.
- `history::Bool = false`: keep a per-parameter record of `alphas` and `etas`, which
  are otherwise overwritten in place. See "History" below.

## Description

At the first parameter every box is sampled with `dense_n` initial conditions and its
prior is initialised as `α_k = c_k + β` from the label counts `c` and atractor k.
At every later parameter:

1. `generate_ics` draws `sparse_n` points per box.
2. `update_sampler!` slices the returned labels per box, decays that box's prior by `λ`,
   and computes `η`. If `η < 0` the box is flagged and its prior left untouched;
   otherwise the posterior `α_k ← λα_k + c_k` is stored.
3. If any box was flagged, `resampling_required` returns `true` and the continuation loop
   calls `generate_ics` again. That round draws `dense_n` points for the flagged boxes
   *only*, and `update_sampler!` rebuilds their priors from scratch and clears the flags.

Because a dense round performs no test, at most one re-sampling round happens per
parameter. Every box starts out flagged, so the dense initialisation of the first
parameter is nothing but step 3 applied to all of them.

## History

`alphas` and `etas` describe the current parameter only: every round overwrites them,
so by the time a continuation returns they say nothing about the sweep that produced
it. With `history = true` the sampler snapshots both at the end of each parameter.

Use [`sampler_history`](@ref)
to read them back.

!!! warning "The history IDs are not those of the continuation output"
    The recorded `alphas` are keyed by the attractor IDs in use while the sweep is
    running. [`global_continuation`](@ref) relabel the IDs to consecutive integers once
    the sweep is over. So the same attractor may be keyed differently in the history.
    The history is meant for computing statistics over the boxes, such as where and when
    the basins changed.
"""
mutable struct BayesianUpdateSampler{D, G} <: InitialConditionsSampler
    boxes::Vector{HRectangle{Float64, SVector{D, Float64}}}
    generators::Vector{G}
    alphas::Vector{Dict{Int, Float64}}
    etas::Vector{Float64}
    sparse_n::Int
    dense_n::Int
    λ::Float64
    β::Float64
    global_reset::Bool
    boxes_flags::Vector{Bool}       # true => this box wants a dense re-sample
    layout::Vector{Pair{Int, Int}}  # (box index => n ics) of the last `generate_ics`
    step_counts::Vector{Dict{Int, Int}}
    did_reset::Bool                 # a global reset was triggered at this parameter
    history::Bool
    history_alphas::Vector{Vector{Dict{Int, Float64}}}
    history_etas::Vector{Vector{Float64}}
    history_resets::Vector{Bool}
end

function BayesianUpdateSampler(region, n_tiles::Int;
        sparse_n::Int, dense_n::Int = sparse_n^2, λ::Real = 0.7, β::Real = 0.5,
        global_reset::Bool = false, seed = abs(rand(Int)), history::Bool = false,
    )
    sparse_n ≥ 1 || throw(ArgumentError("`sparse_n` must be ≥ 1, got $sparse_n"))
    dense_n ≥ 1 || throw(ArgumentError("`dense_n` must be ≥ 1, got $dense_n"))
    0 < λ ≤ 1 || throw(ArgumentError("`λ` must be in (0, 1], got $λ"))
    β > 0 || throw(ArgumentError("`β` must be > 0, got $β"))

    boxes = _tile_region(_to_hrectangle(region), n_tiles)
    # each box gets its own generator, so that `statespace_sampler` does the actual
    # point picking and we only decide how many points come from where
    rng = Xoshiro(seed)
    generators = [statespace_sampler(box, abs(rand(rng, Int)))[1] for box in boxes]
    n = length(boxes)
    return BayesianUpdateSampler(
        boxes, generators,
        [Dict{Int, Float64}() for _ in 1:n], zeros(n),
        sparse_n, dense_n, Float64(λ), Float64(β), global_reset,
        # every box starts flagged
        fill(true, n), Pair{Int, Int}[], [Dict{Int, Int}() for _ in 1:n],
        false, history, Vector{Dict{Int, Float64}}[], Vector{Float64}[], Bool[],
    )
end

_to_hrectangle(r::HRectangle) = r
_to_hrectangle(r) = HRectangle(SVector(minimum.(r)), SVector(maximum.(r)))

# Split `region` into `n_tiles` parts per dimension, in the column major order of
# `CartesianIndices`.
function _tile_region(region::HRectangle, n_tiles::Int)
    n_tiles ≥ 1 || throw(ArgumentError("`n_tiles` must be ≥ 1, got $n_tiles"))
    mins, maxs = region.mins, region.maxs
    all(mins .< maxs) || throw(ArgumentError("`region` must have `mins .< maxs`"))
    D = length(mins)
    edges = ntuple(d -> range(Float64(mins[d]), Float64(maxs[d]); length = n_tiles + 1), D)
    boxes = HRectangle{Float64, SVector{D, Float64}}[]
    for idx in CartesianIndices(ntuple(_ -> n_tiles, D))
        lo = SVector{D, Float64}(ntuple(d -> edges[d][idx[d]], D))
        hi = SVector{D, Float64}(ntuple(d -> edges[d][idx[d] + 1], D))
        push!(boxes, HRectangle(lo, hi))
    end
    return boxes
end

n_boxes(s::BayesianUpdateSampler) = length(s.boxes)

# Number of initial conditions the *next* `generate_ics` call will produce.
function Base.length(s::BayesianUpdateSampler)
    resampling_required(s) && return s.dense_n * count(s.boxes_flags)
    return s.sparse_n * n_boxes(s)
end

function Base.show(io::IO, s::BayesianUpdateSampler{D}) where {D}
    println(io, "BayesianUpdateSampler in $(D)D")
    println(io, "  boxes:    ", n_boxes(s))
    println(io, "  sparse_n: ", s.sparse_n)
    println(io, "  dense_n:  ", s.dense_n)
    println(io, "  λ, β:     ", s.λ, ", ", s.β)
    println(io, "  global reset: ", s.global_reset ? "on" : "off")
    print(io,   "  history:  ",
          s.history ? "$(length(s.history_etas)) parameter(s) recorded" : "not kept")
end

resampling_required(s::BayesianUpdateSampler) = any(s.boxes_flags)

function generate_ics(s::BayesianUpdateSampler{D}, args...) where {D}
    resample = resampling_required(s)
    # a round that is not a re-sample is the first one of a new parameter
    resample || (foreach(empty!, s.step_counts); s.did_reset = false)
    empty!(s.layout)
    total = 0 
    if resample
        for i in eachindex(s.boxes)
            s.boxes_flags[i] || continue
            push!(s.layout, i => s.dense_n)
            total += s.dense_n
        end
    else
        for i in eachindex(s.boxes)
            push!(s.layout, i => s.sparse_n)
            total += s.sparse_n
        end
    end
    ics = Vector{SVector{D, Float64}}(undef, total)
    j = 0
    for (i, n) in s.layout
        gen = s.generators[i]
        for _ in 1:n
            ics[j += 1] = SVector{D, Float64}(gen())
        end
    end
    return StateSpaceSet(ics)
end


"""
    update_sampler!(sampler::BayesianUpdateSampler, labels)
"""
function update_sampler!(s::BayesianUpdateSampler, labels, args...)
    expected = isempty(s.layout) ? 0 : sum(last, s.layout)
    length(labels) == expected || throw(DimensionMismatch(
        "got $(length(labels)) labels but the recorded layout expects $expected; " *
        "`update_sampler!` must be called once per `generate_ics` call"
    ))
    # if it fires, every box is about to be re-learned, so the priors are left alone below
    # and this round contributes nothing but its counts
    global_reset = _global_reset_detection!(s, labels)
    cursor = 1
    for (i, n) in s.layout
        counts = _count_labels(view(labels, cursor:(cursor + n - 1)))
        cursor += n
        mergewith!(+, s.step_counts[i], counts)
        global_reset && continue
        if s.boxes_flags[i]
            # Dense round: relearn this box from scratch, no test.
            s.alphas[i] = Dict{Int, Float64}(k => c + s.β for (k, c) in counts)
            s.boxes_flags[i] = false
        else
            α = Dict{Int, Float64}(k => s.λ * v for (k, v) in s.alphas[i]) # decay priors
            η = s.etas[i] = _log_bayes_factor(counts, α, s.β)
            if η < 0 # alarm: ask for a dense re-sample
                s.boxes_flags[i] = true
            else # posterior update
                for (k, c) in counts
                    α[k] = get(α, k, s.β) + c
                end
                s.alphas[i] = α
            end
        end
    end
    s.history && !resampling_required(s) && _push_history!(s)
    return nothing
end

# An attractor born or dead changes the set of labels, which no sampling noise can fake.
# If the keyword global_reset is true, then all the boxes are flagged when attractors 
# changes from one parameter slice to the next.
function _global_reset_detection!(s::BayesianUpdateSampler, labels)
    (s.global_reset && !any(s.boxes_flags)) || return false
    prior_labels = mapreduce(keys, union!, s.alphas; init = Set{Int}())
    Set{Int}(Int(l) for l in labels) == prior_labels && return false
    fill!(s.boxes_flags, true)
    fill!(s.etas, 0.0)      # no box was tested at this parameter
    s.did_reset = true
    return true
end

function _push_history!(s::BayesianUpdateSampler)
    push!(s.history_alphas, [copy(α) for α in s.alphas])
    push!(s.history_etas, copy(s.etas))
    push!(s.history_resets, s.did_reset)
    return nothing
end

"""
    sampler_history(sampler::BayesianUpdateSampler) → NamedTuple

Return the history of alphas and etas per box, and of the global resets (`resets[i]` is
`true` if the whole tiling was re-learned at parameter `i` because the label set changed),
provided the sampler was created with `history = true`. Mind that the `alphas` are keyed by the attractor IDs of the running
sweep, which need not be those of the continuation output; see the "History" section of
[`BayesianUpdateSampler`](@ref).
"""
sampler_history(s::BayesianUpdateSampler) =
    (; alphas = s.history_alphas, etas = s.history_etas, resets = s.history_resets)

"""
    weighted_fractions(sampler::BayesianUpdateSampler, counts)

The boxes all have the same volume, so the fraction of the region belonging to basin `k`
is the average over the boxes of the fraction of each box belonging to it.

The function takes into account the boxes that have been resampled.

Note that the argument `counts` is not used in the function, the internal counts
for each box are used but the argument is necessary for the function signature. 
"""
function weighted_fractions(s::BayesianUpdateSampler, counts)
    fs = Dict{Int, Float64}()
    sampled = count(!isempty, s.step_counts)
    sampled == 0 && return fs
    for c in s.step_counts
        isempty(c) && continue
        nᵢ = sum(values(c))
        for (k, v) in c
            fs[k] = get(fs, k, 0.0) + v / (nᵢ * sampled)
        end
    end
    return fs
end

function _count_labels(labels)
    counts = Dict{Int, Int}()
    for l in labels
        k = Int(l)
        counts[k] = get(counts, k, 0) + 1
    end
    return counts
end

# Log Bayes factor η comparing the evidence for the counts under the historical prior
# `α` against the evidence under an uninformative one (every label has weight `β`).
# Both are Dirichlet-multinomial log marginal likelihoods,
# `lnΓ(α₀) - lnΓ(N + α₀) + Σₖ [lnΓ(cₖ + αₖ) - lnΓ(αₖ)]`.
function _log_bayes_factor(counts::Dict{Int, Int}, α::Dict{Int, Float64}, β::Real)
    L_hist = log_marginal_likelihood(counts, α, β)

    # Log-evidence under reset (uninformative) prior: α_i = β for all categories
    all_keys = union(keys(counts), keys(α))
    α_reset = Dict{Int, Float64}(k => β for k in all_keys)
    L_reset = log_marginal_likelihood(counts, α_reset, β)
    return L_hist - L_reset
end

"""
    log_marginal_likelihood(counts, alpha, β)

Computes the log-marginal likelihood (log-evidence) of observed counts `c`
under a Dirichlet-Multinomial model with prior `alpha`.

    L(α) = lnΓ(α₀) − lnΓ(N_s + α₀) + Σᵢ [lnΓ(cᵢ + αᵢ) − lnΓ(αᵢ)]

where α₀ = Σ αᵢ and N_s = Σ cᵢ. Categories not present in `alpha` get
the base prior `β`.
"""
function log_marginal_likelihood(new_counts::Dict{Int, Int}, alpha::Dict{Int, Float64}, β::Float64)
    all_keys = union(keys(new_counts), keys(alpha))

    alpha_0 = 0.0
    N_s = sum(values(new_counts))
    log_lik = 0.0

    for k in all_keys
        a_k = get(alpha, k, β)
        c_k = get(new_counts, k, 0)
        alpha_0 += a_k
        log_lik += loggamma(c_k + a_k) - loggamma(a_k)
    end

    log_lik += loggamma(alpha_0) - loggamma(N_s + alpha_0)
    return log_lik
end

