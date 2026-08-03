export InitialConditionSampler
export RandomICSampler, PerParameterICs, PerParameterInitialConditions
export BayesianUpdateSampler, sampler_history

using Random: Xoshiro
using SpecialFunctions: loggamma

"""
    InitialConditionSampler

Data structure deciding how to sample initial conditions during
[`global_continuation`](@ref).
Concerete subtypes are:

- [`RandomICSampler`]
- [`PrescribedICs`]
- [`PerParameterICs`](@ref)

`InitialConditionSampler` defines a currently experimental extendable interface
based on the internal functions
`generate_ics, update_sampler!, resampling_required, weighted_fractions`.
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
resampling_required(sampler) = false

"""
    weighted_fractions(sampler::InitialConditionSampler, counts) → Dict{Int, Float64}

Basin fractions of the region the `sampler` covers at the parameter whose sampling rounds
have just finished. `counts` maps attractor ID to the number of initial conditions that
landed on it, pooled over every round of that parameter. This is what
[`global_continuation`](@ref) can see from the outside, and it weighs every initial
condition equally; its total is the number of initial conditions the parameter used, so
nothing else about the size of the sample needs to be passed around.

The default implementation is therefore just `counts ./ sum(counts)`, which is right
whenever the sampler covers the region uniformly and in one round. A sampler that samples
some parts of the region more densely than others, or that revisits a sub-region in a
re-sampling round, must define its own method and **ignore `counts`**: only the sampler
knows the density its points were drawn from, and no estimate of the fractions is unbiased
without it. It should keep whatever running tally it needs while the rounds go by — see
[`BayesianUpdateSampler`](@ref) for an example — rather than have the continuation hold on
to the labels of every initial condition, which is the one thing here that grows with the
sample size.
"""
function weighted_fractions(sampler, counts)
    n = sum(values(counts); init = 0.0)
    return Dict{Int, Float64}(k => c / n for (k, c) in counts)
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
(elements of `pcurve` which are always a dictionary), and outputs
`N` initial conditions.
The sampler generates overall `N` initial conditions.
"""
struct PerParameterInitialConditions{F} <: InitialConditionSampler
    f::F
    N::Int
end
generate_ics(p::PerParameterInitialConditions, params, args...) = p.f(params, p.N)

"""
    BayesianUpdateSampler(region, n_tiles::Int; sparse_n, kwargs...) <: InitialConditionSampler

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
"""
mutable struct BayesianUpdateSampler{D, G} <: InitialConditionSampler
    # geometry: one box, and one point generator for it, per tile
    boxes::Vector{HRectangle{Float64, SVector{D, Float64}}}
    generators::Vector{G}
    # working memory: Dirichlet pseudo-counts and last log Bayes factor, per box
    alphas::Vector{Dict{Int, Float64}}
    etas::Vector{Float64}
    # configuration
    sparse_n::Int
    dense_n::Int
    λ::Float64
    β::Float64
    # per-round bookkeeping
    boxes_flags::Vector{Bool}       # true => this box wants a dense re-sample
    layout::Vector{Pair{Int, Int}}  # (box index => n ics) of the last `generate_ics`
    # per-parameter bookkeeping: label counts per box, over every round of the current
    # parameter. This is all `weighted_fractions` needs, and unlike a record of the
    # labels themselves its size is bounded by the number of boxes times the number of
    # attractors, no matter how many initial conditions are drawn.
    step_counts::Vector{Dict{Int, Int}}
    # optional record of `alphas` and `etas` as they were at the end of each parameter,
    # both empty unless the sampler was built with `history = true`
    history::Bool
    history_alphas::Vector{Vector{Dict{Int, Float64}}}
    history_etas::Vector{Vector{Float64}}
end

function BayesianUpdateSampler(region, n_tiles::Int;
        sparse_n::Int, dense_n::Int = sparse_n^2, λ::Real = 0.7, β::Real = 0.5,
        seed = abs(rand(Int)), history::Bool = false,
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
        sparse_n, dense_n, Float64(λ), Float64(β),
        # every box starts flagged: the first round is then just an ordinary
        # re-sampling round, which is exactly the dense initialisation we want
        fill(true, n), Pair{Int, Int}[], [Dict{Int, Int}() for _ in 1:n],
        history, Vector{Dict{Int, Float64}}[], Vector{Float64}[],
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
    print(io,   "  history:  ",
          s.history ? "$(length(s.history_etas)) parameter(s) recorded" : "not kept")
end

resampling_required(s::BayesianUpdateSampler) = any(s.boxes_flags)

function generate_ics(s::BayesianUpdateSampler, args...)
    empty!(s.layout)
    # A re-sampling round only visits the boxes that asked for it, densely; every other
    # round visits all of them, sparsely.
    resample = resampling_required(s)
    # ... so a round that is not a re-sample is the first one of a new parameter, and the
    # tally the fractions are computed from starts afresh
    resample || foreach(empty!, s.step_counts)
    for i in eachindex(s.boxes)
        n = s.boxes_flags[i] ? s.dense_n : (resample ? 0 : s.sparse_n)
        n == 0 || push!(s.layout, i => n)
    end
    # `layout` tells `update_sampler!` how to slice the labels between boxes it gets back
    ics = Vector{Vector{Float64}}(undef, length(s))
    j = 1
    for (i, n) in s.layout
        gen = s.generators[i]
        for _ in 1:n
            ics[j] = copy(gen()) # the generator reuses its output buffer
            j += 1
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
    cursor = 1
    for (i, n) in s.layout
        counts = _count_labels(view(labels, cursor:(cursor + n - 1)))
        cursor += n
        # every round of this parameter contributes to the box's tally, whether it ends
        # up flagged or not; `weighted_fractions` reads it once the rounds are over
        mergewith!(+, s.step_counts[i], counts)
        if s.boxes_flags[i]
            # Dense round: relearn this box from scratch, no test.
            s.alphas[i] = Dict{Int, Float64}(k => c + s.β for (k, c) in counts)
            s.boxes_flags[i] = false
        else
            α = Dict{Int, Float64}(k => s.λ * v for (k, v) in s.alphas[i]) # decay priors
            η = s.etas[i] = _log_bayes_factor(counts, α, s.β)
            if η < 0 # the prior cannot explain the data: ask for a dense re-sample
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

function _push_history!(s::BayesianUpdateSampler)
    push!(s.history_alphas, [copy(α) for α in s.alphas])
    push!(s.history_etas, copy(s.etas))
    return nothing
end

"""
    sampler_history(sampler::BayesianUpdateSampler) → NamedTuple
Return the history of alphas and etas per box if the keyword argument history is true
"""
sampler_history(s::BayesianUpdateSampler) =
    (; alphas = s.history_alphas, etas = s.history_etas)

"""
    weighted_fractions(sampler::BayesianUpdateSampler, counts)

The boxes all have the same volume, so the fraction of the region belonging to basin `k`
is the average over the boxes of the fraction of each box belonging to it. Read off
`step_counts`, the per-box tally `update_sampler!` keeps as the rounds of the parameter go
by, so that nothing proportional to the number of initial conditions is ever stored.

The pooled `counts` given by [`global_continuation`](@ref) are ignored, as they weigh each
box by how densely it happened to be covered: a box that panicked and was re-sampled would
count `dense_n/sparse_n` times more than a quiet one, which bends the fractions of the
whole region towards wherever the basins are currently changing.
"""
function weighted_fractions(s::BayesianUpdateSampler, counts)
    fs = Dict{Int, Float64}()
    # a box that received nothing cannot contribute an estimate, and must not count
    # towards the average either, or the fractions would not sum to one
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
    N = sum(values(counts))
    α₀ = 0.0; nk = 0
    L_hist = 0.0; L_reset = 0.0
    for k in union(keys(counts), keys(α))
        nk += 1
        a = get(α, k, β)
        c = get(counts, k, 0)
        α₀ += a
        L_hist += loggamma(c + a) - loggamma(a)
        L_reset += loggamma(c + β) - loggamma(β)
    end
    L_hist += loggamma(α₀) - loggamma(N + α₀)
    L_reset += loggamma(nk * β) - loggamma(N + nk * β)
    return L_hist - L_reset
end
