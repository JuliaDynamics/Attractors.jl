export BayesianUpdateSampler
using SpecialFunctions: loggamma, digamma, trigamma

"""
    BayesianUpdateSampler(region, tiles_n::Int, sparse_n::Int; kwargs...) <: InitialConditionsSampler

Sampler that may re-sample initial conditions wih larger density by employing the
adaptive Bayesian framework of [Haerter2026](@cite).

`region` is anything [`statespace_sampler`](@ref) accepts as a region: an `HRectangle`,
or a tuple of ranges/`(min, max)` pairs, one per dimension.
This is then tiled into `tiles_n^D` equally sized boxes, with `D` the system dimension.
Each box is sampled with `sparse_n` samples by default. However,
when this sampler is used with [`global_continuation`](@ref), boxes may be resampled
with higher density according to the Bayesian framework of [Haerter2026](@cite),
see description below.

## Keyword arguments

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
  are otherwise overwritten in place. They are stored in the fields
  `history_alphas, history_etas`. Note that the IDs of these containers have not been
  matched; you can apply `match_sequentially!` to them if you want though.

## Description

We outline how this sampler works in global continuation.
At the first parameter every box is sampled with `dense_n` initial conditions and its
prior is initialised as `α_k = c_k + β` from the label counts `c` and atractor k.
At every later parameter:

1. `generate_ics` draws `sparse_n` points per box.
2. `update_sampler!` slices the returned labels per box, decays that box's prior by `λ`,
   and computes `η`, the Bayes factor. If `η < 0` the box is flagged and its prior left untouched;
   otherwise the posterior `α_k ← λα_k + c_k` is stored.
3. If any box was flagged, `resampling_required` returns `true` and the continuation loop
   calls `generate_ics` again. That round draws `dense_n` points for the flagged boxes
   *only*, and `update_sampler!` rebuilds their priors from scratch and clears the flags.

Because a dense round performs no test, at most one re-sampling round happens per
parameter. Every box starts out flagged, so the dense initialisation of the first
parameter is nothing but step 3 applied to all of them.

This sampler adds the number of alarms
in the `other` component of the continuation output. If `tiles_n > 1`,
then the basin entropy and basin entropy variance is also attributed.
Note that their accuracy depends on the number of tiles, and the basin entropy does not
equal exactly the KDTree-based estimation that occurs generically through
[`StabilityQuantifiersAccumulator`](@ref).
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
    boxes_flags::Vector{Bool} # true => this box wants a dense re-sample
    prior_alarms::Int # number of alarms raised at the current parameter
    layout::Vector{Pair{Int, Int}}  # (box index => n ics) of the last `generate_ics`
    step_counts::Vector{Dict{Int, Int}}
    did_reset::Bool                 # a global reset was triggered at this parameter
    history::Bool
    history_alphas::Vector{Vector{Dict{Int, Float64}}}
    history_etas::Vector{Vector{Float64}}
end

function BayesianUpdateSampler(region, tiles_n::Int, sparse_n::Int;
        dense_n::Int = sparse_n^2, λ::Real = 0.7, β::Real = 0.5,
        global_reset::Bool = false, seed = abs(rand(Int)), history::Bool = false,
    )
    sparse_n ≥ 1 || throw(ArgumentError("`sparse_n` must be ≥ 1, got $sparse_n"))
    dense_n ≥ 1 || throw(ArgumentError("`dense_n` must be ≥ 1, got $dense_n"))
    0 < λ ≤ 1 || throw(ArgumentError("`λ` must be in (0, 1], got $λ"))
    β > 0 || throw(ArgumentError("`β` must be > 0, got $β"))

    boxes = _tile_region(_to_hrectangle(region), tiles_n)
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
        fill(true, n), 0, Pair{Int, Int}[], [Dict{Int, Int}() for _ in 1:n],
        false, history, Vector{Dict{Int, Float64}}[], Vector{Float64}[],
    )
end

_to_hrectangle(r::HRectangle) = r
_to_hrectangle(r) = HRectangle(SVector(minimum.(r)), SVector(maximum.(r)))

# Split `region` into `tiles_n` parts per dimension
function _tile_region(region::HRectangle, tiles_n::Int)
    tiles_n ≥ 1 || throw(ArgumentError("`tiles_n` must be ≥ 1, got $tiles_n"))
    mins, maxs = region.mins, region.maxs
    all(i -> mins[i] < maxs[i], eachindex(mins)) ||
        throw(ArgumentError("`region` must have `mins .< maxs`"))
    D = length(mins)
    N = tiles_n^D # Careful in high dim, this will explode
    edges = [range(Float64(mins[d]), Float64(maxs[d]); length = tiles_n + 1) for d in 1:D]
    boxes = Vector{HRectangle{Float64, SVector{D, Float64}}}(undef, N)
    tile = ones(Int, D)  # tile index per dimension, the first one varying fastest
    lo, hi = zeros(D), zeros(D)
    for i in 1:N
        for d in 1:D
            lo[d], hi[d] = edges[d][tile[d]], edges[d][tile[d] + 1]
        end
        boxes[i] = HRectangle(SVector{D, Float64}(lo), SVector{D, Float64}(hi))
        for d in 1:D  # odometer: bump the first dimension, carry into the next ones
            tile[d] < tiles_n ? (tile[d] += 1; break) : (tile[d] = 1)
        end
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
    resample || (foreach(empty!, s.step_counts); s.did_reset = false; s.prior_alarms = 0)
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
                s.prior_alarms += 1
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
    return nothing
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

# The following fuctions are for estimating the basin entropy

# Posterior mean of the Shannon entropy of a Dirichlet(α) distribution,
# E[S] = ψ(α₀ + 1) − Σₖ (αₖ/α₀) ψ(αₖ + 1),    α₀ = Σₖ αₖ
function bayes_entropy(alpha::AbstractDict{Int, Float64})
    a0 = sum(values(alpha))
    a0 > 0 || return 0.0
    term2 = 0.0
    for val in values(alpha)
        term2 += (val / a0) * digamma(val + 1)
    end
    return digamma(a0 + 1) - term2
end

# Exact posterior variance of the Shannon entropy under Dirichlet(α), following
# Wolpert & Wolf (1995), Theorem 16 (Eqs. 16.1–16.2). The `i ≠ j` cross terms are
# summed with the `(Σ)² − Σ(·²)` trick, so the cost is O(K) rather than O(K²).
function bayes_entropy_variance(alpha::Union{AbstractVector{Float64}, AbstractDict{Int, Float64}})
    vals = isa(alpha, AbstractDict) ? collect(values(alpha)) : alpha
    a0 = sum(vals)
    a0 > 0 || return 0.0

    # Polygamma terms of the total, needed by every summand
    psi_a0_1 = digamma(a0 + 1)
    psi_a0_2 = digamma(a0 + 2)
    tri_a0_2 = trigamma(a0 + 2)

    E_S = 0.0      # E[S], Eq. 16.1
    sum_A = 0.0    # Σ Aᵢ, for the cross-term trick
    sum_A2 = 0.0   # Σ Aᵢ²
    sum_a_sq = 0.0 # Σ αᵢ²
    D = 0.0        # diagonal (i == j) terms of Eq. 16.2

    for a_i in vals
        a_i > 0 || continue
        E_S -= (a_i / a0) * (digamma(a_i + 1) - psi_a0_1)

        A_i = a_i * (digamma(a_i + 1) - psi_a0_2)
        sum_A += A_i
        sum_A2 += A_i^2
        sum_a_sq += a_i^2

        D += (a_i * (a_i + 1)) / (a0 * (a0 + 1)) *
             ((digamma(a_i + 2) - psi_a0_2)^2 + trigamma(a_i + 2) - tri_a0_2)
    end

    C = ((sum_A^2 - sum_A2) - tri_a0_2 * (a0^2 - sum_a_sq)) / (a0 * (a0 + 1))
    # Var = E[S²] − E[S]²; clamp away float noise around 0
    return max(0.0, (C + D) - E_S^2)
end


# Basin entropy of the whole region: the mean over boxes of their posterior mean entropy.
mean_basin_entropy(alphas::AbstractVector{<:AbstractDict{Int, Float64}}) =
    isempty(alphas) ? 0.0 : sum(bayes_entropy, alphas) / length(alphas)

# Variance of [`mean_entropy`](@ref). The boxes are treated as independent, so the
# variance of their mean is `(1/N²) Σᵢ Var[Sᵢ]`.
function mean_basin_entropy_variance(alphas::AbstractVector{<:AbstractDict{Int, Float64}})
    n = length(alphas)
    n > 0 || return 0.0
    return sum(bayes_entropy_variance, alphas) / n^2
end

function init_other_info(s::BayesianUpdateSampler)
    out = Dict{String, Any}("alarms" => Int)
    if length(s.boxes) > 1
        out["basin_entropy"] = Float64
        out["basin_entropy_variance"] = Float64
    end
    return out
end

function extract_other_info(s::BayesianUpdateSampler)
    out = Dict{String, Any}("alarms" => s.prior_alarms)
    if length(s.boxes) > 1
        out["basin_entropy"] = mean_basin_entropy(s.alphas)
        out["basin_entropy_variance"] = mean_basin_entropy_variance(s.alphas)
    end
    return out
end
