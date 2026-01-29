using LinearAlgebra
using Optim
using Distributions

# helper type for not computing any weights for uniform distribution
struct EverywhereUniform end

Distributions.pdf(::EverywhereUniform, u) = one(eltype(u))

"""
    StabilityMeasuresAccumulator(mapper::AttractorMapper; kwargs...)

A special data structure that allows mapping initial conditions to attractors
while _at the same time_ calculating many stability measures in the most efficient
way possible. `mapper` is any instance of an [`AttractorMapper`](@ref)
that implements the `id = mapper(u0)` syntax. This functionality was developed
as part of [Morr2025](@cite). If you use it, cite this paper along with the Attractors.jl
publication [Datseris2023](@cite).

`StabilityMeasuresAccumulator` can be used as any `AttractorMapper` with library functions
such as [`basins_fractions`](@ref). After mapping all initial conditions to attractors,
the [`finalize_accumulator`](@ref) function should be called which will return a dictionary
of all stability measures estimated by the accumulator.
Each dictionary maps the stability measure description (`String`) to a dictionary
mapping attractor IDs to the stability measure value.
Calling `reset_mapper!(accumulator)` cleans up all accumulated measures.

**Using with [`global_continuation`](@ref)**:
Since `StabilityMeasuresAccumulator` is formally an `AttractorMapper`, it can be
used with [`global_continuation`](@ref). Simply give it as a `mapper` input
to [`AttractorSeedContinueMatch`](@ref) and then call `global_continuation` as normal.
The only difference now is that `global_continuation` will not return just one
measure of stability (the basin fraction). Rather,
now the first return argument of `global_continuation` will be a
`measures_cont`, a dictionary mapping stability measures (strings)
to vectors of dictionaries. Each vector of dictionaries is similar to `fractions_cont`
of the typical [`global_continuation`](@ref): each dictionary maps attractor ID
to the corresponding nonlocal stability measure.

Use [`stability_measures_along_continuation`](@ref) for continuation of stability  measures computed
on the basis of an `AttractorsViaProximity` mapper from already found attractors.
This is useful to do for measures related to the convergence time, which is defined
more rirogously and is estimated more accurately for a proximity mapper.

## Keyword arguments

* `finite_time = 1.0`: Finite time horizon considered for the computation of the
  `finite_time_basin_stability`. Initial coditions with a convergence time
  larger than `finite_time` are not considered to be in the respective finite time basin.
  Convergence time is determined by the `mapper`.
* `weighting_distribution::Distribution`: Distribution of uncertain initial conditions
  used for example in the computation of `basin_stability`. By default it is a uniform
  distribution everywhere in the state space.\
* `distance = Centroid()`: How to compute the distance between an initial condition `u0`
  and an attractor `A`. Estimated via `set_distance([u0], A, distance)`.

## Description

`StabilityMeasuresAccumulator` efficiently uses a single `id = mapper(u0)` call
to accumulate information for many different stability measures corresponding
to each attractor of the dynamical system.
It accumulates all these different measures when different initial conditions
are mapped through it. After enough `u0`s have been given to the accumulator, they
can be finalized (comput maxima or averages) using `finalize!(accumulator)`.

You can extent this functionality by adding new stability measures as long as their
estimation can be done on the basis of the three quantities accumulated:
the basin dictionary mapping initial conditions to attractors,
the distance of each initial condition to each attractor,
and the convergence time of each initial condition.

The following stability measures are estimated for each attractor
(and the returned dictionary maps strings with the names of the measures
to the dictionaries containing the measure values for each attractor):

### Local (fixed point) stability measures

These measures apply only to fixed point attractors.
Their value is `NaN` if an attractor is not a fixed point (`length(A) > 1`).
If an unstable fixed point attractor is recorded (due to an initial condition starting
there for example), a value `Inf` is assigned to all measures.
Currently linear measures for discrete time systems are not computed.

* `characteristic_return_time`: The reciprocal of the largest real part of the
  eigenvalues of the Jacobian matrix at the fixed point.
* `reactivity`: The largest growth rate of the linearized system at the fixed point.
  See also [Krakovska2024ResilienceDynamicalSystems](@cite).
* `maximal_amplification`: The maximal (with respect to disturbances) amplification of the
  linearized system at the attractor over all time.
* `maximal_amplification_time`: The time at which the maximal amplification
  occurs.

### Nonlocal stability measures

The information for these nonlocal stability measures is accumulated while initial
conditions are mapped to attractors. Afterwards it is aggregated according to the
probability density `weighting_distribution` when calling `finalize_accumulator!`.

The word "distance" here refers to the distance established by the `distance` keyword.

* `mean_convergence_time`: The convergence time is determined by the
  `mapper` using [`convergence_time`](@ref). The mean is computed with respect
  to the `weighting_distribution`.
* `maximal_convergence_time`: The maximal convergence time of initial conditions
  to the attractor. Only initial conditions with non-zero probability under
  `weighting_distribution` are considered.
* `median_convergence_time`: The median convergence time of initial conditions under the
  `weighting_distribution`.
* `mean_convergence_pace`: The mean convergence pace of initial conditions to
  the attractor. Similar to the mean convergence time, except that each
  convergence time is divided by the distance of the respective initial
  condition to the attractor.
* `maximal_convergence_pace`: The maximal convergence pace of initial conditions
  to the attractor. Only initial conditions with non-zero probability under
  `weighting_distribution` are considered.
* `median_convergence_pace`: The median convergence pace of initial conditions under the
  `weighting_distribution`.
* `minimal_critical_shock_magnitude`: The minimal distance of the attractor to the
  closest non-zero probability point (under `weighting_distribution`) in a basin of
  attraction of a different attractor. If only a single attractor exists,
  the value `Inf` is assigned.
* `maximal_noncritical_shock_magnitude`: The distance of the attractor to the
  furthest non-zero probability point (under `weighting_distribution`) of its own basin of
  attraction. If only a single attractor exists, the value `Inf` is assigned.
* `mean_noncritical_shock_magnitude`: same as above but computing the mean under
  `weighting_distribution` instead of maximum distance.
* `basin_fraction`: The fraction of initial conditions that converge to the
  attractor.
* `basin_stability`: The fraction of initial conditions that converge to the
  attractor, weighted by `weighting_distribution`. For the default value of
  `weighting_distribution` this is identical to `basin_fraction`.
* `finite_time_basin_stability`: The fraction of initial conditions that
  converge to the attractor within the time horizon `finite_time`, weighted by
  `weighting_distribution`.
"""
mutable struct StabilityMeasuresAccumulator{AM <: AttractorMapper, V <: AbstractVector, F, M, W} <: AttractorMapper
    mapper::AM
    u0s::Vector{V}
    bs::Vector{Int}
    cts::Vector{Float64}
    finite_time::F
    weighting_distribution::W
    distance::M
end

function StabilityMeasuresAccumulator(
        mapper::AttractorMapper;
        finite_time = 1.0, weighting_distribution = EverywhereUniform(), distance = Centroid()
    )
    reset_mapper!(mapper)
    ds = referenced_dynamical_system(mapper)
    AM = typeof(mapper)
    V = typeof(current_state(ds))
    F = typeof(finite_time)
    M = typeof(distance)
    W = typeof(weighting_distribution)
    return StabilityMeasuresAccumulator{AM, V, F, M, W}(
        mapper,
        Vector{V}(),
        Vector{Int}(),
        Vector{Float64}(),
        finite_time,
        weighting_distribution,
        distance
    )
end

# Function to reset the accumulator
function reset_mapper!(a::StabilityMeasuresAccumulator)
    reset_mapper!(a.mapper)
    ds = referenced_dynamical_system(a.mapper)
    V = typeof(current_state(ds))
    a.u0s = Vector{V}()
    a.bs = Vector{Int}()
    a.cts = Vector{Float64}()
    return
end

# extensions
function extract_attractors(accumulator::StabilityMeasuresAccumulator)
    return extract_attractors(accumulator.mapper)
end
function referenced_dynamical_system(accumulator::StabilityMeasuresAccumulator)
    return referenced_dynamical_system(accumulator.mapper)
end
function convergence_time(accumulator::StabilityMeasuresAccumulator)
    return convergence_time(accumulator.mapper)
end

# Function to update the accumulator with a new state
function (accumulator::StabilityMeasuresAccumulator)(u0; show_progress = false)
    id = accumulator.mapper(u0)
    push!(accumulator.u0s, u0)
    push!(accumulator.bs, id)
    push!(accumulator.cts, convergence_time(accumulator.mapper))
    return id
end

"""
    finalize_accumulator(accumulator::StabilityMeasuresAccumulator)

Return a dictionary mapping stability measures (strings) to dictionaries
mapping attractor IDs to corresponding measure values.
See [`StabilityMeasuresAccumulator`](@ref) for more.
"""
function finalize_accumulator(accumulator::StabilityMeasuresAccumulator)
    ds = referenced_dynamical_system(accumulator)
    attractors = extract_attractors(accumulator.mapper)
    ids = vcat(collect(keys(attractors)), -1)
    js = 1:length(ids)
    ids_to_js = Dict(id => j for (j, id) in enumerate(ids))
    u0s = accumulator.u0s
    bs = accumulator.bs
    cts = accumulator.cts
    N = length(u0s)
    N == 0 && error("No initial conditions have been processed. Cannot finalize accumulator.")

    d = zeros(length(u0s), length(js))
    for i in 1:length(u0s)
        for j in js
            if ids[j] == -1
                d[i, j] = Inf
            else
                d[i, j] = set_distance(StateSpaceSet([u0s[i]]), attractors[ids[j]], accumulator.distance)
            end
        end
    end

    if (isa(accumulator.mapper, AttractorsViaProximity) && accumulator.mapper.ε != nothing)
        ε = accumulator.mapper.ε
    else
        ε = 0.0
    end

    cps = zeros(length(cts))
    for i in 1:length(cts)
        j = ids_to_js[bs[i]]
        if d[i, j] > ε
            cps[i] = cts[i] / d[i, j]
        else
            cts[i] = 0.0
            cps[i] = 0.0
        end
    end

    basin_frac = Dict(id => 0.0 for id in ids)
    basin_stab = Dict(id => 0.0 for id in ids)
    finite_time_basin_stab = Dict(id => 0.0 for id in ids)
    mean_convergence_time = Dict(id => (id == -1 ? NaN : 0.0) for id in ids)
    mean_convergence_pace = Dict(id => (id == -1 ? NaN : 0.0) for id in ids)
    maximal_convergence_time = Dict(id => (id == -1 ? NaN : 0.0) for id in ids)
    maximal_convergence_pace = Dict(id => (id == -1 ? NaN : 0.0) for id in ids)
    median_convergence_time = Dict(id => (id == -1 ? NaN : 0.0) for id in ids)
    median_convergence_pace = Dict(id => (id == -1 ? NaN : 0.0) for id in ids)
    mean_noncritical_shock_magnitude = Dict(id => (id == -1 ? NaN : 0.0) for id in ids)
    maximal_noncritical_shock_magnitude = Dict(id => (id == -1 ? NaN : 0.0) for id in ids)

    ws = [pdf(accumulator.weighting_distribution, u0) for u0 in u0s]

    for i in 1:length(u0s)
        id = bs[i]
        j = ids_to_js[id]
        w = ws[i]
        ct = cts[i]
        cp = cps[i]

        basin_frac[id] += 1 / N
        basin_stab[id] += w / N
        if ct <= accumulator.finite_time
            finite_time_basin_stab[id] += w / N
        end

        mean_convergence_time[id] += w * ct / N

        mean_convergence_pace[id] += w * cp / N

        if w > 0.0
            maximal_convergence_time[id] = max(
                maximal_convergence_time[id], ct
            )

            maximal_noncritical_shock_magnitude[id] = max(
                maximal_noncritical_shock_magnitude[id], d[i, j]
            )

            maximal_convergence_pace[id] = max(
                maximal_convergence_pace[id], cp
            )
        end

        mean_noncritical_shock_magnitude[id] += w * d[i, j] / N
    end

    normalization = sum(values(basin_stab))
    if normalization == 0.0
        normalization = 1.0
    end

    for id in ids
        basin_stab[id] /= normalization
        mean_convergence_time[id] /= normalization
        mean_convergence_pace[id] /= normalization
        mean_noncritical_shock_magnitude[id] /= normalization

        cts_id = cts[bs .== id]
        cps_id = cps[bs .== id]
        ws_id = ws[bs .== id]
        median_convergence_time[id] = weighted_median(cts_id, ws_id)
        median_convergence_pace[id] = weighted_median(cps_id, ws_id)
    end

    minimal_critical_shock_magnitude = Dict(
        id => minimum(
                (
                    d[i, ids_to_js[id]] for i in eachindex(accumulator.bs)
                    if accumulator.bs[i] != id && ws[i] > 0
                );
                init = Inf,
            )
            for id in ids
    )
    minimal_critical_shock_magnitude[-1] = NaN # no critical shock for -1 attractor

    # linear measures
    characteristic_return_time = Dict(id => NaN for id in ids)
    reactivity = Dict(id => NaN for id in ids)
    maximal_amplification = Dict(id => NaN for id in ids)
    maximal_amplification_time = Dict(id => NaN for id in ids)
    if !isdiscretetime(ds)
        jac = jacobian(ds)
        for (id, A) in attractors
            if length(A) > 1
                continue
            end
            # Get the Jacobian matrix at the fixed point
            if isinplace(ds)
                # For in-place systems, pre-allocate J and then compute it
                J = Array{Float64}(undef, length(A[1]), length(A[1]))
                jac(J, Array(A[1]), current_parameters(ds), 0)
            else
                # For out-of-place systems, compute J directly
                J = jac(Array(A[1]), current_parameters(ds), 0)
            end

            λ = min(0, maximum(real.(eigvals(J))))
            characteristic_return_time[id] = abs(1 / λ)
            H = (J + J') / 2
            evs = real.(eigvals(H))
            reactivity[id] = maximum(evs)
            if λ == 0
                maximal_amplification[id] = Inf
                maximal_amplification_time[id] = Inf
            else
                f(t) = -opnorm(exp(t * J))
                T = range(0.0, 10 * characteristic_return_time[id], length = 20001)
                step_length = T[2] - T[1]
                t0 = T[argmin(f.(T))]
                if t0 == 10 * characteristic_return_time[id] # maximum is at the end
                    res = Optim.optimize(f, t0, t0 + 100 * characteristic_return_time[id], Brent())
                else
                    res = Optim.optimize(f, max(0.0, t0 - step_length), t0 + step_length, Brent())
                end
                maximal_amplification[id] = (-1) * Optim.minimum(res)
                maximal_amplification_time[id] = Optim.minimizer(res)[1]
            end
        end
    end

    return Dict(
        "characteristic_return_time" => characteristic_return_time,
        "reactivity" => reactivity,
        "maximal_amplification" => maximal_amplification,
        "maximal_amplification_time" => maximal_amplification_time,
        "mean_convergence_time" => mean_convergence_time,
        "maximal_convergence_time" => maximal_convergence_time,
        "median_convergence_time" => median_convergence_time,
        "mean_convergence_pace" => mean_convergence_pace,
        "maximal_convergence_pace" => maximal_convergence_pace,
        "median_convergence_pace" => median_convergence_pace,
        "minimal_critical_shock_magnitude" => minimal_critical_shock_magnitude,
        "mean_noncritical_shock_magnitude" => mean_noncritical_shock_magnitude,
        "maximal_noncritical_shock_magnitude" => maximal_noncritical_shock_magnitude,
        "basin_fraction" => basin_frac,
        "basin_stability" => basin_stab,
        "finite_time_basin_stability" => finite_time_basin_stab,
    )
end

# Weighted median: smallest x with cumulative weight ≥ 0.5.
function weighted_median(
        vals::AbstractVector{<:Real},
        w::AbstractVector{<:Real}
    )
    @assert length(vals) == length(w)
    if isempty(vals)
        return NaN
    end
    s = sum(w)
    if s == 0.0
        return NaN
    end
    wn = w ./ s
    p = sortperm(vals)
    cum = 0.0
    for j in p
        cum += wn[j]
        if cum >= 0.5
            return float(vals[j])
        end
    end
    return float(vals[p[end]])
end
