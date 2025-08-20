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
Each dictionary maps the measure description (`String`) to a dictionary
mapping attractor IDs to the measure value.
Calling `reset_mapper!(accumulator)` cleans up all accumulated measures.

**Using with [`global_continuation`](@ref)**:
Since `StabilityMeasuresAccumulator` is formally an `AttractorMapper`, it can be
used with [`global_continuation`](@ref). Simply give it as a `mapper` input
to [`AttractorSeedContinueMatch`](@ref) and then call `global_continuation` as normal.
The only difference now is that `global_continuation` will not return just one
measure of nonlocal stability (the basin fraction). Rather,
now the first return argument of `global_continuation` will be a
`measures_cont`, a dictionary mapping nonlocal stability measures (strings)
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
  distribution everywhere in the state space.

## Description

`StabilityMeasuresAccumulator` efficiently uses a single `id = mapper(u0)` call
to accumulate information for many differnt stability measures corresponding
to each attractor of the dynamical system.
It accumulates all these different measures when different initial conditions
are mapped through it. After enough `u0`s have been given to the accumulator, they
can be finalized (comput maxima or averages) using `finalize!(accumulator)`.

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

These nonlocal stability measures are accumulated while initial conditions are mapped
to attractors. Afterwards they are averaged according to the probability density
`weighting_distribution` when calling `finalize_accumulator!`. The word "distance" here
refers to the distance established by the `metric` keyword.

* `mean_convergence_time`: The convergence time is determined by the
  `mapper` using [`convergence_time`](@ref). The mean is computed with respect
  to the `weighting_distribution`.
* `maximal_convergence_time`: The maximal convergence time of initial conditions
  to the attractor. Only initial conditions with non-zero probability under
  `weighting_distribution` are considered.
* `median_convergence_time`: The median convergence time of initial conditions.
  Only initial conditions with non-zero probability under `weighting_distribution`
  are considered.
* `mean_convergence_pace`: The mean convergence pace of initial conditions to
  the attractor. Similar to the mean convergence time, except that each
  convergence time is divided by the distance of the respective initial
  condition to the attractor.
* `maximal_convergence_pace`: The maximal convergence pace of initial conditions
  to the attractor. Only initial conditions with non-zero probability under
  `weighting_distribution` are considered.
* `median_convergence_pace`: The median convergence pace of initial conditions.
  Only initial conditions with non-zero probability under `weighting_distribution`
  are considered.
* `minimal_critical_shock_magnitude`: The minimal distance of the attractor to the
  closest non-zero probability point (under `weighting_distribution`) in a basin of
  attraction of a different attractor. If only a single attractor exists,
  the value `Inf` is assigned.
* `maximal_noncritical_shock_magnitude`: The distance of the attractor to the
  furthest non-zero probability point (under `weighting_distribution`) of its own basin of
  attraction. If only a single attractor exists, the value `Inf` is assigned.
* `mean_noncritical_shock_magnitude`: same as above but the mean instead of maximum distance.
* `basin_fraction`: The fraction of initial conditions that converge to the
  attractor.
* `basin_stability`: The fraction of initial conditions that converge to the
  attractor, weighted by `weighting_distribution`. For the default value of
  `weighting_distribution` this is identical to `basin_fraction`.
* `finite_time_basin_stability`: The fraction of initial conditions that
  converge to the attractor within the time horizon `finite_time`, weighted by
  `weighting_distribution`.
"""
mutable struct StabilityMeasuresAccumulator{AM<:AttractorMapper, Dims, Datatype} <: AttractorMapper
    mapper::AM
    basin_points::Dict{Int, StateSpaceSet{Dims, Datatype}}
    finite_time_basin_points::Dict{Int, StateSpaceSet{Dims, Datatype}}
    nonzero_measure_basin_points::Dict{Int, StateSpaceSet{Dims, Datatype}}
    mean_convergence_time::Dict{Int, Float64}
    maximal_convergence_time::Dict{Int, Float64}
    mean_convergence_pace::Dict{Int, Float64}
    maximal_convergence_pace::Dict{Int, Float64}
    minimal_critical_shock_magnitude::Dict{Int, Float64}
    maximal_noncritical_shock_magnitude::Dict{Int, Float64}
    mean_noncritical_shock_magnitude::Dict{Int, Float64}
    convergence_times::Dict{Int, Vector{Float64}}
    convergence_paces::Dict{Int, Vector{Float64}}
    convergence_weights::Dict{Int, Vector{Float64}}
    finite_time::Float64 # Discussed that this should be F and metric::M but this leads to errors...
    weighting_distribution::Union{EverywhereUniform, Distribution}
    metric::Metric
end

function StabilityMeasuresAccumulator(mapper::AttractorMapper;
        finite_time=1.0, weighting_distribution=EverywhereUniform(), metric=Euclidean()
    )
    reset_mapper!(mapper)
    ds = referenced_dynamical_system(mapper)
    Dims = dimension(ds)
    Datatype = eltype(current_state(ds))
    AM = typeof(mapper)
    StabilityMeasuresAccumulator{AM, Dims, Datatype}(
        mapper,
        Dict{Int, StateSpaceSet{Dims, Datatype}}(),
        Dict{Int, StateSpaceSet{Dims, Datatype}}(),
        Dict{Int, StateSpaceSet{Dims, Datatype}}(),
        Dict{Int, Float64}(),
        Dict{Int, Float64}(),
        Dict{Int, Float64}(),
        Dict{Int, Float64}(),
        Dict{Int, Float64}(),
        Dict{Int, Float64}(),
        Dict{Int, Float64}(),
        Dict{Int, Vector{Float64}}(),
        Dict{Int, Vector{Float64}}(),
        Dict{Int, Vector{Float64}}(),
        finite_time,
        weighting_distribution,
        metric
    )
end

# Function to reset the accumulator
function reset_mapper!(a::StabilityMeasuresAccumulator)
    reset_mapper!(a.mapper)
    ds = referenced_dynamical_system(a.mapper)
    Dims = dimension(ds)
    Datatype = eltype(current_state(ds))
    a.basin_points = Dict{Int, StateSpaceSet{Dims,Datatype}}()
    a.finite_time_basin_points = Dict{Int, StateSpaceSet{Dims,Datatype}}()
    a.nonzero_measure_basin_points = Dict{Int, StateSpaceSet{Dims,Datatype}}()
    a.mean_convergence_time = Dict{Int, Float64}()
    a.maximal_convergence_time = Dict{Int, Float64}()
    a.mean_convergence_pace = Dict{Int, Float64}()
    a.maximal_convergence_pace = Dict{Int, Float64}()
    a.convergence_times = Dict{Int, Vector{Float64}}()
    a.convergence_paces = Dict{Int, Vector{Float64}}()
    a.convergence_weights = Dict{Int, Vector{Float64}}()
    a.minimal_critical_shock_magnitude = Dict{Int, Float64}()
    a.maximal_noncritical_shock_magnitude = Dict{Int, Float64}()
    a.mean_noncritical_shock_magnitude = Dict{Int, Float64}()
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
    if (isa(accumulator.mapper, AttractorsViaProximity) && accumulator.mapper.ε != nothing)
      ε = accumulator.mapper.ε
    else
      ε = 0.0
    end
    # Update basin points
    if !(id in keys(accumulator.basin_points))
        accumulator.basin_points[id] = StateSpaceSet([u0])
    else
        push!(accumulator.basin_points[id], u0)
    end

    # Gather convergence time and pace
    ct = convergence_time(accumulator.mapper)
    attractors = extract_attractors(accumulator.mapper)
    if id == -1
      u0_dist = Inf
    else
      u0_dist = set_distance(
        StateSpaceSet([u0]), attractors[id],
        StateSpaceSets.StrictlyMinimumDistance(true, accumulator.metric)
      )
    end
    ct = u0_dist > ε ? ct : 0

    this_weight = pdf(accumulator.weighting_distribution, u0)
    has_nonzero_weight = this_weight > 0.0
    if !(id in keys(accumulator.convergence_times))
        accumulator.convergence_times[id] = Float64[]
        accumulator.convergence_paces[id] = Float64[]
        accumulator.convergence_weights[id] = Float64[]
    end
    has_nonzero_weight && push!(accumulator.convergence_times[id], ct)
    has_nonzero_weight && push!(accumulator.convergence_paces[id], ct/u0_dist)
    has_nonzero_weight && push!(accumulator.convergence_weights[id], this_weight)

    if !(id in keys(accumulator.mean_convergence_time)) # initialize times/paces
        accumulator.mean_convergence_time[id] = 0.0
        accumulator.maximal_convergence_time[id] = 0.0
        accumulator.mean_convergence_pace[id] = 0.0
        accumulator.maximal_convergence_pace[id] = 0.0
    end
    accumulator.mean_convergence_time[id] += this_weight*ct
    accumulator.mean_convergence_pace[id] += this_weight*ct/u0_dist
    if has_nonzero_weight
        accumulator.maximal_convergence_time[id] = max(
            accumulator.maximal_convergence_time[id], ct
        )
        accumulator.maximal_convergence_pace[id] = max(
            accumulator.maximal_convergence_pace[id], ct/u0_dist
        )
    end

    # Update finite time basin points
    if ct <= accumulator.finite_time
        if !(id in keys(accumulator.finite_time_basin_points))
            accumulator.finite_time_basin_points[id] = StateSpaceSet([u0])

        else
            push!(accumulator.finite_time_basin_points[id], u0)
        end
    end

    # Update nonzero measure basin points
    if pdf(accumulator.weighting_distribution, u0) > 0.0
        if !(id in keys(accumulator.nonzero_measure_basin_points))
            accumulator.nonzero_measure_basin_points[id] = StateSpaceSet([u0])
        else
            push!(accumulator.nonzero_measure_basin_points[id], u0)
        end

        for id2 in keys(attractors) 
          if id2 != id
            thisdist = set_distance(
                StateSpaceSet([u0]), attractors[id2],
                StateSpaceSets.StrictlyMinimumDistance(true, accumulator.metric)
            )
            if !(id2 in keys(accumulator.minimal_critical_shock_magnitude))
              accumulator.minimal_critical_shock_magnitude[id2] = thisdist
            else
              accumulator.minimal_critical_shock_magnitude[id2] = min(
              accumulator.minimal_critical_shock_magnitude[id2], thisdist
              )
            end
          end
        end

        if !(id in keys(accumulator.maximal_noncritical_shock_magnitude))
            accumulator.maximal_noncritical_shock_magnitude[id] = u0_dist
        else
            accumulator.maximal_noncritical_shock_magnitude[id] = max(
                accumulator.maximal_noncritical_shock_magnitude[id], u0_dist
            )
        end

        if !(id in keys(accumulator.mean_noncritical_shock_magnitude))
            accumulator.mean_noncritical_shock_magnitude[id] = u0_dist*pdf(accumulator.weighting_distribution, u0)
        else
            accumulator.mean_noncritical_shock_magnitude[id] = (
                accumulator.mean_noncritical_shock_magnitude[id] +
                u0_dist*pdf(accumulator.weighting_distribution, u0)
            )
        end

        
    end

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
    Dims = dimension(ds)
    Datatype = eltype(current_state(ds))
    attractors = extract_attractors(accumulator.mapper)
    bpoints = accumulator.basin_points
    foreach(key -> !(key in keys(bpoints)) && (bpoints[key] = StateSpaceSet{Dims, Datatype}()), keys(attractors))
    foreach(key -> !(key in keys(accumulator.finite_time_basin_points)) && (accumulator.finite_time_basin_points[key] = StateSpaceSet{Dims, Datatype}()), keys(attractors))
    foreach(key -> !(key in keys(accumulator.nonzero_measure_basin_points)) && (accumulator.nonzero_measure_basin_points[key] = StateSpaceSet{Dims, Datatype}()), keys(attractors))
    normalization = sum(sum(pdf(accumulator.weighting_distribution, point) for point in v) for v in values(bpoints) if !isempty(v))
    N = sum(length(v) for v in values(bpoints))

    # basin stability and finite time basin stability
    basin_frac = Dict(k => length(v)/N for (k, v) in bpoints)
    basin_stab = Dict(k => (isempty(v) ? 0.0 : sum(pdf(accumulator.weighting_distribution, point) for point in v)) for (k, v) in bpoints)
    finite_time_basin_stab = Dict(k => (isempty(v) ? 0.0 : sum(pdf(accumulator.weighting_distribution, point) for point in v)) for (k, v) in accumulator.finite_time_basin_points)
    normalization != 0.0 && (foreach(k -> basin_stab[k] /= normalization, keys(basin_stab)); foreach(k -> finite_time_basin_stab[k] /= normalization, keys(finite_time_basin_stab)))
    foreach(k -> !(k in keys(basin_stab)) && (basin_stab[k] = 0.0), keys(attractors))
    foreach(k -> !(k in keys(finite_time_basin_stab)) && (finite_time_basin_stab[k] = 0.0), keys(attractors))

    # normalize all mean values
    for key in keys(accumulator.mean_convergence_time)
        if normalization != 0.0 && basin_stab[key] != 0.0
            accumulator.mean_convergence_time[key] /= (normalization*basin_stab[key])
            accumulator.mean_convergence_pace[key] /= (normalization*basin_stab[key])
            accumulator.mean_noncritical_shock_magnitude[key] /= (normalization*basin_stab[key])
        end
    end
    foreach(key -> !(key in keys(accumulator.mean_convergence_time)) && (accumulator.mean_convergence_time[key] = Inf), keys(attractors))
    foreach(key -> !(key in keys(accumulator.maximal_convergence_time)) && (accumulator.maximal_convergence_time[key] = Inf), keys(attractors))
    foreach(key -> !(key in keys(accumulator.mean_convergence_pace)) && (accumulator.mean_convergence_pace[key] = Inf), keys(attractors))
    foreach(key -> !(key in keys(accumulator.maximal_convergence_pace)) && (accumulator.maximal_convergence_pace[key] = Inf), keys(attractors))
    foreach(key -> !(key in keys(accumulator.mean_noncritical_shock_magnitude)) && (accumulator.mean_noncritical_shock_magnitude[key] = Inf), keys(attractors))
    foreach(key -> !(key in keys(accumulator.maximal_noncritical_shock_magnitude)) && (accumulator.maximal_noncritical_shock_magnitude[key] = Inf), keys(attractors))
    foreach(key -> !(key in keys(accumulator.minimal_critical_shock_magnitude)) && (accumulator.minimal_critical_shock_magnitude[key] = Inf), keys(attractors))

    # median convergence time and pace
    median_convergence_time = Dict{Int, Float64}()
    median_convergence_pace = Dict{Int, Float64}()
    for key in keys(accumulator.convergence_times)
        w = get(accumulator.convergence_weights, key, Float64[])
        median_convergence_time[key] = weighted_median(accumulator.convergence_times[key], w)
    end

    for key in keys(accumulator.convergence_paces)
        w = get(accumulator.convergence_weights, key, Float64[])
        median_convergence_pace[key] = weighted_median(accumulator.convergence_paces[key], w)
    end

    # linear measures
    characteristic_return_time = Dict(k => NaN for k in keys(attractors))
    reactivity = Dict(k => NaN for k in keys(attractors))
    maximal_amplification = Dict(k => NaN for k in keys(attractors))
    maximal_amplification_time = Dict(k => NaN for k in keys(attractors))
    if !isdiscretetime(ds)
        jac = jacobian(ds)
        for (id, A) in attractors
            if length(A) > 1
                characteristic_return_time[id] = NaN
                reactivity[id] = NaN
                maximal_amplification[id] = NaN
                maximal_amplification_time[id] = NaN
                continue
            end
            J = isinplace(ds) ? Array{Float64}(undef, length(A[1]), length(A[1])) : jac(Array(A[1]), initial_parameters(ds), 0)
            isinplace(ds) && jac(J, Array(A[1]), initial_parameters(ds), 0)

            thisλ = min(0, maximum(real.(eigvals(J))))
            characteristic_return_time[id] = abs(1 / thisλ)
            H = (J + J') / 2
            evs = real.(eigvals(H))
            reactivity[id] = maximum(evs)
            if thisλ == 0
                maximal_amplification[id] = Inf
                maximal_amplification_time[id] = Inf
            else
                #res = Optim.optimize(t -> (-1) * opnorm(exp(t[1] * J)), [0.0], [Inf], [1.0])
                #res = Optim.optimize(t -> (-1) * opnorm(exp(t * J)), 0.0, 1000.0, Brent())

                f(t) = -opnorm(exp(t*J))        # same objective
                T = range(0.0, 1000*characteristic_return_time[id], length=2001)  # coarse grid every 0.5
                t0 = argmin(f.(T))                 # best coarse point
                res = Optim.optimize(f, max(0.0, t0-5), min(1000.0, t0+5), Brent())
                maximal_amplification[id] = exp((-1) * Optim.minimum(res))
                maximal_amplification_time[id] = Optim.minimizer(res)[1]
            end
        end
    end

    return Dict(
        "characteristic_return_time" => characteristic_return_time,
        "reactivity" => reactivity,
        "maximal_amplification" => maximal_amplification,
        "maximal_amplification_time" => maximal_amplification_time,
        "mean_convergence_time" => accumulator.mean_convergence_time,
        "maximal_convergence_time" => accumulator.maximal_convergence_time,
        "median_convergence_time" => median_convergence_time,
        "mean_convergence_pace" => accumulator.mean_convergence_pace,
        "maximal_convergence_pace" => accumulator.maximal_convergence_pace,
        "median_convergence_pace" => median_convergence_pace,
        "minimal_critical_shock_magnitude" => accumulator.minimal_critical_shock_magnitude,
        "mean_noncritical_shock_magnitude" => accumulator.mean_noncritical_shock_magnitude,
        "maximal_noncritical_shock_magnitude" => accumulator.maximal_noncritical_shock_magnitude,
        "basin_fraction" => basin_frac,
        "basin_stability" => basin_stab,
        "finite_time_basin_stability" => finite_time_basin_stab,
    )
end

# Weighted median: smallest x with cumulative weight ≥ 0.5.
function weighted_median(vals::AbstractVector{<:Real},
                         w::AbstractVector{<:Real})
    @assert length(vals) == length(w)
    if isempty(vals); return NaN; end
    s = sum(w)
    if s == 0.0; return NaN; end
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