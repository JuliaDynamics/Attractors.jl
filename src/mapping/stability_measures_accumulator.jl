using LinearAlgebra
using Optim
using Distributions
struct EverywhereUniform end # helper type for not computing any weights for uniform distribution

Distributions.pdf(::EverywhereUniform, u) = one(eltype(u))

"""
    StabilityMeasuresAccumulator(mapper::AttractorMapper; kwargs...)

A special data structure that allows mapping initial conditions to attractors
_while_ at the same time calculating many stability measures in the most efficient
way possible. `mapper` is any instance of an [`AttractorMapper`](@ref)
that satisfies the `id = mapper(u0)` syntax.

`StabilityMeasuresAccumulator` can be used as any `AttractorMapper` with library functions such as
[`basins_fractions`](@ref). After mapping all initial conditions to attractors,
the [`finalize!`](@ref) function should be called which will return a dictionary
of all stability measures estimated by the accumulator.
Each dictionary maps the measure description (`Symbol`) to a dictionary
mapping attractor IDs to the measure value.

## Keyword arguments

* `T = 1.0`: Finite time horizon considered for the computation of the
  `finite_time_basin_stability`. Initial coditions with a convergence time
  larger than `T` are not considered to be in the respective finite time basin.
  Convergence time is determined by the `mapper`.
* `d::Distribution`: Distribution of uncertain initial conditions
  used for example in the computation of `basin_stability`. By default it is a uniform
  distribution everywhere in the state space.
* `p = initial_parameters(referenced_dynamical_system(mapper))`: Alternative
  parameter setting of the dynamical system considered for computing the
  `persistence` of the attractors.

## Description

`StabilityMeasuresAccumulator` efficiently uses a single `id = mapper(u0)` call
to accumulate information for many differnt stability measures corresponding
to each attractor of the dynamical system.
### Finalizing the stability measures
After calling the accumulator on all initial conditions in the state space sample of the `mapper`, the
stability measures need to be finalized. To obtain the final stability measures,
call `Attractors.finalize(accumulator::StabilityMeasuresAccumulator)`. This
function returns a dictionary with the stability measures. Each value of this
dictionary is again a dictionary mapping the attractor ID to the respective
stability measure.

The following stability measures are estimated for each attractor:

### Local (fixed point) stability measures

These measures apply only to fixed point attractors.
Their value is `NaN` if an attractor is not a fixed point (`length(A) > 1`).
If an unstable fixed point attractor is recorded (due to an initial condition starting
there for example), a value `Inf` is assigned to all measures.

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
to attractors. Afterwards they are averaged according to the probability density `d`
when calling `finalize!`. The word "distance" here refers to the distance established
by the `metric` keyword.

* `mean_convergence_time`: The convergence time is determined by the
  `mapper` using [`convergence_time`](@ref). The mean is computed with respect
  to the probability distribution `d`.
* `maximal_convergence_time`: The maximal convergence time of initial conditions
  to the attractor. Only initial conditions with non-zero probability under `d` 
  are considered.
* `mean_convergence_pace`: The mean convergence pace of initial conditions to
  the attractor. Similar to the mean convergence time, except that each
  convergence time is divided by the distance of the respective initial
  condition to the attractor.
* `maximal_convergence_pace`: Same as above but maximum instead of mean.
* `minimal_critical_shock_magnitude`: The minimal distance of the attractor to the
  closest non-zero probability point (under `d`) in a basin of attraction of a 
  different attractor.
* `maximal_noncritical_shock_magnitude`: The distance of the attractor to the
  furthest non-zero probability point (under `d`) of its own basin of attraction.
* `basin_stability`: The fraction of initial conditions that converge to the
  attractor, weighted by `d`.
* `finite_time_basin_stability`: The fraction of initial conditions that
  converge to the attractor within the time horizon `T`, weighted by `d`.
* `persistence`: Trajectories from all points of the attractor are evolved under
  the alternative parameter setting `p`. The persistence is the time at which
  one of the trajectories first leaves the original basin of attraction of the
  attractor. If the trajectories do not leave the basin of attraction, the
  persistence is set to `Inf`.
"""
mutable struct StabilityMeasuresAccumulator{AM<:AttractorMapper, Dims, Datatype} <: AttractorMapper
    mapper::AM
    basin_points::Dict{Int64, StateSpaceSet{Dims, Datatype}}
    finite_time_basin_points::Dict{Int64, StateSpaceSet{Dims, Datatype}}
    nonzero_measure_basin_points::Dict{Int64, StateSpaceSet{Dims, Datatype}}
    mean_convergence_time::Dict{Int64, Float64}
    maximal_convergence_time::Dict{Int64, Float64}
    mean_convergence_pace::Dict{Int64, Float64}
    maximal_convergence_pace::Dict{Int64, Float64}
    convergence_times::Dict{Int64, Vector{Float64}}
    convergence_paces::Dict{Int64, Vector{Float64}}
    # TODO: It is wrong to say that T is a float. it should be type parameterized.
    T::Float64
    d::Union{EverywhereUniform, Distribution}
    p::AbstractArray
    metric::Metric
    function StabilityMeasuresAccumulator(mapper::AttractorMapper; T=1.0::Float64, d=EverywhereUniform(), p=initial_parameters(referenced_dynamical_system(mapper))::Vector, metric=Euclidean())
        reset_mapper!(mapper)
        ds = referenced_dynamical_system(mapper)
        Dims = dimension(ds)
        Datatype = eltype(current_state(ds))
        AM = typeof(mapper)
        new{AM,Dims,Datatype}(
            mapper,
            Dict{Int64, StateSpaceSet{Dims, Datatype}}(),
            Dict{Int64, StateSpaceSet{Dims, Datatype}}(),
            Dict{Int64, StateSpaceSet{Dims, Datatype}}(),
            Dict{Int64, Float64}(),
            Dict{Int64, Float64}(),
            Dict{Int64, Float64}(),
            Dict{Int64, Float64}(),
            Dict{Int64, Vector{Float64}}(),
            Dict{Int64, Vector{Float64}}(),
            T,
            d,
            p,
            metric
        )
    end
end

# Function to reset the accumulator
function reset_mapper!(a::StabilityMeasuresAccumulator)
    reset_mapper!(a.mapper)
    ds = referenced_dynamical_system(a.mapper)
    Dims = dimension(ds)
    Datatype = eltype(current_state(ds))
    reinit!(ds)
    a.basin_points = Dict{Int64, StateSpaceSet{Dims,Datatype}}()
    a.finite_time_basin_points = Dict{Int64, StateSpaceSet{Dims,Datatype}}()
    a.nonzero_measure_basin_points = Dict{Int64, StateSpaceSet{Dims,Datatype}}()
    a.mean_convergence_time = Dict{Int64, Float64}()
    a.maximal_convergence_time = Dict{Int64, Float64}()
    a.mean_convergence_pace = Dict{Int64, Float64}()
    a.maximal_convergence_pace = Dict{Int64, Float64}()
    a.convergence_times = Dict{Int64, Vector{Float64}}()
    a.convergence_paces = Dict{Int64, Vector{Float64}}()
end

# Function to extract attractors from the accumulator
function extract_attractors(accumulator::StabilityMeasuresAccumulator)
    return extract_attractors(accumulator.mapper)
end

# Function to get the referenced dynamical system from the accumulator
function referenced_dynamical_system(accumulator::StabilityMeasuresAccumulator)
    return referenced_dynamical_system(accumulator.mapper)
end

# Function to update the accumulator with a new state
function (accumulator::StabilityMeasuresAccumulator)(u0; show_progress = false)
    id = accumulator.mapper(u0)

    # Update basin points
    if !(id in keys(accumulator.basin_points))
        accumulator.basin_points[id] = StateSpaceSet([u0])
    else
        push!(accumulator.basin_points[id], u0)
    end

    # Gather convergence time and pace
    ct = convergence_time(accumulator.mapper)
    attractors = extract_attractors(accumulator.mapper) ### AM more elegent way to do this?
    u0_dist = id == -1 ? Inf64 : set_distance(StateSpaceSet([u0]), attractors[id], StateSpaceSets.StrictlyMinimumDistance(true, accumulator.metric))
    accumulator.convergence_times[id] = pdf(accumulator.d, u0) > 0.0 ? push!(get(accumulator.convergence_times, id, []), ct) : get(accumulator.convergence_times, id, [])
    accumulator.convergence_paces[id] = pdf(accumulator.d, u0) > 0.0 ? push!(get(accumulator.convergence_paces, id, []), ct/u0_dist) : get(accumulator.convergence_paces, id, [])

    # Update mean and maximal convergence time and pace
    accumulator.mean_convergence_time[id] = get(accumulator.mean_convergence_time, id, 0.0) + pdf(accumulator.d, u0)*ct
    accumulator.maximal_convergence_time[id] = pdf(accumulator.d, u0) > 0.0 ? max(get(accumulator.maximal_convergence_time, id, 0.0), ct) : get(accumulator.maximal_convergence_time, id, 0.0)
    accumulator.mean_convergence_pace[id] = get(accumulator.mean_convergence_pace, id, 0.0) + pdf(accumulator.d, u0)*ct/u0_dist
    accumulator.maximal_convergence_pace[id] = pdf(accumulator.d, u0) > 0.0 ? max(get(accumulator.maximal_convergence_pace, id, 0.0), ct/u0_dist) : get(accumulator.maximal_convergence_pace, id, 0.0)
    
    # Update finite time basin points
    if ct <= accumulator.T
        if !(id in keys(accumulator.finite_time_basin_points))
            accumulator.finite_time_basin_points[id] = StateSpaceSet([u0])
        else
            push!(accumulator.finite_time_basin_points[id], u0)
        end
    end

    # Update nonzero measure basin points
    if pdf(accumulator.d, u0) > 0.0
        if !(id in keys(accumulator.nonzero_measure_basin_points))
            accumulator.nonzero_measure_basin_points[id] = StateSpaceSet([u0])
        else
            push!(accumulator.nonzero_measure_basin_points[id], u0)
        end
    end

    return id
end

# Function to finalize the stability measures
function finalize(accumulator::StabilityMeasuresAccumulator)
    attractors = extract_attractors(accumulator.mapper)
    normalization = sum([sum(pdf(accumulator.d, point) for point in accumulator.basin_points[key]) for key in keys(accumulator.basin_points)])

    # Calculate mean convergence time and pace
    normalization != 0.0 && foreach(key -> accumulator.mean_convergence_time[key] /= normalization, keys(accumulator.mean_convergence_time))
    normalization != 0.0 && foreach(key -> accumulator.mean_convergence_pace[key] /= normalization, keys(accumulator.mean_convergence_pace))
    foreach(key -> !(key in keys(accumulator.mean_convergence_time)) && (accumulator.mean_convergence_time[key] = Inf64), keys(attractors))
    foreach(key -> !(key in keys(accumulator.maximal_convergence_time)) && (accumulator.maximal_convergence_time[key] = Inf64), keys(attractors))
    foreach(key -> !(key in keys(accumulator.mean_convergence_pace)) && (accumulator.mean_convergence_pace[key] = Inf64), keys(attractors))
    foreach(key -> !(key in keys(accumulator.maximal_convergence_pace)) && (accumulator.maximal_convergence_pace[key] = Inf64), keys(attractors))

    # Calculate median convergence time and pace
    median_convergence_time = Dict{Int64, Float64}()
    median_convergence_pace = Dict{Int64, Float64}()
    foreach(key -> median_convergence_time[key] = median(accumulator.convergence_times[key]), keys(accumulator.convergence_times))
    foreach(key -> median_convergence_pace[key] = median(accumulator.convergence_paces[key]), keys(accumulator.convergence_paces))

    # Calculate minimal fatal shock and maximal non-fatal shock magnitude
    minimal_fatal_shock_magnitudes = Dict{Int64, Float64}()
    maximal_nonfatal_shock_magnitudes = Dict{Int64, Float64}()
    for key1 in keys(attractors)
        minimal_fatal_shock_magnitudes[key1] = Inf64
        maximal_nonfatal_shock_magnitudes[key1] = Inf64
        (length(keys(accumulator.nonzero_measure_basin_points)) == 1 && key1 in keys(accumulator.nonzero_measure_basin_points)) && continue
        minimal_fatal_shock_magnitudes[key1] = minimum([set_distance(attractors[key1], accumulator.nonzero_measure_basin_points[key2], StateSpaceSets.StrictlyMinimumDistance(true, accumulator.metric)) for key2 in keys(accumulator.nonzero_measure_basin_points) if key1 != key2])
        maximal_nonfatal_shock_magnitudes[key1] = maximum([accumulator.metric(a, b) for a in attractors[key1] for b in accumulator.nonzero_measure_basin_points[key1]])
    end

    # Calculate persistence
    ds = referenced_dynamical_system(accumulator.mapper)
    persistence = Dict{Int64, Float64}()
    if accumulator.p == initial_parameters(ds)
        for key in keys(attractors)
            persistence[key] = Inf64
        end
    else
        p_init = initial_parameters(ds)
        set_parameters!(ds, accumulator.p)
        if typeof(accumulator.mapper)<:AttractorsViaRecurrences
            new_mapper = AttractorsViaRecurrences(ds, accumulator.mapper.grid.grid)
        elseif typeof(accumulator.mapper)<:AttractorsViaProximity
            new_mapper = AttractorsViaProximity(ds, extract_attractors(accumulator.mapper), getfield(accumulator.mapper, :ε))
        else
            new_mapper = nothing
            println("Unsupported mapper type")
        end
        if new_mapper != nothing
            for attr_key in keys(attractors)
                persistence[attr_key] = Inf64
                for u0 in attractors[attr_key]
                    id_new = new_mapper(u0)
                    ct = convergence_time(new_mapper)
                    X, t = trajectory(ds, ct, u0, Δt=0.01)
                    X_dict = Dict(zip(t, vec(X)))
                    for this_t in t
                        distances = setsofsets_distances(Dict([(attr_key, StateSpaceSet([X_dict[this_t]]))]), accumulator.basin_points, StateSpaceSets.StrictlyMinimumDistance(true, accumulator.metric))
                        #println("Distance of attractor $attr_key at time $this_t: $distances")
                        if findmin(distances[attr_key])[2] != attr_key
                            persistence[attr_key] = min(persistence[attr_key], this_t)
                            break
                        end
                    end
                end
            end
        else
            for key in keys(attractors)
                persistence[key] = Inf64
            end
        end
        set_parameters!(ds, p_init)
    end

    # Calculate basin stability and finite time basin stability
    basin_stab = Dict(key => sum(pdf(accumulator.d, point) for point in accumulator.basin_points[key]) for key in keys(accumulator.basin_points))
    finite_time_basin_stab = Dict(key => sum(pdf(accumulator.d, point) for point in accumulator.finite_time_basin_points[key]) for key in keys(accumulator.finite_time_basin_points))
    normalization != 0.0 && (foreach(key -> basin_stab[key] /= normalization, keys(basin_stab)); foreach(key -> finite_time_basin_stab[key] /= normalization, keys(finite_time_basin_stab)))
    foreach(key -> !(key in keys(basin_stab)) && (basin_stab[key] = 0.0), keys(attractors))
    foreach(key -> !(key in keys(finite_time_basin_stab)) && (finite_time_basin_stab[key] = 0.0), keys(attractors))

    # Initialize dictionaries for linear measures
    characteristic_return_time = Dict(id => NaN for id in keys(attractors))
    reactivity = Dict(id => NaN for id in keys(attractors))
    maximal_amplification = Dict(id => NaN for id in keys(attractors))
    maximal_amplification_time = Dict(id => NaN for id in keys(attractors))

    # Calculate linear measures
    ds = referenced_dynamical_system(accumulator.mapper)
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
            maximal_amplification[id] = Inf64
            maximal_amplification_time[id] = Inf64
        else
            res = Optim.optimize(t -> (-1) * opnorm(exp(t[1] * J)), [0.0], [Inf64], [1.0])
            maximal_amplification[id] = (-1) * Optim.minimum(res)
            maximal_amplification_time[id] = Optim.minimizer(res)[1]
        end
    end

    # Return the final stability measures
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
        "minimal_fatal_shock_magnitude" => minimal_fatal_shock_magnitudes,
        "maximal_nonfatal_shock_magnitude" => maximal_nonfatal_shock_magnitudes,
        "basin_stability" => basin_stab,
        "finite_time_basin_stability" => finite_time_basin_stab,
        "persistence" => persistence
    )
end
