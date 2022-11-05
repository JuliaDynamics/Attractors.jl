export AttractorsViaFeaturizing, group_features, extract_features
using ProgressMeter

#####################################################################################
# Structs and documentation string
#####################################################################################
abstract type GroupingConfig end

struct AttractorsViaFeaturizing{I, G<:GroupingConfig, T, F} <: AttractorMapper
    integ::I
    featurizer::F
    group_config::G
    Ttr::T
    Δt::T
    total::T
    threaded::Bool
end

"""
    AttractorsViaFeaturizing(
        ds::GeneralizedDynamicalSystem, featurizer::Function,
        grouping_config = GroupViaClustering(); kwargs...
    )

Initialize a `mapper` that maps initial conditions to attractors using a featurizing and
grouping approach. This is a supercase of the featurizing and clustering approach that
is utilized by bSTAB[^Stender2021] and MCBB[^Gelbrecht2021].
See [`AttractorMapper`](@ref) for how to use the `mapper`.

`featurizer` is a function that takes as an input an integrated trajectory `A::Dataset` and
the corresponding time vector `t` and returns a `Vector{<:Real}` of features describing the
trajectory.

See also the intermediate functions [`extract_features`](@ref) and [`group_features`](@ref),
which can be utilized when wanting to work directly with features.

## Keyword arguments
* `T=100, Ttr=100, Δt=1, diffeq=NamedTuple()`: Propagated to [`trajectory`](@ref).
* `treaded = true`: Whether to run the generation of features over threads.

## Description
The trajectory `X` of each initial condition is transformed into a vector of features. Each
feature is a number useful in _characterizing the attractor_ the initial condition ends up
at, and distinguishing it from other attractors. Example features are the mean or standard
deviation of some the dimensions of the trajectory, the entropy of some of the
dimensions, the fractal dimension of `X`, or anything else you may fancy.

The vectors of features are then grouped using one of the sevaral available grouping
configurations. Each group is assumed to be a unique attractor, and hence each initial
condition is labelled according to the group it is part of. The method thus relies on the
user having at least some basic idea about what attractors to expect in order to pick the
right features, and the right way to group them,
in contrast to [`AttractorsViaRecurrences`](@ref).

The following configuration structs can be used to decide how the features are grouped:
- [`GroupViaClustering`](@ref)
- [`GroupViaNearestFeature`](@ref)
- [`GroupViaHistogram`](@ref)

Generally speaking, the [`AttractorsViaProximity`](@ref) mapper is superior to this.
However, if the dynamical system has extremely high-dimensionality, there may be reasons
to use `AttractorsViaFeaturizing` algorithm instead, as it projects the trajectories
into a much lower dimensional representation of features.

[^Stender2021]:
    Stender & Hoffmann, [bSTAB: an open-source software for computing the basin
    stability of multi-stable dynamical systems](https://doi.org/10.1007/s11071-021-06786-5)
"""
function AttractorsViaFeaturizing(ds::GeneralizedDynamicalSystem, featurizer::Function,
    group_config::GroupingConfig = ClusteringGrouping(); T=100, Ttr=100, Δt=1,
    diffeq = NamedTuple(), threaded = false)
    if ds isa ContinuousDynamicalSystem
        T, Ttr, Δt = float.((T, Ttr, Δt))
    end
    # We put an integrator into the mapper. For parallelization, this is deepcopied.
    # Notice that if `ds` is already an integrator, both `integrator` and `trajectory`
    # still work as expected.
    return AttractorsViaFeaturizing(
        integrator(ds; diffeq), featurizer, group_config, Ttr, Δt, T, threaded
    )
end

"""
    group_features(features, group_config) → labels
Group the given features vector according to the configuration and return
the labels (vector of equal length as `features`).
See [`AttractorsViaFeaturizing`](@ref) for possible configurations.
"""
function group_features end

DynamicalSystemsBase.get_rule_for_print(m::AttractorsViaFeaturizing) =
get_rule_for_print(m.integ)

function Base.show(io::IO, mapper::AttractorsViaFeaturizing)
    ps = generic_mapper_print(io, mapper)
    println(io, rpad(" type: ", ps), nameof(typeof(mapper.integ)))
    println(io, rpad(" Ttr: ", ps), mapper.Ttr)
    println(io, rpad(" Δt: ", ps), mapper.Δt)
    println(io, rpad(" T: ", ps), mapper.total)
    return
end

ValidICS = Union{AbstractDataset, Function}

#####################################################################################
# Extension of `AttractorMapper` API
#####################################################################################
# We only extend the general `basins_fractions`, because the clustering method
# cannot map individual initial conditions to attractors
function basins_fractions(mapper::AttractorsViaFeaturizing, ics::ValidICS;
        show_progress = true, N = 1000
    )
    features = extract_features(mapper, ics; show_progress, N)
    cluster_labels  = group_features(features, mapper.group_config)
    fs = basins_fractions(cluster_labels) # Vanilla fractions method with Array input
    if typeof(ics) <: AbstractDataset
        attractors = extract_attractors(mapper, cluster_labels, ics)
        return fs, cluster_labels, attractors
    else
        return fs
    end
end

#####################################################################################
# featurizing and grouping source code
#####################################################################################
# TODO: This functionality should be a generic parallel evolving function...

"""
    extract_features(mapper, ics; N = 1000, show_progress = true)
Return a vector of the features of each initial condition in `ics` (as in
[`basins_fractions`](@ref)), using the configuration of `mapper::AttractorsViaFeaturizin`.
Keyword `N` is ignored if `ics isa Dataset`.
"""
function extract_features(mapper::AttractorsViaFeaturizing, args...; kwargs...)
    if !(mapper.threaded)
        extract_features_single(mapper, args...; kwargs...)
    else
        extract_features_threaded(mapper, args...; kwargs...)
    end
end

function extract_features_single(mapper, ics; show_progress = true, N = 1000)
    N = (typeof(ics) <: Function)  ? N : size(ics, 1) # number of actual ICs
    feature_vector = Vector{Vector{Float64}}(undef, N)
    progress = ProgressMeter.Progress(N; desc = "Integrating trajectories:", enabled=show_progress)
    for i ∈ 1:N
        ic = _get_ic(ics,i)
        feature_vector[i] = extract_feature(mapper.integ, ic, mapper)
        ProgressMeter.next!(progress)
    end
    return feature_vector
end

# TODO: We need an alternative to deep copying integrators that efficiently
# initialzes integrators for any given kind of system. But that can be done
# later in the DynamicalSystems.jl 3.0 rework.
function extract_features_threaded(mapper, ics; show_progress = true, N = 1000)
    N = (typeof(ics) <: Function)  ? N : size(ics, 1) # number of actual ICs
    integs = [deepcopy(mapper.integ) for i in 1:(Threads.nthreads() - 1)]
    pushfirst!(integs, mapper.integ)
    feature_vector = Vector{Vector{Float64}}(undef, N)
    progress = ProgressMeter.Progress(N; desc = "Integrating trajectories:", enabled=show_progress)
    Threads.@threads for i ∈ 1:N
        integ = integs[Threads.threadid()]
        ic = _get_ic(ics, i)
        feature_vector[i] = extract_feature(integ, ic, mapper)
        ProgressMeter.next!(progress)
    end
    return feature_vector
end

function extract_feature(integ, u0::AbstractVector{<:Real}, mapper)
    # Notice that this uses the low-level interface of `trajectory` that works
    # given an integrator. In DynamicalSystems 3.0 this will a part of the API
    A = trajectory(integ, mapper.total, u0; Ttr = mapper.Ttr, Δt = mapper.Δt)
    t = (mapper.Ttr):(mapper.Δt):(mapper.total+mapper.Ttr)
    return mapper.featurizer(A, t)
end

function extract_attractors(mapper::AttractorsViaFeaturizing, labels, ics)
    uidxs = unique(i -> labels[i], 1:length(labels))
    return Dict(labels[i] => trajectory(mapper.integ, mapper.total, ics[i];
    Ttr = mapper.Ttr, Δt = mapper.Δt) for i in uidxs if i ≠ -1)
end