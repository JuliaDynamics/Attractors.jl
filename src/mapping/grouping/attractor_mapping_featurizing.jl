export AttractorsViaFeaturizing, extract_features

# Flexible mapping of initial conditions into "attractors" by featurizing
# and grouping using arbitrary grouping configurations.

#####################################################################################
# Structs and documentation string
#####################################################################################
include("grouping/all_grouping_configs.jl")

struct AttractorsViaFeaturizing{DS<:DynamicalSystem, G<:GroupingConfig, F, T, D, V} <: AttractorMapper
    ds::DS
    featurizer::F
    group_config::G
    Ttr::T
    Δt::T
    total::T
    threaded::Bool
    attractors::Dict{Int, StateSpaceSet{D, V}}
end

"""
    AttractorsViaFeaturizing(
        ds::DynamicalSystem, featurizer::Function,
        grouping_config = GroupViaClustering(); kwargs...
    )

Initialize a `mapper` that maps initial conditions to attractors using a featurizing and
grouping approach. This is a supercase of the featurizing and clustering approach that
is utilized by bSTAB [Stender2021](@cite) and MCBB [Gelbrecht2021](@cite).
See [`AttractorMapper`](@ref) for how to use the `mapper`.
This `mapper` also allows the syntax `mapper(u0)` but only if the `grouping_config`
is _not_ `GroupViaClustering`.

`featurizer` is a function `f(A, t) that takes as an input an integrated trajectory
`A::StateSpaceSet` and the corresponding time vector `t` and returns a vector
`v` of features describing the trajectory.
For better performance, it is strongly recommended that `v  isa SVector{<:Real}`.

`grouping_config` is an instance of any subtype of [`GroupingConfig`](@ref) and decides
how features will be grouped into attractors, see below.

See also the intermediate functions [`extract_features`](@ref) and [`group_features`](@ref),
which can be utilized when wanting to work directly with features.

## Keyword arguments

* `T=100, Ttr=100, Δt=1`: Propagated to [`trajectory`](@ref).
* `threaded = true`: Whether to run the generation of features over threads by integrating
  trajectories in parallel.

## Description

The trajectory `X` of an initial condition is transformed into features. Each
feature is a number useful in _characterizing the attractor_ the initial condition ends up
at, and distinguishing it from other attractors. Example features are the mean or standard
deviation of some the dimensions of the trajectory, the entropy of some of the
dimensions, the fractal dimension of `X`, or anything else you may fancy.

All feature vectors (each initial condition = 1 vector) are then grouped using one of the
sevaral available grouping configurations. Each group is assumed to be a unique attractor,
and hence each initial condition is labelled according to the group it is part of.
The method thus relies on the user having at least some basic idea about what attractors
to expect in order to pick the right features, and the right way to group them,
in contrast to [`AttractorsViaRecurrences`](@ref).
"""
function AttractorsViaFeaturizing(ds::DynamicalSystem, featurizer::Function,
        group_config::GroupingConfig = GroupViaClustering();
        T=100, Ttr=100, Δt=1, threaded = true,
    )
    D = dimension(ds)
    V = eltype(current_state(ds))
    # For parallelization, the dynamical system is deepcopied.
    return AttractorsViaFeaturizing(
        ds, featurizer, group_config, Ttr, Δt, T, threaded, Dict{Int, StateSpaceSet{D,V}}(),
    )
end

DynamicalSystemsBase.rulestring(m::AttractorsViaFeaturizing) = DynamicalSystemsBase.rulestring(m.ds)

function Base.show(io::IO, mapper::AttractorsViaFeaturizing)
    ps = generic_mapper_print(io, mapper)
    println(io, rpad(" type: ", ps), nameof(typeof(mapper.ds)))
    println(io, rpad(" Ttr: ", ps), mapper.Ttr)
    println(io, rpad(" Δt: ", ps), mapper.Δt)
    println(io, rpad(" T: ", ps), mapper.total)
    println(io, rpad(" group via: ", ps), nameof(typeof(mapper.group_config)))
    return
end

ValidICS = Union{AbstractStateSpaceSet, Function}

#####################################################################################
# Extension of `AttractorMapper` API
#####################################################################################
# We only extend the general `basins_fractions`, because the clustering method
# cannot map individual initial conditions to attractors
function basins_fractions(mapper::AttractorsViaFeaturizing, ics::ValidICS;
        show_progress = true, N = 1000
    )
    features = extract_features(mapper, ics; show_progress, N)
    group_labels = group_features(features, mapper.group_config)
    fs = basins_fractions(group_labels) # Vanilla fractions method with Array input
    if typeof(ics) <: AbstractStateSpaceSet
        # TODO: If we could somehow extract the used initial conditions from `ics` 7
        # in case `ics` was a function, that would be cool...
        attractors = extract_attractors(mapper, group_labels, ics)
        overwrite_dict!(mapper.attractors, attractors)
        return fs, group_labels
    else
        return fs
    end
end

#####################################################################################
# featurizing and grouping source code
#####################################################################################
import ProgressMeter

# TODO: This functionality should be a generic parallel evolving function...
"""
    extract_features(mapper, ics; N = 1000, show_progress = true)

Return a vector of the features of each initial condition in `ics` (as in
[`basins_fractions`](@ref)), using the configuration of `mapper::AttractorsViaFeaturizing`.
Keyword `N` is ignored if `ics isa StateSpaceSet`.
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
    first_feature = extract_feature(mapper.ds, _get_ic(ics, 1), mapper)
    feature_vector = Vector{typeof(first_feature)}(undef, N)
    feature_vector[1] = first_feature
    progress = ProgressMeter.Progress(N; desc = "Integrating trajectories:", enabled=show_progress)
    for i ∈ 1:N
        ic = _get_ic(ics,i)
        feature_vector[i] = extract_feature(mapper.ds, ic, mapper)
        ProgressMeter.next!(progress)
    end
    return feature_vector
end

# TODO: We need an alternative to deep copying integrators that efficiently
# initialzes integrators for any given kind of system. But that can be done
# later in the DynamicalSystems.jl 3.0 rework.
function extract_features_threaded(mapper, ics; show_progress = true, N = 1000)
    N = (typeof(ics) <: Function)  ? N : size(ics, 1) # number of actual ICs
    systems = [deepcopy(mapper.ds) for _ in 1:(Threads.nthreads() - 1)]
    pushfirst!(systems, mapper.ds)
    first_feature = extract_feature(mapper.ds, _get_ic(ics, 1), mapper)
    feature_vector = Vector{typeof(first_feature)}(undef, N)
    feature_vector[1] = first_feature
    progress = ProgressMeter.Progress(N; desc = "Integrating trajectories:", enabled=show_progress)
    Threads.@threads for i ∈ 2:N
        ds = systems[Threads.threadid()]
        ic = _get_ic(ics, i)
        feature_vector[i] = extract_feature(ds, ic, mapper)
        ProgressMeter.next!(progress)
    end
    return feature_vector
end

function extract_feature(ds::DynamicalSystem, u0::AbstractVector{<:Real}, mapper)
    A, t = trajectory(ds, mapper.total, u0; Ttr = mapper.Ttr, Δt = mapper.Δt)
    return mapper.featurizer(A, t)
end

function extract_attractors(mapper::AttractorsViaFeaturizing, labels, ics)
    uidxs = unique(i -> labels[i], eachindex(labels))
    return Dict(labels[i] => trajectory(mapper.ds, mapper.total, ics[i];
    Ttr = mapper.Ttr, Δt = mapper.Δt)[1] for i in uidxs if i ≠ -1)
end

extract_attractors(mapper::AttractorsViaFeaturizing) = mapper.attractors