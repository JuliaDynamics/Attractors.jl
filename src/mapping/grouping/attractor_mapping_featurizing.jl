export AttractorsViaFeaturizing, extract_features

# Flexible mapping of initial conditions into "attractors" by featurizing
# and grouping using arbitrary grouping configurations.

#####################################################################################
# Structs and documentation string
#####################################################################################
include("all_grouping_configs.jl")

struct AttractorsViaFeaturizing{DS <: DynamicalSystem, G <: GroupingConfig, F, T, SSS <: AbstractStateSpaceSet} <: AttractorMapper
    ds::DS
    featurizer::F
    group_config::G
    Ttr::T
    Δt::T
    total::T
    threaded::Bool
    attractors::Dict{Int, SSS}
end

"""
    AttractorsViaFeaturizing(
        ds::DynamicalSystem, featurizer::Function,
        grouping_config = GroupViaClustering(); kwargs...
    )

Initialize a `mapper` that maps initial conditions to attractors using a featurizing and
grouping approach. This is a supercase of the featurizing and clustering approach that
is utilized by bSTAB [Stender2021](@cite) and MCBB [Gelbrecht2020](@cite).
See [`AttractorMapper`](@ref) for how to use the `mapper`.
This `mapper` also allows the syntax `mapper(u0)` if the `grouping_config` is
is _not_ `GroupViaClustering` or `GroupViaPairwiseComparison`.

`featurizer` is a function `f(A, t)` that takes as an input an integrated trajectory
`A::StateSpaceSet` and the corresponding time vector `t` and returns a vector
`v` of features describing the trajectory.
For better performance, it is strongly recommended that `v isa SVector{<:Real}`.

`grouping_config` is an instance of any subtype of [`GroupingConfig`](@ref) and decides
how features will be grouped into attractors, see below.

See also the intermediate functions [`extract_features`](@ref) and [`group_features`](@ref),
which can be utilized when wanting to work directly with features.

## Keyword arguments

* `T=100, Ttr=100, Δt=1`: Propagated to `DynamicalSystems.trajectory` for integrating
  an initial condition to yield `A, t`.
* `threaded = true`: Whether to run the generation of features over threads by integrating
  trajectories in parallel.

## Description

The trajectory `X` of an initial condition is transformed into features. Each
feature is a number useful in _characterizing the attractor_ the initial condition ends up
at, and **distinguishing it from other attractors**. Example features are the mean or standard
deviation of some the dimensions of the trajectory, the entropy of some of the
dimensions, the fractal dimension of `X`, or anything else you may fancy.

All feature vectors (each initial condition = 1 feature vector) are then grouped using one of the
sevaral available grouping configurations. Each group is assumed to be a unique attractor,
and hence each initial condition is labelled according to the group it is part of.
The method thus relies on the user having at least some basic idea about what attractors
to expect in order to pick the right features, and the right way to group them,
in contrast to [`AttractorsViaRecurrences`](@ref).

Attractors are stored and can be accessed with [`extract_attractors`](@ref),
however it should be clear that this mapper never actually finds attractors.
They way we store attractors is by picking the first initial condition that belongs
to the corresponding "attractor group", and then recording its trajectory
with the same arguments `T, Ttr, Δt`. This is stored as the attractor,
but of course there is no guarantee that this is actually an attractor.
"""
function AttractorsViaFeaturizing(
        ds::DynamicalSystem, featurizer::Function,
        group_config::GroupingConfig = GroupViaClustering();
        T = 100, Ttr = 100, Δt = 1, threaded = true,
    )
    D = dimension(ds)
    V = eltype(current_state(ds))
    T, Ttr, Δt = promote(T, Ttr, Δt)
    # For parallelization, the dynamical system is deepcopied.
    return AttractorsViaFeaturizing(
        ds, featurizer, group_config, Ttr, Δt, T, threaded, Dict{Int, StateSpaceSet{D, V}}(),
    )
end

# AttractorMapper API:
function reset_mapper!(mapper::AttractorsViaFeaturizing)
    empty!(mapper.attractors)
    return
end
function (mapper::AttractorsViaFeaturizing)(u0)
    f = extract_feature(referrenced_dynamical_system(mapper), u0, mapper)
    return feature_to_group(f, mapper.group_config)
end

DynamicalSystemsBase.rulestring(m::AttractorsViaFeaturizing) = DynamicalSystemsBase.rulestring(m.ds)

function Base.show(io::IO, mapper::AttractorsViaFeaturizing)
    ps = generic_mapper_print(io, mapper)
    println(io, rpad(" Ttr: ", ps), mapper.Ttr)
    println(io, rpad(" Δt: ", ps), mapper.Δt)
    println(io, rpad(" T: ", ps), mapper.total)
    println(io, rpad(" group via: ", ps), nameof(typeof(mapper.group_config)))
    println(io, rpad(" featurizer: ", ps), nameof(mapper.featurizer))
    return
end


#####################################################################################
# Extension of `AttractorMapper` API: basins_fractions_grouped
#####################################################################################
function basins_fractions_grouped(mapper::AttractorsViaFeaturizing, ics, progress, labels)
    features = extract_features(mapper, ics; progress)
    glabels = group_features(features, mapper.group_config)
    if !isempty(labels)
        labels .= glabels
    end
    extract_attractors!(mapper, glabels, ics)
    fs = basins_fractions(glabels) # Vanilla fractions method with Array input
    return fs
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
function extract_features(
        mapper::AttractorsViaFeaturizing, ics;
        show_progress = true, progress = nothing, N = 1000
    )
    if isnothing(progress)
        progress = ProgressMeter.Progress(
            N;
            desc = "Integrating trajectories:", enabled = show_progress
        )
    end
    return if !(mapper.threaded)
        extract_features_single(mapper, ics; N, progress)
    else
        extract_features_threaded(mapper, ics; N, progress)
    end
end

function extract_features_single(mapper, ics; progress, N)
    N = (typeof(ics) <: Function) ? N : size(ics, 1) # number of actual ICs
    first_feature = extract_feature(mapper.ds, _get_ic(ics, 1), mapper)
    feature_vector = Vector{typeof(first_feature)}(undef, N)
    feature_vector[1] = first_feature
    for i in 1:N
        ic = _get_ic(ics, i)
        feature_vector[i] = extract_feature(mapper.ds, ic, mapper)
        ProgressMeter.next!(progress)
    end
    return feature_vector
end

using OhMyThreads: @tasks, @local
function extract_features_threaded(mapper, ics; progress, N)
    N = (typeof(ics) <: Function) ? N : size(ics, 1) # number of actual ICs
    # systems = [deepcopy(mapper.ds) for _ in 1:(Threads.nthreads() - 1)]
    # pushfirst!(systems, mapper.ds)
    first_feature = extract_feature(mapper.ds, _get_ic(ics, 1), mapper)
    feature_vector = Vector{typeof(first_feature)}(undef, N)
    feature_vector[1] = first_feature
    @tasks for i in 2:N
        @local ds = deepcopy(referenced_dynamical_system(mapper))
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

function extract_attractors!(mapper::AttractorsViaFeaturizing, labels, ics)
    uidxs = unique(i -> labels[i], eachindex(labels))
    attractors = Dict(
        labels[i] => trajectory(
                mapper.ds, mapper.total, ics[i];
                Ttr = mapper.Ttr, Δt = mapper.Δt
            )[1] for i in uidxs if i > 0
    )
    overwrite_dict!(mapper.attractors, attractors)
    return
end

_extract_attractors(mapper::AttractorsViaFeaturizing) = mapper.attractors
