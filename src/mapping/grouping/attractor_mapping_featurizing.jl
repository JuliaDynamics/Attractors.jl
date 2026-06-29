export BasinMapFeaturizeGroup, extract_features

# Flexible mapping of initial conditions into "attractors" by featurizing
# and grouping using arbitrary grouping configurations.

#####################################################################################
# Structs and documentation string
#####################################################################################
include("all_grouping_configs.jl")

struct BasinMapFeaturizeGroup{DS <: DynamicalSystem, G <: GroupingConfig, F, T, C, SSS <: AbstractStateSpaceSet} <: BasinMap
    ds::DS
    featurizer::F
    group_config::G
    Ttr::T
    Δt::T
    total::T
    container::C
    threaded::Bool
    attractors::Dict{Int, SSS}
end

"""
    BasinMapFeaturizeGroup(
        ds::DynamicalSystem, featurizer::Function,
        grouping_config = GroupViaClustering(); kwargs...
    )

Initialize a `bmap` that maps initial conditions to attractors using a featurizing and
grouping approach. This is a supercase of the featurizing and clustering approach that
is utilized by bSTAB [Stender2021](@cite) and MCBB [Gelbrecht2020](@cite).
See [`BasinMap`](@ref) for how to use the `bmap`.
This `bmap` also allows the syntax `bmap(u0)` if the `grouping_config` is
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

* `T=100, Ttr=100, Δt=1`: Propagated to `DynamicalSystems.trajectory` for
  integrating an initial condition to yield `A, t`.
* `container`: Also propagated to `DynamicalSystems.trajectory`, but here its default value
  is `Vector` if `dimension(ds) ≥ 32`.
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
in contrast to [`BasinMapRecurrences`](@ref).

Attractors are stored and can be accessed with [`extract_attractors`](@ref),
however it should be clear that this bmap never actually finds attractors.
They way we store attractors is by picking the first initial condition that belongs
to the corresponding "attractor group", and then recording its trajectory
with the same arguments `T, Ttr, Δt`. This is stored as the attractor,
but of course there is no guarantee that this is actually an attractor.

Note that the convergence time of this bmap is always the constant `T + Ttr`.
"""
function BasinMapFeaturizeGroup(
        ds::DynamicalSystem, featurizer::Function,
        group_config::GroupingConfig = GroupViaClustering(); threaded::Bool = true,
        T = 100, Ttr = 100, Δt = 1, container = dimension(ds) ≥ 32 ? Vector : SVector,
    )
    D = dimension(ds)
    V = eltype(current_state(ds))
    T, Ttr, Δt = promote(T, Ttr, Δt)
    # For parallelization, the dynamical system is deepcopied.
    return BasinMapFeaturizeGroup(
        ds, featurizer, group_config, Ttr, Δt, T, container, threaded, Dict{Int, StateSpaceSet{D, V}}(),
    )
end

# BasinMap API:
function reset_mapper!(bmap::BasinMapFeaturizeGroup)
    empty!(bmap.attractors)
    return
end
function (bmap::BasinMapFeaturizeGroup)(u0)
    ds = referenced_dynamical_system(bmap)
    A, t = trajectory(ds, bmap.total, u0; container = bmap.container, Ttr = bmap.Ttr, Δt = bmap.Δt)
    f = bmap.featurizer(A, t)
    id = feature_to_group(f, bmap.group_config)
    # store attractor if this is the first id
    if !haskey(bmap.attractors, id) && id > 0
        bmap.attractors[id] = A
    end
    return id
end
convergence_time(bmap::BasinMapFeaturizeGroup) = bmap.Ttr + bmap.total

DynamicalSystemsBase.rulestring(m::BasinMapFeaturizeGroup) = DynamicalSystemsBase.rulestring(m.ds)

function Base.show(io::IO, bmap::BasinMapFeaturizeGroup)
    ps = generic_mapper_print(io, bmap)
    println(io, rpad(" Ttr: ", ps), bmap.Ttr)
    println(io, rpad(" Δt: ", ps), bmap.Δt)
    println(io, rpad(" T: ", ps), bmap.total)
    println(io, rpad(" group via: ", ps), nameof(typeof(bmap.group_config)))
    println(io, rpad(" featurizer: ", ps), nameof(bmap.featurizer))
    return
end


#####################################################################################
# Extension of `BasinMap` API: basins_fractions_grouped
#####################################################################################
function basins_fractions_grouped(bmap::BasinMapFeaturizeGroup, ics, progress, labels)
    features = extract_features(bmap, ics; progress)
    glabels = group_features(features, bmap.group_config)
    if !isempty(labels)
        labels .= @view(glabels[1:length(labels)])
    end
    extract_attractors!(bmap, glabels, ics)
    fs = basins_fractions(glabels) # Vanilla fractions method with Array input
    return fs
end

#####################################################################################
# featurizing and grouping source code
#####################################################################################
import ProgressMeter

"""
    extract_features(bmap::BasinMapFeaturizeGroup, ics; N = 1000, show_progress = true)

Return a vector of the features of each initial condition in `ics` (as in
[`basins_fractions`](@ref)), using the configuration of `bmap`.
Keyword `N` is ignored if `ics isa StateSpaceSet`.
"""
function extract_features(
        bmap::BasinMapFeaturizeGroup, ics;
        show_progress = true, progress = nothing, N = 1000
    )
    if isnothing(progress)
        progress = ProgressMeter.Progress(
            N;
            desc = "Integrating trajectories:", enabled = show_progress
        )
    end
    return if !(bmap.threaded)
        extract_features_single(bmap, ics; N, progress)
    else
        extract_features_threaded(bmap, ics; N, progress)
    end
end

function extract_features_single(bmap, ics; progress, N)
    N = (typeof(ics) <: Function) ? N : size(ics, 1) # number of actual ICs
    ds = referenced_dynamical_system(bmap)
    first_feature = extract_feature(ds, _get_ic(ics, 1), bmap)
    feature_vector = Vector{typeof(first_feature)}(undef, N)
    feature_vector[1] = first_feature
    for i in 2:N
        ic = _get_ic(ics, i)
        feature_vector[i] = extract_feature(ds, ic, bmap)
        ProgressMeter.next!(progress)
    end
    return feature_vector
end

using OhMyThreads: @tasks, @local
function extract_features_threaded(bmap, ics; progress, N)
    N = (typeof(ics) <: Function) ? N : size(ics, 1) # number of actual ICs
    # systems = [deepcopy(bmap.ds) for _ in 1:(Threads.nthreads() - 1)]
    # pushfirst!(systems, bmap.ds)
    first_feature = extract_feature(bmap.ds, _get_ic(ics, 1), bmap)
    feature_vector = Vector{typeof(first_feature)}(undef, N)
    feature_vector[1] = first_feature
    @tasks for i in 2:N
        @local ds = deepcopy(referenced_dynamical_system(bmap))
        ic = _get_ic(ics, i)
        feature_vector[i] = extract_feature(ds, ic, bmap)
        ProgressMeter.next!(progress)
    end
    return feature_vector
end

function extract_feature(ds::DynamicalSystem, u0::AbstractVector{<:Real}, bmap)
    A, t = trajectory(ds, bmap.total, u0;
        container = bmap.container, Ttr = bmap.Ttr, Δt = bmap.Δt
    )
    return bmap.featurizer(A, t)
end

function extract_attractors!(bmap::BasinMapFeaturizeGroup, labels, ics)
    uidxs = unique(i -> labels[i], eachindex(labels))
    attractors = Dict(
        labels[i] => trajectory(
                bmap.ds, bmap.total, ics[i];
                Ttr = bmap.Ttr, Δt = bmap.Δt
            )[1] for i in uidxs if i > 0
    )
    overwrite_dict!(bmap.attractors, attractors)
    return
end

_extract_attractors(bmap::BasinMapFeaturizeGroup) = bmap.attractors
