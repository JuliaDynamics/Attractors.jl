export aggregate_continuation, aggregate_attractors, aggregate_fractions

"""
    aggregate_continuation(attractors_cont, featurizer, group_config; kw...)

Aggregate the attractors of a [`global_continuation`](@ref) result.

This is the post-processing companion of [`global_continuation`](@ref).
After a `global_continuation` with [`AttractorSeedContinueMatch`](@ref),
call this function to "group" attractors that have
similar features into a single merged attractor. The
merged attractors carry IDs that stay consistent along the parameter curve.

## Arguments

1. `attractors_cont`: a vector of dictionaries mapping IDs to attractors (`StateSpaceSet`s),
   exactly as returned by [`global_continuation`](@ref) (or [`basins_fractions`](@ref)).
2. `featurizer`: a 1-argument function mapping an attractor to a feature vector.
   Features expected by [`GroupingConfig`](@ref) are typically `SVector`s.
3. `group_config`: a subtype of [`GroupingConfig`](@ref).

All keyword arguments are propagated into [`MatchBySSSetDistance`](@ref). This matcher
is used to match feature groups along the continuation. This is possible because
a group of features is a vector of vectors, just like a StateSpaceSet is.

## Return

1. `agg_attractors_cont`: like `attractors_cont`, but mapping each group ID to the merged
   `StateSpaceSet` of its members. (merged attractors are unions of the individual ones)
2. `centroids_cont`: a vector of dictionaries mapping each group ID to the
   centroid (mean) of its members' feature vectors at that step.
3. `members_cont`: a vector of dictionaries mapping each group ID to the
   vector of original attractor IDs (the keys of `attractors_cont[i]`) merged into it at that step.

All three share consistent group IDs along the parameter axis.

## Description

At each parameter step the attractors are turned into feature vectors with `featurizer` and
partitioned into groups with `group_config` (a [`GroupingConfig`](@ref)); the members of a group
are merged into a single `StateSpaceSet`. To keep group IDs consistent along the parameter curve,
groups are matched between consecutive parameter steps *in feature space*: the set of member
feature vectors of each group is collected into a `StateSpaceSet` and these are matched with
[`MatchBySSSetDistance`](@ref) using the `distance` keyword. By default `distance = Centroid()`,
i.e. groups are matched by the distance between their feature-vector centroids; any other
set-distance accepted by [`setsofsets_distances`](@ref) (e.g. [`Hausdorff`](@ref),
[`StrictlyMinimumDistance`](@ref), or a custom `f(A, B)`) instead compares the full feature-vector
sets. Because matching is in feature space, `featurizer` must return *numeric* feature vectors
(e.g. `SVector`s).

## Aggregating stability measures

To obtain stability measures for the aggregated groups, pass
`agg_attractors_cont` to [`stability_measures_along_continuation`](@ref). Each merged group is
then treated as a single attractor, so every measure — including those that need the raw basin
data, such as medians and critical shock magnitudes — is computed correctly for the group.

See the [aggregation example](@ref aggregate_continuation_example) for an illustration.

!!! note "Aggregating basin fractions only"
    If you only care about aggregating the basin fractions, there is no reason to
    go through the route of `stability_measures_along_continuation`.
    Simply give the returned `members_cont`, together with the `fractions_cont` from the
    continuation, to [`aggregate_fractions`](@ref).
    The denser the sampling for the original continuation was, the more accurate
    the aggregated fractions will be, and there is no reason to re-run a continuation.
"""
function aggregate_continuation(attractors_cont::Vector, featurizer, group_config; kw...)
    P = length(attractors_cont)
    agg_attractors_cont = typeof(attractors_cont)()
    feature_sets_cont = Dict[]
    members_cont = Dict{Int, Vector{Int}}[]
    for i in 1:P
        agg_attractors, feature_sets, members = _aggregate_one_step(
            attractors_cont[i], featurizer, group_config
        )
        push!(agg_attractors_cont, agg_attractors)
        push!(feature_sets_cont, feature_sets)
        push!(members_cont, members)
    end
    # match groups across steps by the set distance of their feature-vector sets (i.e. in feature
    # space), then propagate the same ID remapping to the merged attractors and membership dicts
    if P > 1
        rmaps = match_sequentially!(feature_sets_cont, MatchBySSSetDistance(; kw...))
        match_sequentially!(agg_attractors_cont, rmaps)
        match_sequentially!(members_cont, rmaps)
    end
    centroids_cont = [_feature_centroids(feature_sets) for feature_sets in feature_sets_cont]
    return agg_attractors_cont, centroids_cont, members_cont
end

"""
    aggregate_attractors(attractors::Dict, featurizer, group_config) → agg_attractors, centroids, members

Aggregate the `attractors` of a *single* parameter (a dictionary mapping IDs to
`StateSpaceSet`s, as returned by [`basins_fractions`](@ref) or one slice of a
[`global_continuation`](@ref) output) by merging similar ones into groups.

This is the single-parameter building block of [`aggregate_continuation`](@ref); use it directly
when you only want to group the attractors found at one parameter value and there is nothing to
continue.

Each attractor is turned into a feature vector with `featurizer` and the attractors are
partitioned into groups with `group_config` (a [`GroupingConfig`](@ref)). The members of a group
are merged into a single `StateSpaceSet`.

## Return

1. `agg_attractors`: like `attractors`, but mapping each group ID to the merged `StateSpaceSet`
   of its members.
2. `centroids`: a dictionary mapping each group ID to the centroid (mean) of its members' feature
   vectors. This summarises where each group sits in feature space; the per-step version of this
   is the second output of [`aggregate_continuation`](@ref).
3. `members`: a dictionary mapping each group ID to the vector of original attractor IDs (the keys
   of `attractors`) that were merged into it.

To compute stability measures for the merged groups, pass `agg_attractors` to
[`stability_measures_along_continuation`](@ref) (wrapping it and the parameter in length-1
vectors if you have a single parameter).
"""
function aggregate_attractors(attractors::Dict, featurizer, group_config)
    agg_attractors, feature_sets, members = _aggregate_one_step(attractors, featurizer, group_config)
    centroids = _feature_centroids(feature_sets)
    return agg_attractors, centroids, members
end

# Group one parameter step's attractors. Returns `(agg_attractors, feature_sets, members)`:
# the members of a group are merged into one `StateSpaceSet` (`agg_attractors`), the group's
# member feature vectors are collected into a `StateSpaceSet` (`feature_sets`, used to match
# groups across steps in feature space), and `members` maps each group to the original IDs.
function _aggregate_one_step(attractors, featurizer, group_config)
    K = keytype(attractors)
    ids = filter(!isequal(-1), collect(keys(attractors)))
    if isempty(ids)
        return empty(attractors), Dict{K, Any}(), Dict{K, Vector{K}}()
    end
    features = [featurizer(attractors[id]) for id in ids]
    labels = group_features(features, group_config)
    agg_attractors = Dict{K, valtype(attractors)}()
    feature_vecs = Dict{K, Vector{eltype(features)}}()
    members = Dict{K, Vector{K}}()
    for (id, label, f) in zip(ids, labels, features)
        if haskey(agg_attractors, label)
            agg_attractors[label] = StateSpaceSet(
                vcat(collect(agg_attractors[label]), collect(attractors[id]))
            )
        else
            agg_attractors[label] = attractors[id]
        end
        push!(get!(feature_vecs, label, eltype(features)[]), f)
        push!(get!(members, label, K[]), id)
    end
    feature_sets = Dict(label => StateSpaceSet(vecs) for (label, vecs) in feature_vecs)
    return agg_attractors, feature_sets, members
end

# Centroid (mean feature vector) of each group's feature-vector set.
_feature_centroids(feature_sets) = Dict(g => sum(s) / length(s) for (g, s) in feature_sets)

"""
    aggregate_fractions(fractions_cont, members_cont)

Simple convenience function that aggregates the basin fractions resulting from a
[`global_continuation`](@ref), given the members of each aggregated group.
The latter is the output of [`aggregate_continuation`](@ref).
"""
function aggregate_fractions(fractions_cont, members_cont)
    agg_fractions = map(eachindex(members_cont)) do i
        fs = fractions_cont[i]
        members = members_cont[i]
        Dict(k => sum(fs[mk] for mk in m) for (k, m) in members)
    end
    return agg_fractions
end
