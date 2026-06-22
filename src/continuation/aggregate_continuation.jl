export aggregate_continuation

"""
    aggregate_continuation(
        attractors_cont, featurizer, group_config [, info_extraction]
    ) → agg_attractors_cont, aggregated_info

Aggregate the attractors of a [`global_continuation`](@ref) by merging similar ones into groups
across the parameter range.

This is the post-processing companion of [`global_continuation`](@ref): run the continuation as
usual to obtain `attractors_cont`, then call this function to "group" similar attractors (for
example all states in which a given species is extinct) into a single merged attractor. The
merged attractors carry IDs that stay consistent along the parameter axis, so they can be fed
directly to [`stability_measures_along_continuation`](@ref) to obtain stability measures — the
basin fractions of the groups included — for the aggregated groups.

## How grouping works

At each parameter step the attractors are turned into feature vectors with `featurizer` and
partitioned into groups with `group_config` (a [`GroupingConfig`](@ref)); the members of a group
are merged into a single `StateSpaceSet`. Each group is summarised by the centroid of its
members' feature vectors, and these centroids are matched between consecutive parameter steps
with [`MatchByFeatureDistance`](@ref), so that a group keeps the same ID along the parameter axis
(it is tracked by the continuity of its feature centroid). Because matching is centroid-based,
`featurizer` must return *numeric* feature vectors (e.g. `SVector`s).

## Input

1. `attractors_cont`: a vector of dictionaries mapping IDs to attractors (`StateSpaceSet`s),
   exactly as returned by [`global_continuation`](@ref) (or [`basins_fractions`](@ref)).
2. `featurizer`: a 1-argument function mapping an attractor to a numeric feature vector.
   Features expected by [`GroupingConfig`](@ref) are typically `SVector`s.
3. `group_config`: a subtype of [`GroupingConfig`](@ref).
4. `info_extraction`: a function mapping the vector of a group's per-step centroids to a
   descriptor of the group. Optional; defaults to the mean centroid.

## Return

1. `agg_attractors_cont`: like `attractors_cont`, but mapping each group ID to the merged
   `StateSpaceSet` of its members.
2. `aggregated_info`: a dictionary mapping each group ID to its `info_extraction` descriptor.

Both share consistent group IDs along the parameter axis.

## Aggregating stability measures

To obtain stability measures (including the basin fractions) for the aggregated groups, pass
`agg_attractors_cont` to [`stability_measures_along_continuation`](@ref). Each merged group is
then treated as a single attractor, so every measure — including those that need the raw basin
data, such as medians and critical shock magnitudes — is computed correctly for the group.
"""
function aggregate_continuation(
        attractors_cont::Vector, featurizer, group_config,
        info_extraction = mean_across_features, # function from grouping continuation
    )
    P = length(attractors_cont)
    agg_attractors_cont = Dict[]
    centroids_cont = Dict[]
    for i in 1:P
        agg_attractors, centroids = _aggregate_attractors_step(
            attractors_cont[i], featurizer, group_config
        )
        push!(agg_attractors_cont, agg_attractors)
        push!(centroids_cont, centroids)
    end
    # match groups across steps by the Euclidean distance of their feature centroids
    if P > 1
        rmaps = match_sequentially!(centroids_cont, MatchByFeatureDistance())
        match_sequentially!(agg_attractors_cont, rmaps)
    end
    aggregated_info = _aggregated_info(centroids_cont, info_extraction)
    return agg_attractors_cont, aggregated_info
end

# Group one parameter step's attractors. Returns `(agg_attractors, centroids)`: the members of a
# group are merged into one `StateSpaceSet`, and `centroids` maps each group to the mean of its
# members' feature vectors (used to match groups across parameter steps).
function _aggregate_attractors_step(attractors, featurizer, group_config)
    ids = filter(!isequal(-1), collect(keys(attractors)))
    if isempty(ids)
        return empty(attractors), Dict{Int, Any}()
    end
    features = [featurizer(attractors[id]) for id in ids]
    labels = group_features(features, group_config)
    agg_attractors = Dict{keytype(attractors), valtype(attractors)}()
    feature_sums = Dict{Int, eltype(features)}()
    counts = Dict{Int, Int}()
    for (id, label, f) in zip(ids, labels, features)
        if haskey(agg_attractors, label)
            agg_attractors[label] = StateSpaceSet(
                vcat(collect(agg_attractors[label]), collect(attractors[id]))
            )
        else
            agg_attractors[label] = attractors[id]
        end
        feature_sums[label] = get(feature_sums, label, zero(f)) + f
        counts[label] = get(counts, label, 0) + 1
    end
    centroids = Dict(label => feature_sums[label] / counts[label] for label in keys(feature_sums))
    return agg_attractors, centroids
end

# Summarise each group across the continuation by applying `info_extraction` to the vector of
# the group's per-step centroids.
function _aggregated_info(centroids_cont, info_extraction)
    ids = setdiff(unique_keys(centroids_cont), [-1])
    return Dict(
        id => info_extraction([
            centroids_cont[i][id]
            for i in eachindex(centroids_cont) if haskey(centroids_cont[i], id)
        ])
        for id in ids
    )
end
