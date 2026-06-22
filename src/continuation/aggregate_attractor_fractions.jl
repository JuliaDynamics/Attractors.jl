export aggregate_attractor_fractions, aggregate_continuation

"""
    aggregate_continuation(
        fractions_cont, attractors_cont, featurizer, group_config [, info_extraction]
    )

Aggregate the output of a [`global_continuation`](@ref) by merging similar attractors into
groups across the parameter range, returning the grouped basin fractions, the merged
attractors, and a per-group descriptor.

This is the post-processing companion of [`global_continuation`](@ref): run the continuation
as usual to obtain `fractions_cont` and `attractors_cont`, then call this function to "group"
similar attractors (for example all states in which a given species is extinct) into a single
aggregated attractor whose basin fraction is the sum of its members'. It is the single entry
point for aggregation; it serves both basin-fraction continuations and, via its merged
attractors, stability-measure continuations (see below).

## How grouping works

At each parameter step the attractors are turned into feature vectors with `featurizer` and
partitioned into groups with `group_config` (a [`GroupingConfig`](@ref)); the members of a
group are merged into a single `StateSpaceSet` and their basin fractions are summed. Each
group is summarised by the centroid of its members' feature vectors, and these centroids are
matched between consecutive parameter steps with [`MatchByFeatureDistance`](@ref), so that a
group keeps the same ID along the parameter axis (it is tracked by the continuity of its
feature centroid). Because matching is centroid-based, `featurizer` must return *numeric*
feature vectors (e.g. `SVector`s).

## Input

1. `fractions_cont`: a vector of dictionaries mapping IDs to basin fractions.
2. `attractors_cont`: a vector of dictionaries mapping IDs to attractors (`StateSpaceSet`s).
   The first two arguments are exactly the return value of [`global_continuation`](@ref)
   with [`AttractorSeedContinueMatch`](@ref) (or of [`basins_fractions`](@ref)).
3. `featurizer`: a 1-argument function mapping an attractor to a numeric feature vector.
   Features expected by [`GroupingConfig`](@ref) are typically `SVector`s.
4. `group_config`: a subtype of [`GroupingConfig`](@ref).
5. `info_extraction`: a function mapping the vector of a group's per-step centroids to a
   descriptor of the group. Optional; defaults to the mean centroid.

## Return

1. `agg_fractions_cont`: like `fractions_cont`, but holding the fractions of the aggregated
   groups.
2. `agg_attractors_cont`: like `attractors_cont`, but mapping each group ID to the merged
   `StateSpaceSet` of its members.
3. `aggregated_info`: a dictionary mapping each group ID to its `info_extraction` descriptor.

All three share consistent group IDs along the parameter axis.

## Aggregating stability measures

To obtain stability measures for the aggregated groups, pass `agg_attractors_cont` to
[`stability_measures_along_continuation`](@ref). Each merged group is then treated as a
single attractor, so every measure — including those that need the raw basin data, such as
medians and critical shock magnitudes — is computed correctly for the group.
"""
function aggregate_continuation(
        fractions_cont::Vector, attractors_cont::Vector, featurizer, group_config,
        info_extraction = mean_across_features, # function from grouping continuation
    )
    P = length(fractions_cont)
    agg_fractions_cont = Vector{Dict{Int, Float64}}(undef, P)
    agg_attractors_cont = Dict[]
    centroids_cont = Dict[]
    for i in 1:P
        agg_fractions_cont[i], agg_attractors, centroids = _aggregate_step(
            fractions_cont[i], attractors_cont[i], featurizer, group_config
        )
        push!(agg_attractors_cont, agg_attractors)
        push!(centroids_cont, centroids)
    end
    # match groups across steps by the Euclidean distance of their feature centroids
    if P > 1
        rmaps = match_sequentially!(centroids_cont, MatchByFeatureDistance())
        match_sequentially!(agg_fractions_cont, rmaps)
        match_sequentially!(agg_attractors_cont, rmaps)
    end
    remove_minus_1_if_possible!(agg_fractions_cont)
    aggregated_info = _aggregated_info(centroids_cont, info_extraction)
    return agg_fractions_cont, agg_attractors_cont, aggregated_info
end

"""
    aggregate_attractor_fractions(
        fractions_cont, attractors_cont, featurizer, group_config [, info_extraction]
    )

Like [`aggregate_continuation`](@ref), but returns only the aggregated basin fractions and the
per-group info, without the merged attractors:

    agg_fractions_cont, aggregated_info = aggregate_attractor_fractions(...)

See [`aggregate_continuation`](@ref) for the description of the arguments and of how attractors
are grouped at each parameter step and matched across the parameter axis.
"""
function aggregate_attractor_fractions(
        fractions_cont::Vector, attractors_cont::Vector, featurizer, group_config,
        info_extraction = mean_across_features,
    )
    agg_fractions_cont, _, aggregated_info = aggregate_continuation(
        fractions_cont, attractors_cont, featurizer, group_config, info_extraction
    )
    return agg_fractions_cont, aggregated_info
end
# convenience wrapper for only single input
function aggregate_attractor_fractions(fractions::Dict, attractors::Dict, args...)
    agg_fractions_cont, aggregated_info = aggregate_attractor_fractions(
        [fractions], [attractors], args...
    )
    return agg_fractions_cont[1], aggregated_info
end

# Group one parameter step's attractors. Returns `(agg_fractions, agg_attractors, centroids)`:
# the members of a group are merged into one `StateSpaceSet`, their fractions are summed, and
# `centroids` maps each group to the mean of its members' feature vectors (used to match groups
# across parameter steps).
function _aggregate_step(fractions, attractors, featurizer, group_config)
    ids = filter(!isequal(-1), collect(keys(attractors)))
    if isempty(ids)
        return copy(fractions), empty(attractors), Dict{Int, Any}()
    end
    features = [featurizer(attractors[id]) for id in ids]
    labels = group_features(features, group_config)
    agg_fractions = Dict{Int, Float64}()
    agg_attractors = Dict{keytype(attractors), valtype(attractors)}()
    feature_sums = Dict{Int, eltype(features)}()
    counts = Dict{Int, Int}()
    for (id, label, f) in zip(ids, labels, features)
        agg_fractions[label] = get(agg_fractions, label, 0.0) + get(fractions, id, 0.0)
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
    if haskey(fractions, -1)
        agg_fractions[-1] = fractions[-1]
    end
    centroids = Dict(label => feature_sums[label] / counts[label] for label in keys(feature_sums))
    return agg_fractions, agg_attractors, centroids
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

function remove_minus_1_if_possible!(afs)
    isthere = false
    for fs in afs
        isthere = get(fs, -1, 0) > 0
        isthere && return
    end
    # no `-1`, we remove from everywhere
    for fs in afs
        delete!(fs, -1)
    end
    return
end
