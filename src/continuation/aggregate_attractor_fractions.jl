export aggregate_attractor_fractions

"""
    aggregate_attractor_fractions(
        fractions_cont, attractors_cont, featurizer, group_config [, info_extraction]
    )

Aggregate the already-estimated curves of fractions of basins of attraction of similar
attractors using the same pipeline used by [`GroupingConfig`](@ref).
The most typical application of this function is to transform the output of a
[`global_continuation`](@ref) with [`AttractorSeedContinueMatch`](@ref) so that similar
attractors, even across parameter space, are grouped into one "attractor".
Thus, the fractions of their basins are aggregated (joined).

You could also use this function to aggregate attractors and their fractions even in
a single parameter configuration, i.e., using the output of [`basins_fractions`](@ref).

This function is useful in cases where you want the accuracy and performance of
[`AttractorsViaRecurrences`](@ref), but you also want the convenience of "grouping"
similar attractors like in [`AttractorsViaFeaturizing`](@ref) for presentation or
analysis purposes. For example, a high dimensional model of competition dynamics
with multiple species may have extreme multistability. After finding this multistability
however, one may care about aggregating all attractors into two groups: where a given
species is extinct or not. This is the example highlighted in our documentation,
in [Extinction of a species in a multistable competition model](@ref aggregation_example).

## Input

1. `fractions_cont`: a vector of dictionaries mapping labels to basin fractions.
2. `attractors_cont`: a vector of dictionaries mapping labels to attractors.
   1st and 2nd argument are exactly like the return values of
   [`global_continuation`](@ref) with [`RecurrencesFindAndMatch`](@ref)
   (or, they can be the return of [`basins_fractions`](@ref)).
3. `featurizer`: a 1-argument function to map an attractor into an appropriate feature
   to be grouped later. Features expected by [`GroupingConfig`](@ref) are `SVector`.
4. `group_config`: a subtype of [`GroupingConfig`](@ref).
5. `info_extraction`: a function accepting a vector of features and returning a description
   of the features. I.e., exactly as in [`FeaturizeGroupAcrossParameter`](@ref).
   The 5th argument is optional and defaults to the centroid of the features.

## Return

1. `aggregated_fractions`: same as `fractions_cont` but now contains the fractions of the
   aggregated attractors.
2. `aggregated_info`: dictionary mapping the new labels of `aggregated_fractions` to the
   extracted information using `info_extraction`.

## Clustering attractors directly

_(this is rather advanced)_

You may also use the DBSCAN clustering approach here to group attractors
based on their state space distance (the [`set_distance`](@ref)) by making a distance
matrix as expected by the DBSCAN implementation.
For this, use `identity` as `featurizer`, and choose [`GroupViaClustering`](@ref)
as the `group_config` with `clust_distance_metric = set_distance` and provide a numerical
value for `optimal_radius_method` when initializing the [`GroupViaClustering`](@ref),
and also, for the `info_extraction` argument, you now need to provide a function that
expects a _vector of `StateSpaceSet`s_ and outputs a descriptor.
E.g., `info_extraction = vector -> mean(mean(x) for x in vector)`.
"""
function aggregate_attractor_fractions(
        fractions_cont::Vector, attractors_cont::Vector, featurizer, group_config,
        info_extraction = mean_across_features # function from grouping continuation
    )
    P = length(fractions_cont)

    # Aggregate at each parameter step independently using featurizer + group_config.
    # centroids_cont[i] maps group ID → centroid of the constituent attractors' features.
    agg_fractions_cont = Vector{Dict{Int, Float64}}(undef, P)
    centroids_cont = Dict[]
    for i in 1:P
        agg_fractions_cont[i], centroids_i = _aggregate_fractions_at_step(
            fractions_cont[i], attractors_cont[i], featurizer, group_config
        )
        push!(centroids_cont, centroids_i)
    end

    # Match group labels across steps by minimising centroid distance in feature space.
    if P > 1
        rmaps = match_sequentially!(centroids_cont, MatchByFeatureDistance(identity))
        match_sequentially!(agg_fractions_cont, rmaps)
    end

    remove_minus_1_if_possible!(agg_fractions_cont)

    # Compute representative descriptors as info_extraction over per-step centroids.
    ids = setdiff(unique_keys(centroids_cont), [-1])
    aggregated_info = Dict(
        id => info_extraction([
            centroids_cont[i][id]
            for i in 1:P if haskey(centroids_cont[i], id)
        ])
        for id in ids
    )

    return agg_fractions_cont, aggregated_info
end

# convenience wrapper for only single input
function aggregate_attractor_fractions(fractions::Dict, attractors::Dict, args...)
    aggregated_fractions, aggregated_info = aggregate_attractor_fractions(
        [fractions], [attractors], args...
    )
    return aggregated_fractions[1], aggregated_info
end

# Aggregate fractions at a single parameter step and compute per-group feature centroids.
# Returns (agg_fractions, centroids) where centroids maps group ID → mean of constituent
# attractors' feature vectors. Does NOT build merged StateSpaceSets.
function _aggregate_fractions_at_step(fractions, attractors, featurizer, group_config)
    ids = filter(!isequal(-1), collect(keys(attractors)))
    if isempty(ids)
        return copy(fractions), Dict{Int, Any}()
    end
    features = [featurizer(attractors[id]) for id in ids]
    grouped_labels = group_features(features, group_config)
    agg_fractions = Dict{Int, Float64}()
    feature_sums = Dict{Int, typeof(features[1])}()
    counts = Dict{Int, Int}()
    for (old_id, new_id, f) in zip(ids, grouped_labels, features)
        agg_fractions[new_id] = get(agg_fractions, new_id, 0.0) + get(fractions, old_id, 0.0)
        feature_sums[new_id] = get(feature_sums, new_id, zero(f)) + f
        counts[new_id] = get(counts, new_id, 0) + 1
    end
    if haskey(fractions, -1)
        agg_fractions[-1] = fractions[-1]
    end
    centroids = Dict(label => feature_sums[label] / counts[label] for label in keys(feature_sums))
    return agg_fractions, centroids
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
