export aggregate_attractor_fractions

"""
    aggregate_attractor_fractions(
        fractions_curves, attractors_info, featurizer, group_config [, info_extraction]
    )

Aggregate the already-estimated curves of fractions of basins of attraction of similar
attractors using the same pipeline used by [`GroupingConfig`](@ref).
The most typical application of this function is to transform the output of
[`RecurrencesFindAndMatch`](@ref) so that similar attractors, even across parameter
space, are grouped into one "attractor". Thus, the fractions of their basins are aggregated.

You could also use this function to aggregate attractors and their fractions even in
a single parameter configuration, i.e., using the output of [`basins_fractions`](@ref).

This function is useful in cases where you want the accuracy and performance of
[`AttractorsViaRecurrences`](@ref), but you also want the convenience of "grouping"
similar attractrors like in [`AttractorsViaFeaturizing`](@ref) for presentation or
analysis purposes. For example, a high dimensional model of competition dynamics
across multispecies may have extreme multistability. After finding this multistability
however, one may care about aggregating all attractors into two groups: where a given
species is extinct or not. This is the example highlighted in our documentation,
in [Extinction of a species in a multistable competition model](@ref).

## Input

1. `fractions_curves`: a vector of dictionaries mapping labels to basin fractions.
2. `attractors_info`: a vector of dictionaries mapping labels to attractors.
   1st and 2nd argument are exactly like the return values of
   [`continuation`](@ref) with [`RecurrencesFindAndMatch`](@ref)
   (or, they can be the return of [`basins_fractions`](@ref)).
3. `featurizer`: a 1-argument function to map an attractor into an appropriate feature
   to be grouped later. Features expected by [`GroupingConfig`](@ref) are `SVector`.
4. `group_config`: a subtype of [`GroupingConfig`](@ref).
5. `info_extraction`: a function accepting a vector of features and returning a description
   of the features. I.e., exactly as in [`FeaturizeGroupAcrossParameter`](@ref).
   The 5th argument is optional and defaults to the centroid of the features.

## Return

1. `aggregated_fractions`: same as `fractions_curves` but now contains the fractions of the
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
        fractions_curves::Vector, attractors_info::Vector, featurizer, group_config,
        info_extraction = mean_across_features # function from grouping continuation
    )

    original_labels, unlabeled_fractions, parameter_idxs, features =
    refactor_into_sequential_features(fractions_curves, attractors_info, featurizer)

    grouped_labels = group_features(features, group_config)

    aggregated_fractions = reconstruct_joint_fractions(
        fractions_curves, original_labels, grouped_labels, parameter_idxs, unlabeled_fractions
    )
    remove_minus_1_if_possible!(aggregated_fractions)
    aggregated_info = info_of_grouped_features(features, grouped_labels, info_extraction)
    return aggregated_fractions, aggregated_info
end
# convenience wrapper for only single input
function aggregate_attractor_fractions(fractions::Dict, attractors::Dict, args...)
    aggregated_fractions, aggregated_info = aggregate_attractor_fractions(
        [fractions], [attractors], args...
    )
    return aggregated_fractions[1], aggregated_info[1]
end


function refactor_into_sequential_features(fractions_curves, attractors_info, featurizer)
    # Set up containers
    P = length(fractions_curves)
    example_feature = featurizer(first(values(attractors_info[1])))
    features = typeof(example_feature)[]
    original_labels = keytype(first(fractions_curves))[]
    parameter_idxs = Int[]
    unlabeled_fractions = zeros(P)
    # Transform original data into sequential vectors
    for i in eachindex(fractions_curves)
        fs = fractions_curves[i]
        ai = attractors_info[i]
        A = length(ai)
        append!(parameter_idxs, (i for _ in 1:A))
        unlabeled_fractions[i] = get(fs, -1, 0.0)
        for k in keys(ai)
            push!(original_labels, k)
            push!(features, featurizer(ai[k]))
        end
    end
    # `parameter_idxs` is the indices of the parameter a given feature maps to.
    # necessary because for some parameter values we may have more or less attractors
    # and hence more or less features
    return original_labels, unlabeled_fractions, parameter_idxs, features
end

function reconstruct_joint_fractions(
        fractions_curves, original_labels, grouped_labels, parameter_idxs, unlabeled_fractions
    )
    aggregated_fractions = [Dict{Int,Float64}() for _ in 1:length(fractions_curves)]
    current_p_idx = 0
    for j in eachindex(grouped_labels)
        new_label = grouped_labels[j]
        p_idx = parameter_idxs[j]
        if p_idx > current_p_idx
            current_p_idx += 1
            aggregated_fractions[current_p_idx][-1] = unlabeled_fractions[current_p_idx]
        end
        d = aggregated_fractions[current_p_idx]
        orig_frac = get(fractions_curves[current_p_idx], original_labels[j], 0)
        d[new_label] = get(d, new_label, 0) + orig_frac
    end
    return aggregated_fractions
end

function info_of_grouped_features(features, grouped_labels, info_extraction)
    ids = sort!(unique(grouped_labels))
    # remove -1 if it's there
    if ids[1] == -1; popfirst!(ids); end
    # extract the infos
    Dict(
        id => info_extraction(
            view(features, findall(isequal(id), grouped_labels))
        ) for id in ids
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
