export aggregate_attractor_fractions

"""
    aggregate_attractor_fractions(
        fractions_curves, attractors_info, featurizer, group_config
    )

Aggregate the already-estimated curves of fractions of basins of attraction of similar
attractors using the same pipeline used by [`GroupAcrossParameterContinuation`](@ref).
The most typical application of this function is to transform the output of
[`RecurrencesSeedingContinuation`](@ref) so that similar attractors, even across parameter
space, are grouped into one attractor. Thus, the fractions of their basins are aggregated.

For example... (add here Kalels example for ecosystem dynamics).

## Input
The function accepts the following input arguments:

## Return
joint_fractions
"""
function aggregate_attractor_fractions(
        fractions_curves, attractors_info, featurizer, group_config
    )

    original_labels, unlabeled_fractions, parameter_idxs, features =
    refactor_into_sequential_features(fractions_curves, attractors_info, featurizer)

    grouped_labels = group_features(features, group_config)

    joint_fractions = reconstruct_joint_fractions(
        fractions_curves, original_labels, grouped_labels, parameter_idxs, unlabeled_fractions
    )
    # TODO: Remove label -1 if it is empty in all configurations.
    return joint_fractions
end

# TODO: I also need to return information about the grouped features. Same way as in
# the continuation method.

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
    joint_fractions = [Dict{Int,Float64}() for _ in 1:length(fractions_curves)]
    current_p_idx = 0
    for j in eachindex(grouped_labels)
        new_label = grouped_labels[j]
        p_idx = parameter_idxs[j]
        if p_idx > current_p_idx
            current_p_idx += 1
            joint_fractions[current_p_idx][-1] = unlabeled_fractions[current_p_idx]
        end
        d = joint_fractions[current_p_idx]
        orig_frac = get(fractions_curves[current_p_idx], original_labels[j], 0)
        d[new_label] = get(d, new_label, 0) + orig_frac
    end
    return joint_fractions
end