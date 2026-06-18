export MatchByFeatureDistance

using LinearAlgebra: norm

"""
    MatchByFeatureDistance(featurizer; threshold = Inf)

A matcher that matches attractors by their Euclidean distance in feature space.
`featurizer` is a 1-argument function mapping a `StateSpaceSet` to a feature vector
(the same featurizer used when aggregating attractors, e.g. in [`finalize_accumulator`](@ref)).

The distance between two attractors `A` and `B` is `norm(featurizer(A) - featurizer(B))`.
Matching then follows the same greedy minimum-distance assignment as
[`MatchBySSSetDistance`](@ref): all pairs of (old, new) attractors are sorted by their
feature-space distance, the closest pair is matched first, and matched IDs are removed
from the pool to ensure a unique one-to-one assignment.

This is the natural matcher for grouped/aggregated attractors produced by
[`finalize_accumulator`](@ref) with `featurizer` and `group_config` keywords, since those
groups are precisely characterised by their feature vectors.
[`stability_measures_along_continuation`](@ref) uses this matcher automatically
when `featurizer` and `group_config` are provided.

## Keyword arguments

- `threshold = Inf`: if the feature-space distance between the best-matching pair exceeds
  `threshold`, the attractor in `a₊` is guaranteed to receive a new unique ID rather than
  being matched to an existing one.
"""
struct MatchByFeatureDistance{F, T <: Real} <: IDMatcher
    featurizer::F
    threshold::T
end

MatchByFeatureDistance(featurizer; threshold = Inf) = MatchByFeatureDistance(featurizer, threshold)

_use_vanished(m::MatchByFeatureDistance) = !isinf(m.threshold)

function matching_map(
        a₊::AbstractDict, a₋::AbstractDict, matcher::MatchByFeatureDistance;
        kw...
    )
    isempty(a₊) || isempty(a₋) && return Dict{keytype(a₊), keytype(a₋)}()
    f = matcher.featurizer
    feature_dist = (A, B) -> norm(f(A) - f(B))
    distances = setsofsets_distances(a₊, a₋, feature_dist)
    keys₊, keys₋ = sort.(collect.(keys.((a₊, a₋))))
    return _matching_map_distances(keys₊, keys₋, distances, matcher.threshold; kw...)
end
