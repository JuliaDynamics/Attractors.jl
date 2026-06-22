export MatchByFeatureDistance

using LinearAlgebra: norm

"""
    MatchByFeatureDistance(; threshold = Inf)

A matcher that matches IDs by the Euclidean distance between feature vectors.

Like any [`IDMatcher`](@ref), it is used as `match_sequentially!(dicts, matcher)`, where
`dicts` is a vector holding one dictionary per parameter step, each mapping an integer ID to
a value. What the value is depends on the matcher: for [`MatchBySSSetDistance`](@ref) the
values are attractors (`StateSpaceSet`s), whereas for `MatchByFeatureDistance` **the values
must be feature vectors** `v` (e.g. `SVector`s). The distance between two IDs is then simply
`norm(v₊ - v₋)`.

In practice the vectors you pass are the per-group feature centroids produced when
aggregating attractors: each group is summarised by the mean of `featurizer(A)` over its
member attractors, and `match_sequentially!` matches those centroids from one parameter step
to the next. This is how [`aggregate_continuation`](@ref) keeps group IDs consistent along the
parameter axis.
The centroids are computed by the aggregation code and handed in as the dictionary values, so
this matcher never featurizes a merged attractor set itself.

Matching uses the same greedy minimum-distance assignment as [`MatchBySSSetDistance`](@ref):
all pairs of (old, new) IDs are sorted by their feature-space distance, the closest pair is
matched first, and matched IDs are removed from the pool to ensure a unique one-to-one
assignment.

## Keyword arguments

- `threshold = Inf`: if the feature-space distance between the best-matching pair exceeds
  `threshold`, the ID in `a₊` is guaranteed to receive a new unique ID rather than being
  matched to an existing one.
"""
@kwdef struct MatchByFeatureDistance{T <: Real} <: IDMatcher
    threshold::T = Inf
end

_use_vanished(m::MatchByFeatureDistance) = !isinf(m.threshold)

function matching_map(
        a₊::AbstractDict, a₋::AbstractDict, matcher::MatchByFeatureDistance;
        kw...
    )
    (isempty(a₊) || isempty(a₋)) && return Dict{keytype(a₊), keytype(a₋)}()
    distances = setsofsets_distances(a₊, a₋, (v₊, v₋) -> norm(v₊ - v₋))
    keys₊, keys₋ = sort.(collect.(keys.((a₊, a₋))))
    return _matching_map_distances(keys₊, keys₋, distances, matcher.threshold; kw...)
end
