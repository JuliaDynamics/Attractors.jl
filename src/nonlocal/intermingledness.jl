export intermingledness

"""
    intermingledness(points::StatesSpaceSet, labels [, distance]; kw...)

Return the intermingledness [Datseris2026](@cite) of the `points`
which have been divided into groups (typically attractors) as dictated by the `labels`.
Return a dictionary mapping unique labels to their intermingledness.

The optional `distance = Euclidean()` argument dictates how to estimate distances
between points.
A vector of distance functions can also be given as `distance`, in which case a vector of
intermingledness is returned corresponding to each distance version.
In [Datseris2026](@cite) intermingledness was estimated individually per dimension
of `points`, which you can achieve by e.g.,
```julia
using Distances: WeightedEuclidean
D = dimension(points)
weights = [(1:D .== i) for i in 1:D]
distances = WeightedEuclidean.(weights)
```

The `summarizer = mean` keyword argument dictates how
to summarize the intermingedness statistic across other groups.

## Description

Intermingledness is a way to quantify similarity or dissimilarity between the different
groups that `points` was grouped in. For example, `points`
can be initial conditions fed into [`basin_fractions`](@ref), and `labels` the output.
Or, `points` can be feature vectors and `labels` the output of [`group_features`](@ref).

Intermingledness is effectively the ratio of the pairwise-averaged inter-group distance
divided by the pairwise-averaged intra-group distance.
See [Datseris2026](@cite) for examples using intermingledness and the detailed definition
or honestly, just look at the source code, it is only 10 lines!

!!! note "Expensive!"
    This function becomes quite expensive to compute for many points
    because it scales as ~ `length(unique(labels))^2 * length(points)^2`
"""
function intermingledness(
        us::AbstractVector{<:AbstractArray}, labels::AbstractVector{<:Int},
        distance = Euclidean(); summarizer = mean
    )
    ukeys = unique(labels)
    groups = [us[findall(isequal(gi), labels)] for gi in ukeys]
    return _intermingledness(ukeys, groups, distance, summarizer)
end
function intermingledness(
        us::AbstractVector{<:AbstractArray}, labels::AbstractVector{<:Int},
        distances::AbstractVector; summarizer = mean
    )
    ukeys = unique(labels)
    groups = [us[findall(isequal(gi), labels)] for gi in ukeys]
    return map(d -> _intermingledness(ukeys, groups, d, summarizer), distances)
end

function _intermingledness(ukeys, groups, distance, summarizer)
    # for each group...
    imetric = map(eachindex(groups)) do gi
        g = groups[gi]
        # mean distance of current group to all groups (including self)
        ds = mean_distance.(Ref(g), groups, Ref(distance))
        # the metric is now the distance of own group divided
        # by distance to any other cluster
        imetrics = ds[gi] ./ ds
        # Right, but now we still need to return a summarizing number across
        # all other groups beyond self group.
        # First we drop the same group entry (which is 1 by definition)
        deleteat!(imetrics, gi)
        # and then summarize
        value = isempty(imetrics) ? NaN : summarizer(imetrics) # NaN if only 1 group
        return ukeys[gi] => value
    end
    return Dict(imetric) # make sure this is a dictionary so that labels are respected
end

mean_distance(xs, ys, distance) = mean(distance(x, y) for x in xs for y in ys)
