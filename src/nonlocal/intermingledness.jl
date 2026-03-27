export intermingledness

"""
    intermingledness(us::StatesSpaceSet, labels [, distance]; kw...)

Return the intermingledness [Datseris2026](@cite) of the points in `us`
which have been divided into groups (typically attractors) as dictated by the `labels`.

The optional `distance = Euclidean()` argument dictates how to estimate distances
between points in `us`.
A vector of distances can also be given as `distance`, in which case a vector of
intermingledness is return corresponding to each distance version.

The `summarizer = maximum` keyword argument ditactes how
to summarize the intermingedness statistic across other groups (see description below).

## Description

For example,
`us` can be initial conditions fed into [`basin_fractions`](@ref), and `labels` the output.
Or, `us` can be feature vectors and `labels` the output of the [`group_features`](@ref) function.
"""
function intermingledness(
        us::AbstractStateSpaceSet, labels::AbstractVector{<:Int},
        distance = Euclidean(); summarizer = maximum
    )
    ukeys = unique(labels)
    groups = [us[findall(isequal(gi), labels)] for gi in ukeys]
    return _intermingledness(ukeys, groups, distance, summarizer)
end
function intermingledness(
        us::AbstractStateSpaceSet, labels::AbstractVector{<:Int},
        distances::AbstractVector; summarizer = maximum
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
        return ukeys[gi] => summarizer(imetrics)
    end
    return Dict(imetric) # make sure this is a dictionary so that labels are respected
end

mean_distance(xs, ys, distance) = mean(distance(x, y) for x in xs for y in ys)
