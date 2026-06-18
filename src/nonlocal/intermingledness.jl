export intermingledness

"""
    intermingledness(points::StatesSpaceSet, labels [, distance]; kw...)

Calculate the intermingledness [Datseris2026](@cite) of the `points`
which have been divided into groups (typically attractors) indicated by `labels`.
Return a dictionary mapping unique labels to their intermingledness.

The optional `distance = Euclidean()` argument is a function dictating
how to estimate distances between points.
A vector of distance functions can also be given as `distance`, in which case a vector of
intermingledness is returned corresponding to each distance function.
In [Datseris2026](@cite) intermingledness was estimated individually per dimension
of `points`. You can achieve this by simply passing in `points[:, x]` with `x`
the dimension's index you care about.

The `summarizer = mean` keyword argument dictates how
to summarize the intermingedness statistic across other groups.

## Description

Intermingledness is a way to quantify similarity or dissimilarity between the different
groups that `points` was grouped in. For example, `points`
can be initial conditions fed into [`basin_fractions`](@ref), and `labels` the output.
Or, `points` can be feature vectors and `labels` the output of [`group_features`](@ref).

Intermingledness is effectively the ratio of the pairwise-averaged intra-group distance
divided by the pairwise-averaged inter-group distance. If it is 1, points are as
close to points in their own group as they are to points in other groups.
See [Datseris2026](@cite) for examples using intermingledness and the detailed definition
or honestly, just look at the source code, it is only 10 lines!

!!! note "Expensive!"
    This function becomes expensive to compute for many points
    because it scales as ~ `length(unique(labels))^2 * length(points)^2`
"""
function intermingledness(
        us::AbstractArray{<:AbstractArray}, labels::AbstractArray{<:Integer},
        distance = Euclidean(); summarizer = mean
    )
    return only(intermingledness(us, labels, [distance]; summarizer))
end
function intermingledness(
        us::AbstractArray{<:AbstractArray}, labels::AbstractArray{<:Integer},
        distances::AbstractVector; summarizer = mean
    )
    length(us) ≠ length(labels) && error("points and labels must be same length.")
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

# Pairwise average distance between all points in xs and ys
mean_distance(xs, ys, distance) = mean(distance(x, y) for x in xs for y in ys)

function boundary_intermingledness(
        us::AbstractArray{<:AbstractArray}, labels::AbstractArray{<:Integer},
        distance = Euclidean(); summarizer = mean
    )

    # Unique group labels and the point indices belonging to each label.
    ulabels = unique(labels)
    L = length(ulabels)
    group_indices = Dict(l => findall(isequal(l), labels) for l in ulabels)
    bsets = boundary_sets(us, ulabels, group_indices, distance)
    group_sizes = Dict(l => length(v) for (l, v) in group_indices)

    # diagonals by definition zero
    imatrix = zeros(L, L)
    for i in 1:L-1
        la = ulabels[i]
        for j in i+1:L
            lb = ulabels[j]
            v = length(bsets[(la,lb)])/(group_sizes[la] + group_sizes[lb])
            imatrix[i, j] = v
            imatrix[j, i] = v
        end
    end
    # summarize per label
    return Dict(ulabels[i] => summarizer(imatrix[i, setdiff(1:L, i)]) for i in 1:L)
end

import NearestNeighbors
function boundary_sets(points, ulabels, group_indices, distance)

    L = length(ulabels)
    label_type = eltype(ulabels)

    # Precompute per-label indices, point subsets and KD-trees once.
    idxs = [group_indices[l] for l in ulabels]
    groups = [points[idx] for idx in idxs]
    trees = [NearestNeighbors.KDTree(X, distance) for X in groups]

    # For each label pair (la, lb), store a single (symmetric) boundary index set.
    boundary_sets = Dict{Tuple{label_type, label_type}, Vector{Int}}()
    sizehint!(boundary_sets, L * (L - 1) ÷ 2)

    # Iterate over all unique unordered label pairs.
    for i in 1:L-1
        la = ulabels[i]
        idxA = idxs[i]
        XA = groups[i]
        treeA = trees[i]

        for j in i+1:L
            lb = ulabels[j]
            idxB = idxs[j]
            XB = groups[j]
            treeB = trees[j]

            # A -> B -> A:
            # 1) for each A-point, find its nearest neighbor in B
            # 2) keep unique B neighbors
            # 3) map those B neighbors back to their nearest in A
            # Under symmetry, this defines the single boundary set for pair (la, lb).
            b_nns, _ = NearestNeighbors.nn(treeB, XA)
            b_unique = unique!(b_nns)
            a_back, _ = NearestNeighbors.nn(treeA, @view XB[b_unique])
            boundaryA = unique!(idxA[a_back])

            # Include both sides participating in the same boundary relation.
            boundary_sets[(la, lb)] = unique!(vcat(boundaryA, idxB[b_unique]))
        end
    end

    return boundary_sets
end
