# TODO: I don't know how this fits into the new interface...

###########################################################################################
# Matching with basins and possibly overlaps
###########################################################################################
"""
    match_basins_ids!(b₊::AbstractArray, b₋; threshold = Inf)

Similar to [`match_statespacesets!`](@ref) but operate on basin arrays instead
(the arrays typically returned by [`basins_of_attraction`](@ref)).

This method matches IDs of attractors whose basins of attraction before and after `b₋,b₊`
have the most overlap (in pixels). This overlap is normalized in 0-1 (with 1 meaning
100% overlap of pixels). The `threshold` in this case is compared to the inverse
of the overlap (so, for `threshold = 2` attractors that have less than 50% overlap get
different IDs guaranteed).
"""
function match_basins_ids!(b₊::AbstractArray, b₋; threshold = Inf)
    ids₊, ids₋ = unique(b₊), unique(b₋)
    distances = _similarity_from_overlaps(b₊, ids₊, b₋, ids₋)
    rmap = replacement_map(ids₊, ids₋, distances, threshold)
    replace!(b₊, rmap...)
    return rmap
end

function _similarity_from_overlaps(b₊, ids₊, b₋, ids₋)
    size(b₋) == size(b₊) || error("Sizes of b₊ and  b₋ must match")
    distances = Dict{eltype(ids₊), Dict{eltype(ids₋), Float64}}()
    for i in ids₊
        Bi = findall(isequal(i), b₊)
        d = valtype(distances)()
        # Compute normalized overlaps of each basin with each other basis
        for j in ids₋
            Bj = findall(isequal(j), b₋)
            overlap = length(Bi ∩ Bj)/length(Bj)
            d[j] = 1 / overlap # distance is inverse overlap
        end
        distances[i] = d
    end
    return distances
end
