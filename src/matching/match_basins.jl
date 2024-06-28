export MatchByBasinOverlap

"""
    MatchByBasinOverlap(threshold = Inf)

A special matcher that matches IDs given full basins of attraction.
It matches IDs of attractors whose basins of attraction before and after `b₋, b₊`
have the most overlap (in pixels). This overlap is normalized in 0-1 (with 1 meaning
100% overlap of pixels). The `threshold` in this case is compared to the inverse
of the overlap, and basins whose overlap is less than `1/threshold` are guaranteed
to get assined different IDs.
For example: for `threshold = 2` basins that have less than 50% overlap get
different IDs guaranteed.

The input for this matcher in [`replacement_map`](@ref)
should be dictionaries mapping IDs to vectors of cartesian indices,
where the indices mean which parts of the state space belong to which ID

This information can be derived from the basins of attraction (`Array`),
i.e., the direct output of [`basins_of_attraction`](@ref).
For convenience, as well as backwards compatibility, when using
[`replacement_map`](@ref) with this mapper you may provide two `Array`s `b₊, b₋`
representing basins of attraction after and before, and the conversion to dictionaries
will happen internally as it is supposed to.
To replace the `IDs` in `b₊` given the replacement map just call `replace!(b₊, rmap...)`,
or use the in-place version [`replacement_map!`](@ref) directly.
"""
struct MatchByBasinOverlap
    threshold::Float64 = Inf
end

function replacement_map(b₊::AbstractArray, b₋::AbstractArray, matcher::MatchByBasinOverlap; i = nothing)
    a₊, a₋ = _basin_to_dict.((b₊, b₋))
    replacement_map(a₊, a₋, matcher; i)
end

function replacement_map!(b₊::AbstractArray, b₋::AbstractArray, matcher::MatchByBasinOverlap; i = nothing)
    rmap = replacement_map(b₊, b₋, matcher; i)
    replace!(b₊, rmap...)
    return rmap
end

# actual implementation
function replacement_map(a₊::AbstractDict, a₋, matcher::MatchByBasinOverlap; i = nothing, next_id = nothing)
    # input checks
    if !(keytype(a₊) <: Vector{<:CartesianIndex})
        throw(ArgumentError("Incorrect input given. For matcher `MatchByBasinOverlap`,
        the dictionaries values should be vectors of `CartesianIndex`."))
    end
    if sum(length, values(a₊)) ≠ sum(length, values(a₋))
        throw(ArgumentError("The sizes of the two basins to be matched must be the same."))
    end
    # The source code of this matcher is beautiful. It computes a "dissimilarity"
    # metric, which is the inverse of the basin overlaps. This "dissimilarity" is
    # just a "distance" between basins of attraction. Thus, it actually
    # propagates this "distance" to the matching code of `MatchBySSDistance`!
    keys₊, keys₋ = keys.((a₊, a₋))
    nextid = isnothing(next_id) ? max(maximum(keys₊), maximum(keys₋)) + 1 : next_id

    distances = Dict{eltype(keys₊), Dict{eltype(keys₋), Float64}}()
    for i in keys₊
        Bi = a₊[i]
        d = valtype(distances)() # d is a dictionary of distances
        # Compute normalized overlaps of each basin with each other basis
        for j in ids₋
            Bj = b₋[j]
            overlap = length(Bi ∩ Bj)/length(Bj)
            d[j] = 1 / overlap # distance is inverse overlap
        end
        distances[i] = d
    end
    _replacement_map_distances(keys₊, keys₋, distances, matcher.threshold, nextid)
end

function _basin_to_dict(b::AbstractArray{Int})
    ukeys = unique(b)
    d = Dict(k => findall(isequal(k), b) for k in ukeys)
    return d
end

# TODO: test that it works also with vector of basins of attraction
