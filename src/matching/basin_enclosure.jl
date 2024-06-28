export MatchByBasinEnclosure

"""
    MatchByBasinEnclosure

TODO: Kalel.
"""
struct MatchByBasinEnclosure
    mapper
    prange
    pidx
    T
    coflowing_threshold
end

function replacement_map(a₊::AbstractDict, a₋::AbstractDict, matcher::MatchByBasinEnclosure;
        i = nothing, # keyword `i` is not used by this mapper
        next_id = next_free_id(a₊, a₋) # but next_id is propagated
    )
    if isnothing(i)
        throw(ArgumentError("When using `MatchByBasinEnclosure` you need to provide the parameter index
        of the `a₋` dictionary"))
    end
    # set ds to the "new" parameter
    set_parameter!(matcher.mapper.ds, matcher.pidx, matcher.prange[i])
    # It is not clear to me which point of the attractors to integrate forwards
    # so I just pick the last point
    keys₋ = sort(keys(a₋))
    initconds = [a₋[k][end] for k in keys₋]
    # we then use this initial conditions to map them to attractors

    # integrate the initial conditions forwards in the NEW parameter

    integrated_endpoints = Dict(k => trajectory(matcher.ds, matcher.T, A[end])[1][end] for (k, A) in a₋)

end



distances = setsofsets_distances(a₊, a₋, matcher.distance)
keys₊, keys₋ = sort.(collect.(keys.((a₊, a₋))))
 _replacement_map_distances(keys₊, keys₋, distances::Dict, matcher.threshold; kw...)
end
