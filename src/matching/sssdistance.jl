"""
    MatchBySSDistance <: SSSetMatcher

A matcher type for [`match_continuation!`](@ref) that matches state space sets based on
their distance. This matcher also implements the [`match_statespacesets!`](@ref)
function.

## Keyword arguments

- `distance = Centroid()`: distance to match by, given to [`setsofsets_distances`](@ref).
- `threshold = Inf`: sets with distance larger than the `threshold` are guaranteed
  to not be mapped to each other.


## Description

Sets are matched on a p

"""
@kwdef struct MatchBySSDistance{D, T<:Real}
    distance::D
    threshold::T = Inf
    use_vanished::Bool = !isinf(threshold)
    retract_keys::Bool = true
end

