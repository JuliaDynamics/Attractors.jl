export continuation

"""
   AttractorsBasinsContinuation

Supertype of continuation methods given to [`continuation`](@ref).
Possible subtypes are:

- [`RecurrencesContinuation`](@ref)
- [`FeaturizingContinuation`](@ref)

These are given to the main [`continuation`](@ref) function.
These types can also obtain a [`MatchingMethod`](@ref) that instructs how
to match attractor IDs across the different parameters.
"""
abstract type AttractorsBasinsContinuation end

"""
   continuation(abc::AttractorsBasinsContinuation, prange, pidx, ics; kwargs...)

Find and continue attractors and the fractions of their basins of attraction
across a parameter range.

The continuation type `abc` is a subtype of [`AttractorsBasinsContinuation`](@ref)
and contains an [`AttractorMapper`](@ref). The mapper contains information
on how to find the attractors and basins of a dynamical system. Additional
arguments and keyword arguments given when creating `abc` further tune the continuation.

In the `continuation` function, the basin fractions and the attractors (or some representation
of them in the case of featurizing) are continued across the parameter range `prange`,
for the parameter of the system with index `pidx`.
`ics` is as in [`basins_fractions`](@ref), i.e., it is either a function generating
initial conditions or a set containing them.

## Return

1. `fractions_curves :: Vector{Dict{Int, Float64}}`. The fractions of basins of attraction.
   `fractions_curves[i]` is a dictionary mapping attractor IDs to their basin fraction
   at the `i`-th parameter.
2. `attractors_info <: Vector{Dict{Int, <:Any}}`. Information about the attractors.
   `attractors_info[i]` is a dictionary mapping attractor ID to information about the
   attractor at the `i`-th parameter.
   The type of information stored depends on the chosen continuation type.

## Keyword arguments

- `show_progress = true`: display a progress bar of the computation.
* `samples_per_parameter = 100`: amount of initial conditions sampled at each parameter.
"""
function continuation end

"""
   MatchingMethod

Supertype of types specifying how to match attractors/features in [`continuation`](@ref).
Concrete subtypes are given as options when instantiating an
[`AttractorsBasinsContinuation`](@ref) subtype.

Below we outline the possible matching options. To understand them, you need to abstract
an attractor and a vector of feature vectors into the same thing: a [`StateSpaceSet`](@ref).
Attractors are what is found by [`AttractorsViaRecurrences`](@ref).
Vectors of features (which are also essentially [`StateSpaceSet`](@ref)),
is what is found by [`AttractorsViaFeaturizing`](@ref) because all features that
get assigned the same label form a vector of features.

Let a/f/set mean either an attractor or a vector of feature vectors with same label.
So, a/f/set means a [`StateSpaceSet'](@ref).

The matching of a/f/sets across parameter values is dictactated by the following
subtypes of `MatchingMethod`:

- `ParameterSliceCrossDistance()`: At each parameter slice beyond the first, the new
  a/f/sets are matched to the previous a/f/sets found in the previous parameter value
  by a direct call to the [`match_attractor_ids!`](@ref) function. Hence, the matching
  of a/f/sets here works "slice by slice" on the parameter axis and the a/f/sets
  that are closest to each other (in state space, but for two different parameter values)
  get assigned the same label.

- `ClusterDistanceMatrix()`: There is a slight difference of how this works
  The a/f/sets are grouped over the full parameter range
  using a DBSCAN clustering. A distance matrix is created over all a/f/sets across
  parameter values, using the [`set_distance`](@ref) function. This distance matrix
  is given to DBSCAN, and the output is clusterred attractors. Now each cluster may
  include attractors across different parameter values. After the clustering is finished
  the cluster label fractions are distributed to each parameter value they came from.

TODO: This is WRONG!!! The `ClusterDistanceMatrix` is something completely
different in the Recurrences and the Featurizing versions...

"""
abstract type MatchingMethod end

struct ParameterSliceCrossDistance end

struct ClusterDistanceMatrix end


include("match_attractor_ids.jl")
include("continuation_recurrences.jl")
include("continuation_featurizing.jl")
include("aggregate_attractor_fractions.jl")
