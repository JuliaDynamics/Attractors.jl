"""
A module for finding attractors of dynamical systems,
their basins and their boundaris,
and continuing attractors and basins across parameters.

Part of [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystems.jl/dev/).
"""
module Attractors

using Reexport
@reexport using StateSpaceSets
@reexport using DynamicalSystemsBase

include("dict_utils.jl")
include("mapping/attractor_mapping.jl")
include("basins/basins.jl")
include("continuation/basins_fractions_continuation_api.jl")
include("deprecated.jl")

end # module Attractors
