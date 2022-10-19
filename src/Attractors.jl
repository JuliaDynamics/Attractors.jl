"""
A module for finding attractors of dynamical systems,
their basins, and continuing them across parameters.
Part of [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystems.jl/dev/).
"""
module Attractors


include("dict_utils.jl")
include("mapping/attractor_mapping.jl")
include("basins_utilities.jl")
include("fractality_of_basins.jl")
include("tipping.jl")
include("sampler.jl")
include("continuation/basins_fractions_continuation_api.jl")
include("deprecated.jl")

end # module Attractors
