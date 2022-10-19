"""
A module for finding attractors of dynamical systems,
their basins, and continuing them across parameters.
Part of [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystems.jl/dev/).
"""
module Attractors
using DelayEmbeddings # for datasets, will be replaced by StateSpaceSets.jl
using DynamicalSystemsBase

include("dict_utils.jl")
include("mapping/attractor_mapping.jl")
include("basins/basins.jl")
include("continuation/basins_fractions_continuation_api.jl")
include("sampler.jl")
include("deprecated.jl")

end # module Attractors
