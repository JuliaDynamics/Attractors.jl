module Attractors

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end Attractors

using Reexport

@reexport using StateSpaceSets
@reexport using DynamicalSystemsBase

# main files that import other files
include("dict_utils.jl")
include("mapping/attractor_mapping.jl")
include("basins/basins.jl")
include("continuation/basins_fractions_continuation_api.jl")
include("matching/matching_interface.jl")
include("boundaries/edgetracking.jl")
include("deprecated.jl")
include("tipping/tipping.jl")

# Visualization (export names extended in the extension package)
include("plotting.jl")

# minimal fatal shock algo

end # module Attractors
