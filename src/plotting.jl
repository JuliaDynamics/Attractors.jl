##########################################################################################
# Basins
##########################################################################################
"""
    heatmap_basins_attractors(grid, basins, attractors; kwargs...)

Plot a heatmap of found (2-dimensional) `basins` of attraction and corresponding
`attractors`, i.e., the output of [`basins_of_attraction`](@ref).

Keywords:

- `projection_into_2D = (A) -> (A[:, 1], A[:, 2])`: how to extract two dimensions
  from the attractors (the plot is 2-dimensional).
- Also all [common plotting keywords](@ref).
"""
function heatmap_basins_attractors end
function heatmap_basins_attractors! end
export heatmap_basins_attractors, heatmap_basins_attractors!

##########################################################################################
# Attractors
##########################################################################################
function animate_attractors_continuation end