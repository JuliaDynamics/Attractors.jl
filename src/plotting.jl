##########################################################################################
# Basins
##########################################################################################
"""
    heatmap_basins_attractors(grid, basins, attractors; kwargs...)

Plot a heatmap of found (2-dimensional) `basins` of attraction and corresponding
`attractors`, i.e., the output of [`basins_of_attraction`](@ref).

## Keyword arguments

- All the [common plotting keywords](@ref).
"""
function heatmap_basins_attractors end
function heatmap_basins_attractors! end
export heatmap_basins_attractors, heatmap_basins_attractors!

##########################################################################################
# Continuation
##########################################################################################
"""
    animate_attractors_continuation(
        ds::DynamicalSystem, attractors_info, fractions_curves, prange, pidx;
        kwargs...
    )

Animate how the found system attractors and their corresponding basin fractions
change as the system parameter is increased. This function combines the input
and output of the [`continuation`](@ref) function into a video output.

The input dynamical system `ds` is used to evolve initial conditions sampled from the
found attractors, so that the attractors are better visualized.
`attractors_info, fractions_curves` are the output of [`continuation`](@ref)
while `ds, prange, pidx` are the input to [`continuation`](@ref).

## Keyword arguments

- `savename = "test.mp4"`: name of video output file
- `framerate = 4`: framerate of video output
- `markers`: dictionary mapping attractor ids to markers they should be plotted as;
  attractors found are scatter-plotted with the given markers and additional
  trajectories are integrated and plotted on top for better visual coverage
- `markersize = 10`
- `Δt, T`: propagated to `trajectory` for evolving an initial condition sampled
  from an attractor
- Also all [common plotting keywords](@ref).
"""
function animate_attractors_continuation end
export animate_attractors_continuation
