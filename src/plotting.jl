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
- `markersize = 10`
- `Δt, T`: propagated to `trajectory` for evolving an initial condition sampled
  from an attractor
- Also all [common plotting keywords](@ref).
"""
function animate_attractors_continuation end
export animate_attractors_continuation

"""
    plot_basins_curves(fractions_curves, prange = 1:length(); kwargs...)

Plot the fractions of basins of attraction versus a parameter range,
i.e., visualize the output of [`continuation`](@ref).

## Keyword arguments

- `style = :band`: how to visualize the basin fractions. Choices are
  `:cumulative` for a band plot with cumulative sum = 1 or `:lines` for a scatterline
  plot of each basin fraction
- `separatorwidth = 1, separatorcolor = "white"`: adds a line separating the fractions
  if the style is `:band`
- `axislegend_kwargs = (position = :lt,)`: propagated to `axislegend` if a legend is added
- `series_kwargs = NamedTuple()`: propagated to the band or scatterline plot
- Also all [common plotting keywords](@ref).
"""
function plot_basins_curves end
export plot_basins_curves, plot_basins_curves!