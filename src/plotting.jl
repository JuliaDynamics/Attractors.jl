##########################################################################################
# Basins
##########################################################################################
"""
    heatmap_basins_attractors(grid, basins, attractors; kwargs...)

Plot a heatmap of found (2-dimensional) `basins` of attraction and corresponding
`attractors`, i.e., the output of [`basins_of_attraction`](@ref).

## Keyword arguments

- All the [common plotting keywords](@ref common_plot_kwargs).
"""
function heatmap_basins_attractors end
function heatmap_basins_attractors! end
export heatmap_basins_attractors, heatmap_basins_attractors!

"""
    shaded_basins_heatmap(grid, basins, attractors, iterations; kwargs...)

Plot a heatmap of found (2-dimensional) `basins` of attraction and corresponding
`attractors`. A matrix `iterations` with the same size of `basins` must be provided
to shade the color according to the value of this matrix. A small value corresponds
to a light color and a large value to a darker tone. This is useful to represent
the number of iterations taken for each initial condition to converge. See also
[`convergence_time`](@ref) to store this iteration number.

## Keyword arguments

- `show_attractors = true`: shows the attractor on plot
- `maxit = maximum(iterations)`: clip the values of `iterations` to
the value `maxit`. Useful when there are some very long iterations and keep the
range constrained to a given interval.
- All the [common plotting keywords](@ref common_plot_kwargs).
"""
function shaded_basins_heatmap end
function shaded_basins_heatmap! end
export shaded_basins_heatmap, shaded_basins_heatmap!

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

- `savename = "attracont.mp4"`: name of video output file.
- `framerate = 4`: framerate of video output.
- `Δt, T`: propagated to `trajectory` for evolving an initial condition sampled
  from an attractor.
- Also all [common plotting keywords](@ref common_plot_kwargs).
"""
function animate_attractors_continuation end
export animate_attractors_continuation

"""
    plot_basins_curves(fractions_curves, prange = 1:length(); kwargs...)

Plot the fractions of basins of attraction versus a parameter range,
i.e., visualize the output of [`continuation`](@ref).
See also [`plot_basins_attractors_curves`](@ref).

## Keyword arguments

- `style = :band`: how to visualize the basin fractions. Choices are
  `:band` for a band plot with cumulative sum = 1 or `:lines` for a lines
  plot of each basin fraction
- `separatorwidth = 1, separatorcolor = "white"`: adds a line separating the fractions
  if the style is `:band`
- `axislegend_kwargs = (position = :lt,)`: propagated to `axislegend` if a legend is added
- `series_kwargs = NamedTuple()`: propagated to the band or scatterline plot
- Also all [common plotting keywords](@ref common_plot_kwargs).
"""
function plot_basins_curves end
function plot_basins_curves! end
export plot_basins_curves, plot_basins_curves!

"""
    plot_attractors_curves(attractors_info, attractor_to_real, prange = 1:length(); kwargs...)

Same as in [`plot_basins_curves`](@ref) but visualizes the attractor dependence on
the parameter instead of their fraction.
The function `attractor_to_real` takes as input a `StateSpaceSet` (attractor)
and returns a real number so that it can be plotted versus the parameter axis.
See also [`plot_basins_attractors_curves`](@ref).

Same keywords as [`plot_basins_curves`](@ref common_plot_kwargs).
"""
function plot_attractors_curves end
function plot_attractors_curves! end
export plot_attractors_curves, plot_attractors_curves!

"""
    plot_basins_attractors_curves(
        fractions_curves, attractors_info, a2rs [, prange]
        kwargs...
    )

Convenience combination of [`plot_basins_curves`](@ref) and [`plot_attractors_curves`](@ref)
in a multi-panel plot that shares legend, colors, markers, etc.
This function allows `a2rs` to be a `Vector` of functions, each mapping
attractors into real numbers. Below the basins fractions plot, one additional
panel is created for each entry in `a2rs`.
`a2rs` can also be a single function, in which case only one panel is made.
"""
function plot_basins_attractors_curves end
function plot_basins_attractors_curves! end
export plot_basins_attractors_curves, plot_basins_attractors_curves!
