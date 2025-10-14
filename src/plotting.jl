##########################################################################################
# Attractors
##########################################################################################
"""
    plot_attractors(attractors::Dict{Int, StateSpaceSet}; kwargs...)

Plot the attractors as a scatter plot.

## Keyword arguments

- All the [common plotting keywords](@ref common_plot_kwargs).
  Particularly important is the `access` keyword.
- `sckwargs = (strokewidth = 0.5, strokecolor = :black,)`: additional keywords
  propagated to the `Makie.scatter` function that plots the attractors.
"""
function plot_attractors end
function plot_attractors! end
export plot_attractors, plot_attractors!

##########################################################################################
# Basins
##########################################################################################
"""
    heatmap_basins_attractors(grid, basins, attractors; kwargs...)

Plot a heatmap of found (2-dimensional) `basins` of attraction and corresponding
`attractors`, i.e., the output of [`basins_of_attraction`](@ref).

## Keyword arguments

- All the [common plotting keywords](@ref common_plot_kwargs) and
  `sckwargs` as in [`plot_attractors`](@ref).
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
        ds::DynamicalSystem, attractors_cont, fractions_cont, pcurve; kw...
    )

Animate how the found system attractors and their corresponding basin fractions
change as the system parameter is increased. This function combines the input
and output of the [`global_continuation`](@ref) function into a video output.

The input dynamical system `ds` is used to evolve initial conditions sampled from the
found attractors, so that the attractors are better visualized.
`attractors_cont, fractions_cont` are the output of [`global_continuation`](@ref)
while `ds, pcurve` are the input to [`global_continuation`](@ref).

## Keyword arguments

- `savename = "attracont.mp4"`: name of video output file.
- `framerate = 4`: framerate of video output.
- `Δt, T`: propagated to `trajectory` for evolving an initial condition sampled
  from an attractor.
- Also all [common plotting keywords](@ref common_plot_kwargs).
- `figure, axis, fracaxis, legend`: named tuples propagated as keyword arguments to the
  creation of the `Figure`, the `Axis`, the "bar-like" axis containing the fractions,
  and the `axislegend` that adds the legend (if `add_legend = true`).
- `add_legend = true`: whether to display the axis legend.

## Continuation plot animation

In addition to animating the attractors scatterplot,
the plot of [`plot_attractors_curves`](@ref) can be added below the main attractor plot.
It is animated by a dashed line that moves as the parameter is changed.

This is enabled by providing a non-`nothing` keyword for `a2rs`,
and the behavior of plot is controlled by the following keywords:

- `a2rs`: a vector of functions that map attractors to real numbers.
- `prange = 1:length(attractors_cont)`: range of parameter values to plot on the x-axis
  of the additional plot(s). If given directly to `animate_attractors_continuation`
  it is obtained from there.
- `series_kwargs = NamedTuple()`: named tuple of arguments propagated
  to `Makie.scatterlines!` that plots the curves.
- `a2rs_ratio = 0.5`: the height ratio of the additional plot(s) relative to the main
  attractor plot.
- `a2rs_ylabels`: vector of y-axis labels for the additional plot(s).
- `parameter_name = "parameter"`. Name of the x-axis of the additional plot(s).
-  `vline_kwargs = (linestyle = :dash, linewidth = 3, color = "black")`: named tuple of arguments
  propagated to `Makie.vlines!` that plots the vertical dashed line moving
  as the parameter is changed.
"""
function animate_attractors_continuation end
export animate_attractors_continuation

"""
    plot_basins_curves(fractions_cont [, prange]; kw...)

Plot the fractions of basins of attraction versus a parameter range/curve,
i.e., visualize the output of [`global_continuation`](@ref).
See also [`plot_basins_attractors_curves`](@ref) and
[`plot_continuation_curves`](@ref).

## Keyword arguments

- `style = :band`: how to visualize the basin fractions. Choices are
  `:band` for a band plot with cumulative sum = 1 or `:lines` for a lines
  plot of each basin fraction
- `separatorwidth = 1, separatorcolor = "white"`: adds a line separating the fractions
  if the style is `:band`
- `filler = NaN`: filler value to use for basin fractions for attractor IDs that do not exist at a
  continuation step when the style is `:lines` (filler is always `0` for `:band` style).
- `axislegend_kwargs = (position = :lt,)`: propagated to `axislegend` if a legend is added
- `series_kwargs = NamedTuple()`: propagated to the band or scatterline plot
- Also all [common plotting keywords](@ref common_plot_kwargs).
"""
function plot_basins_curves end
function plot_basins_curves! end
export plot_basins_curves, plot_basins_curves!

"""
    plot_attractors_curves(attractors_cont, attractor_to_real [, prange]; kw...)

Same as in [`plot_basins_curves`](@ref) but visualize the attractor dependence on
the parameter(s) instead of their basin fraction.
The function `attractor_to_real` takes as input a `StateSpaceSet` (attractor)
and returns a real number so that it can be plotted versus the parameter axis.
See also [`plot_basins_attractors_curves`](@ref).

Same keywords as [`plot_continuation_curves`](@ref).
"""
function plot_attractors_curves end
function plot_attractors_curves! end
export plot_attractors_curves, plot_attractors_curves!

"""
    plot_continuation_curves(continuation_info [, prange]; kwargs...)

Same as in [`plot_basins_curves`](@ref) but visualize any arbitrary quantity characterizing
the continuation. Hence, the `continuation_info` is of exactly the same format as
`fractions_cont`: a vector of dictionaries, each dictionary mapping attractor IDs to real numbers.
`continuation_info` is meant to accompany `attractors_cont` in [`plot_attractors_curves`](@ref).

To produce `continuation_info` from `attractors_cont` you can do something like:

```julia
continuation_info = map(attractors_cont) do attractors
    Dict(k => f(A) for (k, A) in attractors)
end
```
with `f` your function of interest that returns a real number.

## Keyword arguments

- `series_kwargs`: named tuple of arguments propagated
  to `Makie.scatterlines!` that plots the curves.
- Also all [common plotting keywords](@ref common_plot_kwargs).

"""
function plot_continuation_curves end
function plot_continuation_curves! end
export plot_continuation_curves, plot_continuation_curves!

"""
    plot_basins_attractors_curves(
        fractions_cont, attractors_cont, a2rs [, prange]
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