# Visualization utilities

In this page we document several plotting utility functions that have been created to make the visualization of the output of Attractors.jl seamless. See the examples page for usage of all these plotting functions.

Note that most functions have an out-of-place and an in-place form, the in-place form always taking as a first input a pre-initialized `Axis` to plot in while the out-of-place creates and returns a new figure object.

E.g.,

```julia
heatmap_basins_attractors(grid, basins, attractors; kwargs...)
heatmap_basins_attractors!(ax, grid, basins, attractors; kwargs...)
```

## [Common plotting keywords](@id common_plot_kwargs)
Common keywords for plotting functions in Attractors.jl are:

- `ukeys`: the basin ids (unique keys, vector of integers) to use. By default all existing keys are used.
- `access = [1, 2]`: indices of which dimensions of an attractor to select and visualize in a two-dimensional plot.
  Only these ids will be visualized. By default all are used.
- `colors`: a dictionary mapping basin ids (i.e., including the `-1` key) to a color. By default the JuliaDynamics colorscheme is used if less than 7 ids are present, otherwise random colors from the `:darktest` colormap.
- `markers`: dictionary mapping attractor ids to markers they should be plotted as
- `labels = Dict(ukeys .=> ukeys)`: how to label each attractor.
- `add_legend = length(ukeys) < 7`: whether to add a legend mapping colors to labels.

## Basins related

```@docs
heatmap_basins_attractors
```

## Continuation related

```@docs
plot_basins_curves
plot_attractors_curves
plot_basins_attractors_curves
```

## Video output

```@docs
animate_attractors_continuation
```