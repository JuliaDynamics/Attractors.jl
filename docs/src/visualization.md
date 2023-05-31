# Visualization utilities

In this page we document several plotting utility functions that have been created to make the visualization of the output of Attractors.jl seamless. See the examples page for usage of all these plotting functions.

Note that all functions have an out-of-place and an in-place form, the in-place form always taking as a first input a pre-initialized `Axis` to plot in.

E.g.,

```julia
heatmap_basins_attractors(grid, basins, attractors; kwargs...)
heatmap_basins_attractors!(ax, grid, basins, attractors; kwargs...)
```

## Common plotting keywords
Common keywords for plotting functions in Attractors.jl are:

- `ukeys`: The keys (attractor ids, vector of integers) to use.
  Only these ids will be visualized. By default all are used.
- `colors`: A dictionary mapping attractor id to a color.
  By default the JuliaDynamics colorscheme is used if less than 7 ids are present, otherwise
  random colors form the `:darktest` colormap.
- `labels = Dict(ukeys .=> ukeys)`: How to label each attractor.
- `add_legend = length(ukeys) < 7`: Whether to add a legend mapping colors to labels.

## Basins related

## Continuation related