# API


## Finding attractors

Attractors.jl defines a generic interface for finding attractors of dynamical systems. One first decides the instance of [`DynamicalSystem`](@ref) they need. Then, an instance of [`AttractorMapper`](@ref) is created from this dynamical system. This `mapper` instance can be used to compute e.g., [`basins_of_attraction`](@ref), and the output can be further analyzed to get e.g., the [`basin_entropy`](@ref).

```@docs
AttractorMapper
```


### Recurrences

```@docs
AttractorsViaRecurrences
automatic_Δt_basins
SubdivisionBasedGrid
subdivision_based_grid
```

### Proximity
```@docs
AttractorsViaProximity
```

### Featurizing
```@docs
AttractorsViaFeaturizing
```


## Grouping configurations

Grouping configurations that can be given to [`AttractorsViaFeaturizing`](@ref)
are part of a generic and extendable interface based on the [`group_features`](@ref)
function.
The grouping configuration sets how the features describing the trajectories will be grouped together.
Nevertheless, this grouping infrastructure can also be used and extended completely independently of finding attractors of dynamical systems!

### Grouping interface

```@docs
group_features
GroupingConfig
```

### Grouping types

```@docs
GroupViaClustering
GroupViaHistogram
GroupViaNearestFeature
GroupViaPairwiseComparison
```

### Grouping utils

```@docs
extract_features
```


## Basins of attraction
Calculating basins of attraction, or their state space fractions, can be done with the functions:
- [`basins_fractions`](@ref)
- [`basins_of_attraction`](@ref).

```@docs
basins_fractions
extract_attractors
basins_of_attraction
statespace_sampler
```

## Convergence times

```@docs
convergence_and_basins_fractions
convergence_and_basins_of_attraction
convergence_time
```

## Final state sensitivity / fractal boundaries
Several functions are provided related with analyzing the fractality of the boundaries of the basins of attraction:

- [`basins_fractal_dimension`](@ref)
- [`basin_entropy`](@ref)
- [`basins_fractal_test`](@ref)
- [`uncertainty_exponent`](@ref)
- [`test_wada_merge`](@ref)

```@docs
basins_fractal_dimension
basin_entropy
basins_fractal_test
uncertainty_exponent
test_wada_merge
```

## Edge tracking and edge states
The edge tracking algorithm allows to locate and construct so-called edge states (also referred to as *Melancholia states*) embedded in the basin boundary separating different basins of attraction. These could be saddle points, unstable periodic orbits or chaotic saddles. The general idea is that these sets can be found because they act as attractors when restricting to the basin boundary.

```@docs
edgetracking
EdgeTrackingResults
bisect_to_edge
```

## Tipping points
This page discusses functionality related with tipping points in dynamical systems with known rule. If instead you are interested in identifying tipping points in measured timeseries, have a look at [TransitionIndicators.jl](https://github.com/JuliaDynamics/TransitionIndicators.jl).

```@docs
tipping_probabilities
```

## Minimal Fatal Shock
The algorithm to find minimal perturbation for arbitrary initial condition `u0` which will kick the system into different from the current basin.
```@docs
minimal_fatal_shock
MFSBlackBoxOptim
MFSBruteForce
```

## Continuation API

```@docs
continuation
```

### Recurrences continuation

```@docs
RecurrencesFindAndMatch
```

### Aggregating attractors and fractions

```@docs
aggregate_attractor_fractions
```

### Grouping continuation

```@docs
FeaturizeGroupAcrossParameter
```

## Matching continuation output

One of the most novel, and most powerful features of Attractors.jl is its
ability to "match" the continuation result: i.e., functionality that decides
which attractor gets which label at each parameter value, so that there is a
"continuity" across the parameter axis.

The matching process is entirely orthogonal to the continuation, something
completely novel in continuation software. This means, that if you don't like
the way the matching worked in the first time you estimated the attractors and their
basins, you don't have to re-compute them! You can simply relabel them using the
[`match_continuation!`](@ref) function!

```@docs
match_statespacesets!
Centroid
Hausdorff
StrictlyMinimumDistance
replacement_map
set_distance
setsofsets_distances
match_continuation!
match_basins_ids!
```

## Visualization utilities

Several plotting utility functions have been created to make the visualization of the output of Attractors.jl seamless. See the examples page for usage of all these plotting functions.

Note that all functions have an out-of-place and an in-place form, the in-place form always taking as a first input a pre-initialized `Axis` to plot in while the out-of-place creates and returns a new figure object.

E.g.,

```julia
heatmap_basins_attractors(grid, basins, attractors; kwargs...)
heatmap_basins_attractors!(ax, grid, basins, attractors; kwargs...)
```

### [Common plotting keywords](@id common_plot_kwargs)
Common keywords for plotting functions in Attractors.jl are:

- `ukeys`: the basin ids (unique keys, vector of integers) to use. By default all existing keys are used.
- `access = [1, 2]`: indices of which dimensions of an attractor to select and visualize in a two-dimensional plot (as in [`animate_attractors_continuation`](@ref)).
- `colors`: a dictionary mapping basin ids (i.e., including the `-1` key) to a color. By default the JuliaDynamics colorscheme is used if less than 7 ids are present, otherwise random colors from the `:darktest` colormap.
- `markers`: dictionary mapping attractor ids to markers they should be plotted as
- `labels = Dict(ukeys .=> ukeys)`: how to label each attractor.
- `add_legend = length(ukeys) < 7`: whether to add a legend mapping colors to labels.

### Basins related

```@docs
heatmap_basins_attractors
shaded_basins_heatmap
```

### Continuation related

```@docs
plot_basins_curves
plot_attractors_curves
plot_basins_attractors_curves
```

### Video output

```@docs
animate_attractors_continuation
```