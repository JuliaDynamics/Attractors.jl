# API


## Finding attractors

Attractors.jl defines a generic interface for finding attractors of dynamical systems. One first decides the instance of [`DynamicalSystem`](@ref) they need. Then, an instance of [`AttractorMapper`](@ref) is created from this dynamical system. This `mapper` instance can be used to compute e.g., [`basins_of_attraction`](@ref), and the output can be further analyzed to get e.g., the [`basin_entropy`](@ref).

```@docs
AttractorMapper
extract_attractors
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
The basins of attraction are often represented as an array or vector. We also provide a convenient extendable structure
that contains the basins themselves, attractors, and the domains on which the basins are defined. All standard basin-related 
functions are compatible with this alternate representation.

```@docs
BasinsOfAttraction
ArrayBasinsOfAttraction
SampledBasinsOfAttraction
```

Calculating basins of attraction, or their state space fractions, can be done with the functions:
- [`basins_fractions`](@ref)
- [`basins_of_attraction`](@ref)

```@docs
basins_fractions
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

## Edge tracking, edge states and chaotic saddles
The edge tracking algorithm allows to locate and construct so-called edge states (also referred to as *Melancholia states*) embedded in the basin boundary separating different basins of attraction. These could be saddle points, unstable periodic orbits or chaotic saddles. The general idea is that these sets can be found because they act as attractors when restricting to the basin boundary. Another technique to get a pseudo trajectory close to a saddle is the stagger-and-step method that requires little information on the dynamical system.

```@docs
edgetracking
EdgeTrackingResults
bisect_to_edge
stagger_and_step
```

## Tipping points
This page discusses functionality related with tipping points in dynamical systems with known rule. If instead you are interested in identifying tipping points in measured timeseries, have a look at [TransitionIndicators.jl](https://github.com/JuliaDynamics/TransitionIndicators.jl).

```@docs
tipping_probabilities
```

## Minimal Critical Shock

The algorithm to find minimal perturbation for arbitrary initial condition `u0` which will kick the system into different from the current basin.

```@docs
minimal_critical_shock
MCSBlackBoxOptim
MCSBruteForce
```


## Global continuation

```@docs
global_continuation
GlobalContinuationAlgorithm
continuation_series
```

### General seeding-based continuation

```@docs
AttractorSeedContinueMatch
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

## Matching attractors

Matching attractors follow an extendable interface based on [`IDMatcher`](@ref).
The available matchers are:

```@docs
MatchBySSSetDistance
MatchByBasinEnclosure
MatchByBasinOverlap
```

### Matching interface

```@docs
IDMatcher
matching_map
matching_map!
match_sequentially!
```

### Low-level distance functions

```@docs
Centroid
Hausdorff
StrictlyMinimumDistance
set_distance
setsofsets_distances
```

### Dict utils

```@docs
unique_keys
swap_dict_keys!
next_free_id
```

## Stability measures

```@docs
StabilityMeasuresAccumulator
finalize_accumulator
stability_measures_along_continuation
```

## Visualization utilities

Several plotting utility functions have been created to make the visualization of the output of Attractors.jl seamless. See the examples page for usage of all these plotting functions.

Note that all functions have an out-of-place and an in-place form, the in-place form always taking as a first input a pre-initialized `Axis` to plot in while the out-of-place creates and returns a new figure object.

E.g.,

```julia
fig = heatmap_basins_attractors(grid, basins, attractors; kwargs...)
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
plot_attractors
heatmap_basins_attractors
shaded_basins_heatmap
```

### Continuation related

```@docs
plot_basins_curves
plot_attractors_curves
plot_basins_attractors_curves
plot_continuation_curves
```

### Video output

```@docs
animate_attractors_continuation
```
