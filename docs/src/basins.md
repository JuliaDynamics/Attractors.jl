# Basins of Attraction
This page provides several functions related to the basins of attraction and their boundaries. It requires you to have first understood the [Finding Attractors](@ref) page.


## Basins of attraction
Calculating basins of attraction, or their state space fractions, can be done with the functions:
- [`basins_fractions`](@ref)
- [`basins_of_attraction`](@ref).

```@docs
basins_fractions
extract_attractors
basins_of_attraction
convergence_and_basins_of_attraction
iterations_to_converge
statespace_sampler
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
