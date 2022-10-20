# Basins of Attraction
This page provides several functions related to the basins of attraction and their boundaries. It requires you to have first understood the [Finding Attractors](@ref) page.


## Basins of attraction
Calculating basins of attraction, or their state space fractions, can be done with the functions [`basins_fractions`](@ref) and [`basins_of_attraction`](@ref). See also the convenience [`statespace_sampler`](@ref) and [`match_attractors`](@ref).

```@docs
basins_fractions
basins_of_attraction
statespace_sampler
match_attractors!
```

## Final state sensitivity / fractal boundaries
Several functions are provided related with analyzing the fractality of the boundaries of the basins of attraction:

* [`basins_fractal_dimension`](@ref)
* [`basin_entropy`](@ref)
* [`basins_fractal_test`](@ref)
* [`uncertainty_exponent`](@ref)

```@docs
basins_fractal_dimension
basin_entropy
basins_fractal_test
uncertainty_exponent
```

## Tipping points
```@docs
tipping_probabilities
```