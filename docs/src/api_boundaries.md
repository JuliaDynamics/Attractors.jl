# Basin boundaries

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

The edge tracking algorithm allows to locate and construct so-called edge states embedded in the basin boundary separating different basins of attraction. These could be saddle points, unstable periodic orbits or chaotic saddles. The general idea is that these sets can be found because they act as attractors when restricting to the basin boundary. Another technique to get a pseudo trajectory close to a saddle is the stagger-and-step method that requires little information on the dynamical system.

Two functions provide functionalities for this: [`edgetracking`](@ref) and [`stagger_and_step`](@ref).

```@docs
edgetracking
EdgeTrackingResults
bisect_to_edge
```

```@docs
stagger_and_step
```
