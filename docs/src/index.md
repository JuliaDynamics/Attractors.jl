# Attractors.jl

```@docs
Attractors
```

```@setup MAIN
using CairoMakie, Attractors
```

## Latest news

- **v1.36**: Three enhancements to `StabilityMeasuresAccumulator`:
  (1) local stability measures are now estimated for discrete-time systems too;
  (2) attractors can be aggregated before computing measures via `featurizer`/`group_config`
  keywords in `finalize_accumulator` and `stability_measures_along_continuation`;
  (3) `finalize_accumulator` now returns the (possibly aggregated) attractors as a second
  value so callers always have access to the correct IDs.
  See the CHANGELOG for details.
- Yet another groundbreaking feature has been implemented in Attractors.jl: estimation of extendably-many stability measures of a dynamical system automatically, with a single function call (also works along a global continuation)
- See the CHANGELOG.md (at the GitHub repo) for more!

## Getting started

Start by having a look at the [tutorial](@ref tutorial), after which you can
consult individual library functions in [API](@ref). Many more examples can be found
in the dedicated [examples](@ref examples) page.
