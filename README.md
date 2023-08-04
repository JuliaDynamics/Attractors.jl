# Attractors.jl

[![](https://img.shields.io/badge/docs-dev-lightblue.svg)](https://JuliaDynamics.github.io/Attractors.jl/dev)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaDynamics.github.io/Attractors.jl/stable)
[![Paper](https://img.shields.io/badge/Cite-DOI:10.1063/5.0159675-purple)](https://arxiv.org/abs/2304.12786)
[![CI](https://github.com/JuliaDynamics/Attractors.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/Attractors.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/JuliaDynamics/Attractors.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaDynamics/Attractors.jl)
[![Package Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/Attractors)](https://pkgs.genieframework.com?packages=Attractors)

A Julia module for

- finding attractors of arbitrary dynamical systems
- finding their basins of attraction or the state space fractions of the basins
- "continuing" the attractors and their basins over a parameter range
- finding the basin boundaries and analyzing their fractal properties
- tipping points related functionality for systems with known dynamic rule
- and more!

It can be used as a standalone package, or as part of
[DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystems.jl/stable/).

To install it, run `import Pkg; Pkg.add("Attractors")`.

All further information is provided in the documentation, which you can either find [online](https://juliadynamics.github.io/Attractors.jl/stable/) or build locally by running the `docs/make.jl` file.

_Previously, Attractors.jl was part of ChaosTools.jl_