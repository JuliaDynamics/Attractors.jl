# Attractors.jl

[![](https://img.shields.io/badge/docs-dev-lightblue.svg)](https://JuliaDynamics.github.io/Attractors.jl/dev)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaDynamics.github.io/Attractors.jl/stable)
[![Paper](https://img.shields.io/badge/Cite-arXiv:2304.12786-purple)](https://arxiv.org/abs/2304.12786)
[![CI](https://github.com/JuliaDynamics/Attractors.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/Attractors.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/JuliaDynamics/Attractors.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaDynamics/Attractors.jl)
[![Package Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/Attractors)](https://pkgs.genieframework.com?packages=Attractors)

A Julia module for finding attractors of dynamical systems,
their basins and their boundaries, fractal properties of the boundaries,
as well as continuing attractors and their basins across parameters.
It can be used as a standalone package, or as part of
[DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystems.jl/dev/).

To install it, run `import Pkg; Pkg.add("Attractors")`.

All further information is provided in the documentation, which you can either find [online](https://juliadynamics.github.io/Attractors.jl/dev/) or build locally by running the `docs/make.jl` file.

_Previously, Attractors.jl was part of ChaosTools.jl_