# Attractors.jl

[![docsdev](https://img.shields.io/badge/docs-dev-lightblue.svg)](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/attractors/dev/)
[![docsstable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/attractors/stable/)
[![Paper](https://img.shields.io/badge/Cite-DOI:10.1063/5.0159675-purple)](https://arxiv.org/abs/2304.12786)
[![CI](https://github.com/JuliaDynamics/Attractors.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/Attractors.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/JuliaDynamics/Attractors.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaDynamics/Attractors.jl)
[![Package Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FAttractors&query=total_requests&label=Downloads)](http://juliapkgstats.com/pkg/Attractors)

A Julia package for

- Finding all attractors, and all types of attractors, of arbitrary dynamical systems. An extendable interface allows for new algorithms for finding attractors.
- Finding their basins of attraction or the state space fractions of the basins.
  This includes finding exit basins (divergence to infinity).
- Analyzing nonlocal stability of attractors (also called global stability or  resilience).
- Performing **global continuation** of attractors and their basins (or other measures of stability), over a parameter range. Global continuation is a new, cutting-edge type of continuation that offers several advantages over traditional local continuation (AUTO, MatCont, BifurcationKit.jl, etc.), see the comparison in our docs.
- Finding the basin boundaries and edges states and analyzing their fractal properties.
- Tipping points related functionality for systems with known dynamic rule.
- And more!

It can be used as a standalone package, or as part of
[DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/stable/).

To install it, run `import Pkg; Pkg.add("Attractors")`.

All further information is provided in the documentation, which you can either find [online](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/attractors/stable/) or build locally by running the `docs/make.jl` file.

_Previously, Attractors.jl was part of ChaosTools.jl_