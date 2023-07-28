cd(@__DIR__)
using Attractors
using Attractors.DynamicalSystemsBase
using Attractors.StateSpaceSets
using PredefinedDynamicalSystems

pages = [
    "index.md",
    "dynsysref.md",
    "attractors.md",
    "basins.md",
    "continuation.md",
    "visualization.md",
    "examples.md",
    "references.md",
]

import Downloads
Downloads.download(
    "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/build_docs_with_style.jl",
    joinpath(@__DIR__, "build_docs_with_style.jl")
)
include("build_docs_with_style.jl")

using DocumenterCitations

bib = CitationBibliography(
    joinpath(@__DIR__, "refs.bib");
    style=:authoryear
)

build_docs_with_style(pages, Attractors, DynamicalSystemsBase, StateSpaceSets;
    expandfirst = ["index.md"], bib,
)
