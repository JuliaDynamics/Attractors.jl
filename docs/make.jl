cd(@__DIR__)
using Attractors
using Attractors.DynamicalSystemsBase
using Attractors.StateSpaceSets

import Documenter
using Literate

# literate the tutorial
Literate.markdown(
    joinpath(@__DIR__, "src", "tutorial.jl"), joinpath(@__DIR__, "src");
    credit = false
)
# and the comparison with BFKIT
Literate.markdown(
    joinpath(@__DIR__, "src", "bfkit_comparison.jl"), joinpath(@__DIR__, "src");
    credit = false,
    flavor = Literate.CommonMarkFlavor()
)


pages = [
    "index.md",
    "tutorial.md",
    "api.md",
    "examples.md",
    Documenter.hide("bfkit_comparison.md"),
    "references.md",
    Documenter.hide("recurrences_animation.md"),
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

build_docs_with_style(pages, Attractors, StateSpaceSets;
    expandfirst = ["index.md"], bib, warnonly = [:doctest, :missing_docs, :cross_references],
)
