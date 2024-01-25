# Attractor & Basins Continuation

## A new kind of continuation

If you have heard before the word "continuation", then you are likely aware of the **traditional continuation-based bifurcation analysis (CBA)** offered by many software, such as AUTO, MatCont, and in Julia [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl). Here we offer a completely different kind of continuation called **attractors & basins continuation**.

A direct comparison of the two approaches is not truly possible, because they do different things. The traditional linearized continuation analysis continues the curves of individual fixed points across the joint state-parameter space. The attractor and basins continuation first finds all attractors at all parameter values and then _matches_ appropriately similar attractors across different parameters, giving the illusion of continuing them individually. Additionally, the curves of stable fixed points in the joint parameter space is only a small by-product of the attractor basins continuation, and the main information is the basin fractions and how these change in the parameter space.

A more detailed comparison between these two fundamentally different approaches in is given in high detail in our [paper](https://github.com/JuliaDynamics/Attractors.jl/blob/main/CITATION.bib).

## Continuation API

```@docs
continuation
```

## Recurrences continuation (best)

```@docs
RecurrencesFindAndMatch
```

## Matching attractors

```@docs
match_statespacesets!
Centroid
Hausdorff
StrictlyMinimumDistance
replacement_map
set_distance
setsofsets_distances
match_continuation!
match_basins_ids!
```

## Aggregating attractors and fractions

```@docs
aggregate_attractor_fractions
```

## Grouping continuation

```@docs
FeaturizeGroupAcrossParameter
```
