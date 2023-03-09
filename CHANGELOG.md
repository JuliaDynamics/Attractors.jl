# v1.2
- Add option `force_non_adaptive` to `AttractorsViaRecurrences`. This option is of instrumental importance when the grid is too fine or when limit cycle attractors are involved.

# v1.1

- `basins_fractions` no longer returns the attractors. Instead, a function `extract_attractors(mapper::AttractorMapper)` is provided that gives the attractors. This makes an overall more convenient experience that doesn't depend on the type of the initial conditions used.

# v1
- Added the `continuation` function and the two types `RecurrencesSeededContinuation` and `GroupAcrossParametersContinuation`

# v0.1.0
This is the first release. It continues from ChaosTools.jl v2.9, and hence, the comparison of attractor-related features is w.r.t to that version.

## Finding attractors
- New attractor mapping algorithm `AttractorsViaRecurrencesSparse` that uses sparse arrays to find attractors of arbitrarily high dimensional dynamical systems, eliminating the main drawback of `AttractorsViaRecurrences`.
- Clustering (used in `AttractorsViaFeaturizing`) has been completely overhauled. Now, a `ClusteringConfig` instances must be created and then passed on to `AttractorsViaFeaturizing`.
- `AttractorsViaFeaturizing` no longer has keywords about clustering.
- A new function `cluster_features` is exported to the public API.
- Multithreading is now an option in `AttractorsViaFeaturizing`. It is enabled by default.
- Added a new clause in automatic `ε` estimation in `AttractorsViaProximity` for when there is only a single attractor passed in by the user.

## Basin fractions continuation
- New function `continuation` that calculates basins fractions and how these change versus a parameter (given a continuation method)
- New basins fraction continuation method `RecurrencesSeededContinuation` that utilizes a brand new algorithm to continuate basins fractions of arbitrary systems.
- `match_attractor_ids!` has been fully overhauled to be more flexible, allow more ways to match, and also allow arbitrary user-defined ways to match.
- New function `match_basins_ids!` for matching the output of basins_of_attraction`.
- New exported functions `swap_dict_keys!, unique_keys, replacement_map` used in code that matches attractors and could be useful to front-end users.