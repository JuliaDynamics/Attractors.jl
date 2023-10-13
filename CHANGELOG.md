# v1.13

- The algorithm of `AttractorsViaRecurrences` has been simplified a bit. The action of the keyword `mx_chk_loc_att` has been changed which may lead to different results in some usage cases.
- Additional benefit of this change is that incorrect algorithm behaviour can be caught eagerly. Now an error is thrown when we know in advance the algorithm will fail. (This also affects `continuation` with `RAFM`)
- The documentation of `AttractorsViaRecurrences` has been improved and clarified. Additionally a video illustrating algorithm behaviour has been added.
- A renaming of some of the options (keyword arguments) of `AttractorsViaRecurrences` has been done in line with the clarity increase of the algorithm. The following renames are in place and currently deprecated:
  - `mx_chk_fnd_att -> consecutive_recurrences`
  - `mx_chk_loc_att -> attractor_locate_steps`
  - `mx_chk_att -> consecutive_attractor_steps`
  - `mx_chk_hit_bas -> consecutive_basin_steps`
  - `mx_chk_lost -> consecutive_lost_steps`
  - `mx_chk_safety -> maximum_iterations`


# v1.12
- New algorithm `GroupViaPairwiseComparison` to group features in `AttractorsViaFeaturizing`. Simpler, typically faster and using less memory than DBSCAN, it can be useful in well-behaved systems.

# v1.11
- New algorithm `subdivision_based_grid`. Allows user to automatically constuct the grid which simulates subdivision into regions with different discretization levels in accordance with state space flow speed.

# v1.10
- Added support of irregular grids to `AttractorsViaRecurrences`, now user can provide ranges without fixed step along the same axis.

# v1.9
- Matching attractors during the continuation with `RAFM` has been improved and is done by the function `match_continuation!` which has two options regarding how to handle attractors of previous parameters that have vanished.

# v1.8
- New algorithm `minimal_fatal_shock` that finds the minimal perturbation for arbitrary initial condition `u0` which will kick the system into different from the current basin.

# v1.7.1
- Fixed issue where poincare map was not working with basins of attraction as intended (sampling points directly on the hyperplane)

# v1.7

- Some functions have been renamed for higher level of clarity (deprecations have been put in place):
  - `match_attractor_ids!` -> `match_statespacesets!`
  - `GroupAcrossParameter` -> `FeaturizeGroupAcrossParameter`.

# v1.6
- Attractors.jl moved to Julia 1.9+
- Plotting utility functions are now part of the API using package extensions. They are exported, documented, and used in the examples.

# v1.5
- Our pre-print regarding the global stability analysis framework offered by Attractors.jl is now online: https://arxiv.org/abs/2304.12786
- Package now has a CITATION.bib file.

# v1.4
- New function `rematch!` that can be used after `continuation` is called with `AttractorsViaRecurrences`, if the original matching performed was not ideal.

# v1.3
- More plotting functions have been added and plotting has been exposed as API. After Julia 1.9 it will also be documented and formally included in the built docs.

# v1.2
- Add option `force_non_adaptive` to `AttractorsViaRecurrences`. This option is of instrumental importance when the grid is too fine or when limit cycle attractors are involved.

# v1.1

- `basins_fractions` no longer returns the attractors. Instead, a function `extract_attractors(mapper::AttractorMapper)` is provided that gives the attractors. This makes an overall more convenient experience that doesn't depend on the type of the initial conditions used.

# v1
- Added the `continuation` function and the two types `RecurrencesFindAndMatch` and `GroupAcrossParameterContinuation`

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
- New basins fraction continuation method `RecurrencesFindAndMatch` that utilizes a brand new algorithm to continuate basins fractions of arbitrary systems.
- `match_attractor_ids!` has been fully overhauled to be more flexible, allow more ways to match, and also allow arbitrary user-defined ways to match.
- New function `match_basins_ids!` for matching the output of basins_of_attraction`.
- New exported functions `swap_dict_keys!, unique_keys, replacement_map` used in code that matches attractors and could be useful to front-end users.