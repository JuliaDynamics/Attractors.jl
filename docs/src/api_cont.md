# Global continuation

```@docs
global_continuation
GlobalContinuationAlgorithm
GlobalContinuationOutput
continuation_series
```



## General seeding-based continuation (ASCM)

```@docs
AttractorSeedContinueMatch
```

## Recurrences continuation

```@docs
RecurrencesFindAndMatch
```

## Grouping continuation (FGAP)

```@docs
FeaturizeGroupAcrossParameter
```

## Aggregating attractors

```@docs
aggregate_attractors
aggregate_continuation
aggregate_fractions
```

## Dict utils

```@docs
unique_keys
swap_dict_keys!
next_free_id
```

# Matching attractors

Matching attractors follow an extendable interface based on [`IDMatcher`](@ref).
The available matchers are:

```@docs
MatchBySSSetDistance
MatchByBasinEnclosure
MatchByBasinOverlap
DontMatch
```

## Matching interface

```@docs
IDMatcher
matching_map
matching_map!
match_sequentially!
```

## Low-level distance functions

```@docs
Centroid
Hausdorff
StrictlyMinimumDistance
set_distance
setsofsets_distances
```

