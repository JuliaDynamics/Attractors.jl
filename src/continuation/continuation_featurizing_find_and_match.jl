export FeaturizingFindAndMatch, FFAM
import ProgressMeter
using Random: MersenneTwister

"""
    FeaturizingFindAndMatch <: AttractorsBasinsContinuation
    FeaturizingFindAndMatch(mapper::AttractorsViaFeaturizing; kwargs...)

A method for [`continuation`](@ref) as in [Datseris2023](@cite) that uses the featurizing algorithm for finding attractors ([`AttractorsViaFeaturizing`](@ref))
and the "matching attractors" functionality offered by [`match_continuation!`](@ref). Based heavily on the `RecurrencesFindAndMatch`, which uses the recurrences algorithm ([`AttractorsViaRecurrences`](@ref)).

You can use `FFAM` as an alias.

## Keyword arguments

- `distance = Centroid(), threshold = Inf, use_vanished = !isinf(threshold)`:
  propagated to [`match_continuation!`](@ref).
- `info_extraction = identity`: A function that takes as an input an attractor (`StateSpaceSet`)
  and outputs whatever information should be stored. It is used to return the
  `attractors_info` in [`continuation`](@ref). Note that the same attractors that
  are stored in `attractors_info` are also used to perform the matching in
  [`match_continuation!`](@ref), hence this keyword should be altered with care.
- `seeds_from_attractor`: A function that takes as an input an attractor and the number of
  initial conditions to sample from the attractor, and returns a vector with the sampled initial
  conditions. By default, we just select the initial conditions randomly from the attractor.

## Description

At the first parameter slice of the continuation process, attractors and their fractions
are found as described in the [`AttractorsViaFeaturizing`](@ref) mapper identifying groups
of features. At each subsequent parameter slice, initial conditions are seeded from the
 previously found attractors. They are put together with the pre-generated `ics` given to
 `continuation`. Then, all of them are passed onto [`basins_fractions`](@ref) for
 [`AttractorsViaFeaturizing`](@ref), which featurizes and groups them to identify the
 attractors. 

Naturally, during this step new attractors may be found. Once the basins fractions are
computed, the parameter is incremented again and we perform the steps as before.

This process continues for all parameter values. After all parameters are exhausted,
the found attractors (and their fractions) are "matched" to the previous ones.
I.e., their _IDs are changed_, so that attractors that are "similar" to those at a
previous parameter get assigned the same ID.
Matching is done by the [`match_continuation!`](@ref) function and is an _orthogonal_
step. This means, that if you don't like the initial
outcome of the matching process, you may call [`match_continuation!`](@ref) again
on the outcome with different matching-related keywords.
You do not need to compute attractors and basins again!

Matching is a very sophisticated process that can be understood in detail by reading
the docstrings of [`match_statespacesets!`](@ref) first and then [`match_continuation!`](@ref).
Here is a short summary: attractors from previous and current parameter are matched
based on their "distance". By default this is distance in state space, but any measure of
"distance" may be used, such as the distance between Lyapunov spectra.
Matching prioritizes new->old pairs with smallest distance: once these are matched
the next available new->old pair with smallest distance is matched, until all new/old
attractors have been matched. The `threshold` keyword establishes that attractors with
distance > `threshold` do not get matched.
Additionally, use `use_vanished = true` if you want to include as matching candidates
attractors that have vanished during the continuation process.
"""
struct FeaturizingFindAndMatch{A, M, R<:Real, S, E} <: AttractorsBasinsContinuation
    mapper::A
    distance::M
    threshold::R
    seeds_from_attractor::S
    info_extraction::E
end

"Alias for [`FeaturizingFindAndMatch`](@ref)."
const FFAM = FeaturizingFindAndMatch

function FeaturizingFindAndMatch(
        mapper::AttractorsViaFeaturizing; distance = Centroid(),
        threshold = Inf, seeds_from_attractor = _default_seeding_process_featurizing,
        info_extraction = identity
    )
    return FeaturizingFindAndMatch(
        mapper, distance, threshold, seeds_from_attractor, info_extraction
    )
end

function _default_seeding_process_featurizing(attractor::AbstractStateSpaceSet, number_seeded_ics=10; rng = MersenneTwister(1))
    return [rand(rng, vec(attractor)) for _ in 1:number_seeded_ics] #might lead to repeated ics, which is intended for the continuation
end


function continuation(
        fam::FeaturizingFindAndMatch,
        prange, pidx, ics;
        samples_per_parameter = 100, show_progress = true, keep_track_maximum=true,
    )
    progress = ProgressMeter.Progress(length(prange);
        desc="Continuating basins fractions:", enabled=show_progress
    )

    if ics isa Function
        error("`ics` needs to be a Dataset.")
    end

    (; mapper, distance, threshold) = fam
    reset!(mapper)
    # first parameter is run in isolation, as it has no prior to seed from
    set_parameter!(mapper.ds, pidx, prange[1])
    fs, _ = basins_fractions(mapper, ics; show_progress = false)
    # At each parmaeter `p`, a dictionary mapping attractor ID to fraction is created.
    fractions_curves = [fs]
    # Furthermore some info about the attractors is stored and returned
    prev_attractors = deepcopy(extract_attractors(mapper))
    get_info = attractors -> Dict(
        k => fam.info_extraction(att) for (k, att) in attractors
    )
    info = get_info(prev_attractors)
    attractors_info = [info]
    ProgressMeter.next!(progress; showvalues = [("previous parameter", prange[1]),])
    alltime_maximum_key = maximum(keys(fs))
    # Continue loop over all remaining parameters
    for p in prange[2:end]
        set_parameter!(mapper.ds, pidx, p)
        reset!(mapper)
        
        # Collect ics from previous attractors to pass as additional ics to basins fractions (seeding).
        # To ensure that the clustering will identify them as clusters, we need to guarantee that there
        # are at least `min_neighbors` entries.
        num_additional_ics = typeof(mapper.group_config) <: GroupViaClustering ? 5*mapper.group_config.min_neighbors : 5
        additional_ics = Dataset(vcat(map(att-> 
            fam.seeds_from_attractor(att, num_additional_ics),
            values(prev_attractors))...)) #dataset with ics seeded from previous attractors
        
        # Now perform basin fractions estimation as normal, utilizing found attractors
        fs, _ = basins_fractions(mapper, ics;
            show_progress = false, additional_ics
        )
        
        current_attractors = extract_attractors(mapper)
        push!(fractions_curves, fs)
        push!(attractors_info, get_info(current_attractors))
        overwrite_dict!(prev_attractors, current_attractors)
        ProgressMeter.next!(progress; showvalues = [("previous parameter", p),])
    end
    # Match attractors (and basins)
    match_continuation!(fractions_curves, attractors_info; distance, threshold)
    return fractions_curves, attractors_info
end

function reset!(mapper::AttractorsViaFeaturizing)
    empty!(mapper.attractors)
end