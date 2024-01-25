export RecurrencesFindAndMatch, RAFM
import ProgressMeter
using Random: MersenneTwister

"""
    RecurrencesFindAndMatch <: AttractorsBasinsContinuation
    RecurrencesFindAndMatch(mapper::AttractorsViaRecurrences; kwargs...)

A method for [`continuation`](@ref) as in [Datseris2023](@cite) that is based on the
recurrences algorithm for finding attractors ([`AttractorsViaRecurrences`](@ref))
and the "matching attractors" functionality offered by [`match_continuation!`](@ref).

You can use `RAFM` as an alias.

## Keyword arguments

- `distance = Centroid(), threshold = Inf, use_vanished = !isinf(threshold)`:
  propagated to [`match_continuation!`](@ref).
- `info_extraction = identity`: A function that takes as an input an attractor (`StateSpaceSet`)
  and outputs whatever information should be stored. It is used to return the
  `attractors_info` in [`continuation`](@ref). Note that the same attractors that
  are stored in `attractors_info` are also used to perform the matching in
  [`match_continuation!`](@ref), hence this keyword should be altered with care.
- `seeds_from_attractor`: A function that takes as an input an attractor and returns
  an iterator of initial conditions to be seeded from the attractor for the next
  parameter slice. By default, we sample only the first stored point on the attractor.

## Description

At the first parameter slice of the continuation process, attractors and their fractions
are found as described in the [`AttractorsViaRecurrences`](@ref) mapper using recurrences
in state space. At each subsequent parameter slice,
new attractors are found by seeding initial conditions from the previously found
attractors and then running these initial conditions through the recurrences algorithm
of the `mapper`. Seeding initial conditions close to previous attractors accelerates
the main bottleneck of [`AttractorsViaRecurrences`](@ref), which is finding the attractors.

After the special initial conditions are mapped to attractors, attractor basin fractions
are computed by sampling random initial conditions using the provided `sampler` in
[`continuation`](@ref)) and mapping them to attractors
using the [`AttractorsViaRecurrences`](@ref) mapper.
I.e., exactly as in [`basins_fractions`](@ref).
Naturally, during this step new attractors may be found, besides those found
using the "seeding from previous attractors".
Once the basins fractions are computed,
the parameter is incremented again and we perform the steps as before.

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
the docstrings of [`match_statespacesets!!`](@ref) first and then [`match_continuation!`](@ref).
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
struct RecurrencesFindAndMatch{A, M, R<:Real, S, E} <: AttractorsBasinsContinuation
    mapper::A
    distance::M
    threshold::R
    use_vanished::Bool
    seeds_from_attractor::S
    info_extraction::E
end

"Alias for [`RecurrencesFindAndMatch`](@ref)."
const RAFM = RecurrencesFindAndMatch

function RecurrencesFindAndMatch(
        mapper::AttractorsViaRecurrences; distance = Centroid(),
        threshold = Inf, use_vanished = !isinf(threshold),
        seeds_from_attractor = _default_seeding_process,
        info_extraction = identity
    )
    return RecurrencesFindAndMatch(
        mapper, distance, threshold, use_vanished, seeds_from_attractor, info_extraction
    )
end

# TODO: This is currently not used, and not sure if it has to be.
function _default_seeding_process_10(attractor::AbstractStateSpaceSet; rng = MersenneTwister(1))
    max_possible_seeds = 10
    seeds = round(Int, log(10, length(attractor)))
    seeds = clamp(seeds, 1, max_possible_seeds)
    return (rand(rng, vec(attractor)) for _ in 1:seeds)
end

# This is the one used
function _default_seeding_process(attractor::AbstractStateSpaceSet)
    return (attractor[1],) # must be iterable
end

function continuation(
        rsc::RecurrencesFindAndMatch,
        prange, pidx, ics = _ics_from_grid(rsc);
        samples_per_parameter = 100, show_progress = true,
    )
    if ics isa AbstractStateSpaceSet
        error("`ics` needs to be a function.")
    end
    progress = ProgressMeter.Progress(length(prange);
        desc="Continuating basins fractions:", enabled=show_progress
    )

    mapper = rsc.mapper
    reset!(mapper)
    # first parameter is run in isolation, as it has no prior to seed from
    set_parameter!(mapper.ds, pidx, prange[1])
    fs = basins_fractions(mapper, ics; show_progress = false, N = samples_per_parameter)
    # At each parmaeter `p`, a dictionary mapping attractor ID to fraction is created.
    fractions_curves = [fs]
    # Furthermore some info about the attractors is stored and returned
    prev_attractors = deepcopy(mapper.bsn_nfo.attractors)
    get_info = attractors -> Dict(
        k => rsc.info_extraction(att) for (k, att) in attractors
    )
    info = get_info(prev_attractors)
    attractors_info = [info]
    ProgressMeter.next!(progress; showvalues = [("previous parameter", prange[1]),])
    # Continue loop over all remaining parameters
    for p in prange[2:end]
        set_parameter!(mapper.ds, pidx, p)
        reset!(mapper)
        # Seed initial conditions from previous attractors
        # Notice that one of the things that happens here is some attractors have
        # really small basins. We find them with the seeding process here, but the
        # subsequent random sampling in `basins_fractions` doesn't. This leads to
        # having keys in `mapper.bsn_nfo.attractors` that do not exist in the computed
        # fractions. The fix is easy: we add the initial conditions mapped from
        # seeding to the fractions using an internal argument.
        seeded_fs = Dict{Int, Int}()
        for att in values(prev_attractors)
            for u0 in rsc.seeds_from_attractor(att)
                # We map the initial condition to an attractor, but we don't care
                # about which attractor we go to. This is just so that the internal
                # array of `AttractorsViaRecurrences` registers the attractors
                label = mapper(u0; show_progress = false)
                seeded_fs[label] = get(seeded_fs, label, 0) + 1
            end
        end
        # Now perform basin fractions estimation as normal, utilizing found attractors
        # (the function comes from attractor_mapping.jl)
        fs = basins_fractions(mapper, ics;
            additional_fs = seeded_fs, show_progress = false, N = samples_per_parameter
        )
        # We do not match attractors here; the matching is independent step done at the end
        current_attractors = mapper.bsn_nfo.attractors
        push!(fractions_curves, fs)
        push!(attractors_info, get_info(current_attractors))
        overwrite_dict!(prev_attractors, current_attractors)
        ProgressMeter.next!(progress; showvalues = [("previous parameter", p),])
    end
    # Match attractors (and basins)
    (; use_vanished, distance, threshold) = rsc
    match_continuation!(fractions_curves, attractors_info; distance, threshold, use_vanished)
    return fractions_curves, attractors_info
end

function reset!(mapper::AttractorsViaRecurrences)
    empty!(mapper.bsn_nfo.attractors)
    if mapper.bsn_nfo.basins isa Array
        mapper.bsn_nfo.basins .= 0
    else
        empty!(mapper.bsn_nfo.basins)
    end
    mapper.bsn_nfo.state = :att_search
    mapper.bsn_nfo.consecutive_match = 0
    mapper.bsn_nfo.consecutive_lost = 0
    mapper.bsn_nfo.prev_label = 0
    # notice that we do not reset the following:
    # mapper.bsn_nfo.current_att_label = 2
    # mapper.bsn_nfo.visited_cell_label = 4
    # because we want the next attractor to be labelled differently in case
    # it doesn't actually match to any of the new ones
    return
end

function _ics_from_grid(rsc::RecurrencesFindAndMatch)
    return _ics_from_grid(rsc.mapper.grid)
end

function _ics_from_grid(grid::Tuple)
    sampler, = statespace_sampler(grid)
    return sampler
end

