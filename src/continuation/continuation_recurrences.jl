export RecurrencesFindAndMatch, RAFM
import ProgressMeter
using Random: MersenneTwister

# The recurrences based distance is rather flexible because it works
# in two independent steps: it first finds attractors and then matches them.

"""
    RecurrencesFindAndMatch <: AttractorsBasinsContinuation
    RecurrencesFindAndMatch(mapper::AttractorsViaRecurrences; kwargs...)

A method for [`continuation`](@ref) as in [^Datseris2023] that is based on the
recurrences-based algorithmf or finding attractors ([`AttractorsViaRecurrences`](@ref))
and the "matching attractors" functionality offered by [`match_attractor_ids!`](@ref).

You can use `RAFM` as an alias.

## Description

At the first parameter slice attractors and their fractions are found as described in the
[`AttractorsViaRecurrences`](@ref) mapper using recurrences in state space.
At each subsequent parameter slice,
new attractors are found by seeding initial conditions from the previously found
attractors and then running these initial conditions through the recurrences algorithm
of the `mapper`. Seeding initial conditions close to previous attractors accelerates
the main bottleneck of [`AttractorsViaRecurrences`](@ref), which is finding the attractors.

After the attractors are found, their fractions are computed by sampling new random initial
(using the provided `sampler` in [`continuation`](@ref)) and mapping them to attractors
using the [`AttractorsViaRecurrences`](@ref) mapper.

Then, the newly found attractors (and their fractions) are "matched" to the previous ones.
I.e., their _IDs_ are changed, accoring to the [`match_attractor_ids!`](@ref) function.
Typically, the matching process matches attractor IDs that are closest in state space
distance, but mroe options are possible, see.
] attractors across parameters are matched so that they have
the same ID. The matching process is based on distances between attractors.
The function that computes these distances is
[`setsofsets_distances`](@ref) and the matching function
is [`match_attractor_ids!`](@ref) (please read those docstrings as well).

This process continues until all parameter values are exhausted and for each parameter
value the attractors and their fractions are found.


At each parameter slice beyond the first, the new
attractors are matched to the previous attractors found in the previous parameter value
by a direct call to the [`match_attractor_ids!`](@ref) function. Hence, the matching
of attractors here works "slice by slice" on the parameter axis and the attractors
that are closest to each other (in state space, but for two different parameter values)
get assigned the same label.

## Keyword arguments
- `distance, threshold`: propagated to [`match_attractor_ids!`](@ref).
- `info_extraction = identity`: A function that takes as an input an attractor (`StateSpaceSet`)
  and outputs whatever information should be stored. It is used to return the
  `attractors_info` in [`rsc`](@ref).
- `seeds_from_attractor`: A function that takes as an input an attractor and returns
  an iterator of initial conditions to be seeded from the attractor for the next
  parameter slice. By default, we sample some points from existing attractors according
  to how many points the attractors themselves contain. A maximum of `10` seeds is done
  per attractor.

[^Datseris2023]: Datseris, Rossi & Wagemakers 2023: Framework for global stability analysis
"""
struct RecurrencesFindAndMatch{A, M, R<:Real, S, E} <: AttractorsBasinsContinuation
    mapper::A
    distance::M
    threshold::R
    seeds_from_attractor::S
    info_extraction::E
end

"Alias for [`RecurrencesFindAndMatch`](@ref)"
const RAFM = RecurrencesFindAndMatch

function RecurrencesFindAndMatch(
        mapper::AttractorsViaRecurrences; distance = Centroid(),
        threshold = Inf, seeds_from_attractor = _default_seeding_process,
        info_extraction = identity
    )
    return RecurrencesFindAndMatch(
        mapper, distance, threshold, seeds_from_attractor, info_extraction
    )
end

function _default_seeding_process(attractor::AbstractStateSpaceSet; rng = MersenneTwister(1))
    max_possible_seeds = 10
    seeds = round(Int, log(10, length(attractor)))
    seeds = clamp(seeds, 1, max_possible_seeds)
    return (rand(rng, attractor.data) for _ in 1:seeds)
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

    (; mapper, distance, threshold) = rsc
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
        fs = basins_fractions(mapper, ics;
            additional_fs = seeded_fs, show_progress = false, N = samples_per_parameter
        )
        current_attractors = mapper.bsn_nfo.attractors
        if !isempty(current_attractors) && !isempty(prev_attractors)
            # If there are any attractors,
            # match with previous attractors before storing anything!
            rmap = match_attractor_ids!(
                current_attractors, prev_attractors; distance, threshold
            )
            swap_dict_keys!(fs, rmap)
        end
        # Then do the remaining setup for storing and next step
        push!(fractions_curves, fs)
        push!(attractors_info, get_info(current_attractors))
        overwrite_dict!(prev_attractors, current_attractors)
        ProgressMeter.next!(progress; showvalues = [("previous parameter", p),])
    end
    # Normalize to smaller available integers for user convenience
    rmap = retract_keys_to_consecutive(fractions_curves)
    for (da, df) in zip(attractors_info, fractions_curves)
        swap_dict_keys!(da, rmap)
        swap_dict_keys!(df, rmap)
    end
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
    # mapper.bsn_nfo.visited_cell = 4
    # because we want the next attractor to be labelled differently in case
    # it doesn't actually match to any of the new ones
    return
end

function _ics_from_grid(rsc::RecurrencesFindAndMatch)
    return _ics_from_grid(rsc.mapper.grid)
end

function _ics_from_grid(grid::Tuple)
    sampler, = statespace_sampler(min_bounds = minimum.(grid), max_bounds = maximum.(grid))
    return sampler
end
