export RecurrencesSeedingContinuation
import ProgressMeter
using Random: MersenneTwister

# The recurrences based method is rather flexible because it works
# in two independent steps: it first finds attractors and then matches them.
struct RecurrencesSeedingContinuation{A, M, S, E} <: BasinsFractionContinuation
    mapper::A
    method::M
    threshold::Float64
    seeds_from_attractor::S
    info_extraction::E
end

"""
    RecurrencesSeedingContinuation(mapper::AttractorsViaRecurrences; kwargs...)
A method for [`basins_fractions_continuation`](@ref).
It uses seeding of previous attractors to find new ones, which is the main performance
bottleneck. The method uses [`match_attractor_ids!`](@ref) to match attractors
as the system parameter is increased.

Will write more once we have the paper going.

## Keyword Arguments
- `method, threshold`: Given to [`match_attractor_ids!`](@ref) which is the function
  used to match attractors between each parameter slice.
- `info_extraction = identity`: A function that takes as an input an attractor (`Dataset`)
  and outputs whatever information should be stored. It is used to return the
  `attractors_info` in [`basins_fractions_continuation`](@ref).
- `seeds_from_attractor`: A function that takes as an input an attractor and returns
  an iterator of initial conditions to be seeded from the attractor for the next
  parameter slice. By default, we sample some points from existing attractors according
  to how many points the attractors themselves contain. A maximum of `10` seeds is done
  per attractor.
"""
function RecurrencesSeedingContinuation(
        mapper::AttractorsViaRecurrences; method = Centroid(),
        threshold = Inf, seeds_from_attractor = _default_seeding_process,
        info_extraction = identity
    )
    return RecurrencesSeedingContinuation(
        mapper, method, threshold, seeds_from_attractor, info_extraction
    )
end
@show 
function _default_seeding_process(attractor::AbstractDataset; rng = MersenneTwister(1))
    max_possible_seeds = 10
    seeds = round(Int, log(10, length(attractor)))
    seeds = clamp(seeds, 1, max_possible_seeds)
    return (rand(rng, attractor.data) for _ in 1:seeds)
end

function basins_fractions_continuation(
        continuation::RecurrencesSeedingContinuation,
        prange, pidx, ics = _ics_from_grid(continuation);
        samples_per_parameter = 100, show_progress = true,
    )
    # show_progress && @info "Starting basins fraction continuation."
    # show_progress && @info "p = $(prange[1])"
    progress = ProgressMeter.Progress(length(prange);
        desc="Continuating basins fractions:", enabled=show_progress
    )
    spp = samples_per_parameter
    (; mapper, method, threshold) = continuation
    # first parameter is run in isolation, as it has no prior to seed from
    labels, current_atts = get_attractors_and_labels(mapper, continuation, ics, pidx, prange[1], [], spp)
    prev_atts = deepcopy(current_atts)
    # At each parmaeter `p`, a dictionary mapping attractor ID to fraction is created.
    fractions_curves = [basins_fractions(labels)]
    get_info = attractors -> Dict(
        k => continuation.info_extraction(att) for (k, att) in attractors
    )
    info = get_info(prev_atts)
    attractors_info = [info]
    ProgressMeter.next!(progress; showvalues = [("previous parameter", prange[1]),])

    # Continue loop over all remaining parameters
    for p in prange[2:end]
        labels, current_atts = get_attractors_and_labels(mapper, continuation, ics, pidx, p, prev_atts, spp)
        fs = basins_fractions(labels)
        if !isempty(current_atts) && !isempty(prev_atts)
            # If there are any attractors,
            # match with previous attractors before storing anything!
            rmap = match_attractor_ids!(
                current_atts, prev_atts; method, threshold
            )
            swap_dict_keys!(fs, rmap)
        end
        # Then do the remaining setup for storing and next step
        push!(fractions_curves, fs)
        push!(attractors_info, get_info(current_atts))
        overwrite_dict!(prev_atts, current_atts)
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

function reset!(mapper::AttractorsViaRecurrences; reset_att = true)
    (reset_att == true) && empty!(mapper.bsn_nfo.attractors)
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

function _ics_from_grid(continuation::RecurrencesSeedingContinuation)
    return _ics_from_grid(continuation.mapper.grid)
end

function _ics_from_grid(grid::Tuple)
    sampler, = statespace_sampler(min_bounds = minimum.(grid), max_bounds = maximum.(grid))
    return sampler
end

# Seed initial conditions from previous attractors
# Notice that one of the things that happens here is some attractors have
# really small basins. We find them with the seeding process here, but the
# subsequent random sampling in `basins_fractions` doesn't. This leads to
# having keys in `mapper.bsn_nfo.attractors` that do not exist in the computed
# fractions. The fix is easy: we add the initial conditions mapped from
# seeding to the fractions using an internal argument.
function seed_attractors(prev_attractors, mapper, continuation)
    seeded_fs = Vector{Int}()
    for att in values(prev_attractors)
        for u0 in continuation.seeds_from_attractor(att)
            # We map the initial condition to an attractor, but we don't care
            # about which attractor we go to. This is just so that the internal
            # array of `AttractorsViaRecurrences` registers the attractors
            label = mapper(u0; show_progress = false)
            push!(seeded_fs, label)
        end
    end
    return seeded_fs
end

function get_attractors_and_labels(mapper, continuation, ics, pidx, p, prev_atts, N)
    set_parameter!(mapper.integ, pidx, p)
    reset!(mapper)
    labels = seed_attractors(prev_atts, mapper, continuation)
    # Now perform basin fractions estimation as normal, utilizing found attractors
    for i ∈ 1:N
        ic = _get_ic(ics, i)
        lab = mapper(ic; show_progress = false)
        push!(labels, lab)
    end
    return labels, mapper.bsn_nfo.attractors
end

_get_ic(ics::Function, i) = ics()
_get_ic(ics::AbstractDataset, i) = ics[i]

function basins_fractions_continuation_group(
        continuation::RecurrencesSeedingContinuation,
        prange, pidx, ics = _ics_from_grid(continuation);
        samples_per_parameter = 100, show_progress = true,
    )
    # show_progress && @info "Starting basins fraction continuation."
    # show_progress && @info "p = $(prange[1])"
    progress = ProgressMeter.Progress(length(prange);
        desc="Continuating basins fractions:", enabled=show_progress
    )
    n, spp = length(prange), samples_per_parameter
    (; mapper, method, threshold) = continuation

    
    ProgressMeter.next!(progress; showvalues = [("previous parameter", prange[1]),])

    labels, current_atts = get_attractors_and_labels(mapper, continuation, ics, pidx, prange[1], [], spp)
    att_v = deepcopy(current_atts)
    # Continue loop over all remaining parameters
    for p in prange[2:end]
        lbs, current_atts = get_attractors_and_labels(mapper, continuation, ics, pidx, p, prev_atts, spp)
        push!(att_v,current_atts)
        push!(labels,lbs)
        overwrite_dict!(prev_atts, current_atts)
        ProgressMeter.next!(progress; showvalues = [("previous parameter", p),])
    end
    # Normalize to smaller available integers for user convenience
    return fractions_curves, attractors_info
end
