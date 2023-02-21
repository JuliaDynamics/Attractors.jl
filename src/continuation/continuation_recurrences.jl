export RecurrencesContinuation
import ProgressMeter
using Random: MersenneTwister

# The recurrences based method is rather flexible because it works
# in two independent steps: it first finds attractors and then matches them.
struct RecurrencesContinuation{A, D, S, E, M} <: AttractorsBasinsContinuation
    mapper::A
    distance::D
    threshold::Float64
    seeds_from_attractor::S
    info_extraction::E
    matching_method::M
end

"""
    RecurrencesContinuation <: AttractorsBasinsContinuation
    RecurrencesContinuation(mapper::AttractorsViaRecurrences; kwargs...)

A method for [`continuation`](@ref). TODO: Cite our preprint here.

## Description

At the first parameter slice attractors are found as described in the
[`AttractorsViaRecurrences`](@ref) mapper using recurrences in state space.
At each subsequent parameter slice,
new attractors are found by seeding initial conditions from the previously found
attractors and then piping these initial conditions through the recurrences algorithm
of the `mapper`. Seeding initial conditions close to previous attractors accelerates
the main bottleneck of [`AttractorsViaRecurrences`](@ref), which is finding the attractors.
This process continues until all parameter values are exhausted and for each parameter
value the attractors and their fractions are found.

Then, the different attractors across parameters are matched so that they have
the same ID. The matching process is based on distances attractors (= sets in state space)
have between each other. The function that computes these distances is
[`setsofsets_distances`](@ref) (please read that docstring before continuing).

## Keyword arguments

- `matching_method = ParameterSliceCrossDistance()`: see description above.
- `distance = Centroid()`: the distance used in [`set_distance`](@ref) or
  [`setsofsets_distances`](@ref)` to provide a distance for matching attractors.
- `threshold = Inf`: given to [`match_attractor_ids!`](@ref) if
  `matching_method isa ParameterSliceCrossDistance`.
- `info_extraction = identity`: A function that takes as an input an attractor (`StateSpaceSet`)
  and outputs whatever information should be stored. It is used to return the
  `attractors_info` in [`continuation`](@ref).
- `seeds_from_attractor`: A function that takes as an input an attractor and returns
  an iterator of initial conditions to be seeded from the attractor for the next
  parameter slice. By default, we sample some points from existing attractors according
  to how many points the attractors themselves contain. A maximum of `10` seeds is done
  per attractor.
"""
function RecurrencesContinuation(
        mapper::AttractorsViaRecurrences; distance = Centroid(),
        threshold = Inf, seeds_from_attractor = _default_seeding_process,
        info_extraction = identity, matching_method = ParameterSliceCrossDistance(),
    )
    return RecurrencesContinuation(
        mapper, distance, threshold, seeds_from_attractor, info_extraction,
        matching_method,
    )
end

function _default_seeding_process(attractor::AbstractStateSpaceSet; rng = MersenneTwister(1))
    max_possible_seeds = 10
    seeds = round(Int, log(10, length(attractor)))
    seeds = clamp(seeds, 1, max_possible_seeds)
    return (rand(rng, attractor.data) for _ in 1:seeds)
end




function continuation(
        rc::RecurrencesContinuation,
        prange, pidx, ics = _ics_from_grid(rc);
        samples_per_parameter = 100, show_progress = true, matching_method = :matching
    )
    progress = ProgressMeter.Progress(length(prange);
        desc="Continuating basins fractions:", enabled=show_progress
    )
    n, spp = length(prange), samples_per_parameter
    (; mapper, distance, threshold) = rc
    get_info = attractors -> Dict(
        k => rc.info_extraction(att) for (k, att) in attractors
    )
    matching_method = rc.matching_method

    # Gather labels, fractions and attractors doing the seeding process for each parameter.
    sav_labs = (matching_method isa ClusterDistanceMatrix)
    sav_labs && (labels = Vector{Int}(undef, n*spp))
    fractions_curves = Vector{Dict{Int, Float64}}(undef, n)
    attractors_info = Vector{Dict}(undef, n)
    prev_atts = Dict()
    for (i,p) in enumerate(prange)
        current_atts, fs, lab = get_attractors_and_fractions(mapper, rc,
            ics, pidx, p, prev_atts, spp)
        fractions_curves[i] = fs
        attractors_info[i] =  get_info(current_atts)
        sav_labs && (labels[((i - 1)*spp + 1):i*spp] .= lab)
        overwrite_dict!(prev_atts, current_atts)
        ProgressMeter.next!(progress; showvalues = [("previous parameter", p),])
    end

    if matching_method isa ParameterSliceCrossDistance
        # Do the matching from one parameter to the next.
        match_attractors_forward!(attractors_info, fractions_curves, distance, threshold)
    elseif matching_method isa ClusterDistanceMatrix
        # Group over the all range of parameters
        group_attractors!(attractors_info, labels,
            fractions_curves, n, spp, distance, threshold)
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

function _ics_from_grid(continuation::RecurrencesContinuation)
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

function get_attractors_and_fractions(mapper, continuation, ics, pidx, p, prev_atts, N)
    set_parameter!(mapper.integ, pidx, p)
    reset!(mapper)
    labels = seed_attractors(prev_atts, mapper, continuation)
    # TODO: change basins fractions to ALWAYS return fs, lab, att
    fs, lab, att = basins_fractions(mapper, ics; show_progress = false, N = N, additional_fs = labels,
        return_all_info = true)
    return att, fs, lab
end

# This function matches the attractors from one parameter slice to the next.
# increasing the value of the parameter. Other matching are possible. For example
# backward.
function match_attractors_forward!(attractors, fractions, distance, threshold)
    n = length(attractors)
    for k in 2:n
        if !isempty(attractors[k]) && !isempty(attractors[k-1])
            # If there are any attractors,
            # match with previous attractors before storing anything!
            rmap = match_attractor_ids!(attractors[k], attractors[k-1]; distance = distance, threshold)
            swap_dict_keys!(fractions[k], rmap)
        end
    end
end

# This function groups the attractors using the DBSCAN algorithm. The optimal
# radius for the algorithm is set by threshold.
function group_attractors!(attractors, labels, fractions_curves, n,
    spp, distance,  threshold)
    att = merge(attractors...)
    # Do the clustering with custom threshold
    att_keys, grouped_labels = clustering(att, distance, threshold)

    # Now rename the labels (we ignore -1 in grouped labels),
    # get the fractions and pack attractors.
    postve_lab = findall(grouped_labels .> 0)
    grouped_labels[postve_lab] .+= maximum(att_keys)
    rmap = Dict( att_keys[k] => grouped_labels[k] for k in postve_lab)
    replace!(labels, rmap...)

    # finally we collect/group stuff into their dictionaries
    for i in 1:n
        swap_dict_keys!(attractors[i], rmap)
        current_labels = view(labels, ((i - 1)*spp + 1):i*spp)
        fractions_curves[i] = basins_fractions(current_labels)
    end

end

# DBSCAN algorithm performed using the low-level function
# _cluster_distances_into_labels
function clustering(att::Dict, distance, threshold)
    distances = datasets_sets_distances(att, att, distance)
    dist_mat = [distances[i][j] for i in keys(distances), j in keys(distances)]
    att_keys = collect(keys(distances))
    grouped_labels = _cluster_distances_into_labels(dist_mat, threshold, 1)
    return att_keys, grouped_labels
end
