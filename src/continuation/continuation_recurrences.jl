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
        samples_per_parameter = 100, show_progress = true, group_method = :matching
    )
    progress = ProgressMeter.Progress(length(prange);
        desc="Continuating basins fractions:", enabled=show_progress
    )
    n, spp = length(prange), samples_per_parameter
    (; mapper, method, threshold) = continuation
    get_info = attractors -> Dict(
        k => continuation.info_extraction(att) for (k, att) in attractors
    )

    sav_labs = (group_method == :grouping) 
    sav_labs && (labels = Vector{Int}(undef, n*spp))
    # first parameter is run in isolation, as it has no prior to seed from
    current_atts, fs, lab = get_attractors_and_fractions(mapper, continuation, ics, pidx, prange[1], [], spp)
    prev_atts = deepcopy(current_atts)
    fractions_curves = [fs]
    attractors_info = [get_info(prev_atts)]
    ProgressMeter.next!(progress; showvalues = [("previous parameter", prange[1]),])
    sav_labs && (labels[1:spp] .= lab)
    # Continue loop over all remaining parameters
    for (i,p) in enumerate(prange[2:end])
        current_atts, fs, lab = get_attractors_and_fractions(mapper, continuation, ics, pidx, p, prev_atts, spp)
        push!(fractions_curves, fs)
        push!(attractors_info, get_info(current_atts))
        sav_labs && (labels[((i - 1)*spp + 1):i*spp] .= lab)
        overwrite_dict!(prev_atts, current_atts)
        ProgressMeter.next!(progress; showvalues = [("previous parameter", p),])
    end

    # Do the matching from one parameter to the next.
    match_attractors_forward!(attractors_info, fractions_curves, method, threshold)

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

function get_attractors_and_fractions(mapper, continuation, ics, pidx, p, prev_atts, N)
    set_parameter!(mapper.integ, pidx, p)
    reset!(mapper)
    labels = seed_attractors(prev_atts, mapper, continuation)
    fs, lab, att = basins_fractions(mapper, ics; show_progress = false, N = N, additional_fs = labels,
        return_all_info = true)
    return att, fs, lab
end


function match_attractors_forward!(attractors, fractions, method, threshold)
    n = length(attractors)
    for k in 2:n
        if !isempty(attractors[k]) && !isempty(attractors[k-1])
            # If there are any attractors,
            # match with previous attractors before storing anything!
            rmap = match_attractor_ids!(attractors[k], attractors[k-1]; method, threshold)
            swap_dict_keys!(fractions[k], rmap)
        end
    end
    # Normalize to smaller available integers for user convenience
    rmap = retract_keys_to_consecutive(fractions)
    for (da, df) in zip(attractors, fractions)
        swap_dict_keys!(da, rmap)
        swap_dict_keys!(df, rmap)
    end
end


function group_attractors(attractors, label, threshold)

    # Compute distances. 
    distances = datasets_sets_distances(att_v, att_v) 
    dist_mat = zeros(length(distances),length(distances))
    att_keys = Vector{Int}()
    for i in keys(distances)
        push!(att_keys, i)
        for j in keys(distances[i])
            dist_mat[i,j] =  distances[i][j]
        end
    end
    group_labels = _cluster_distances_into_labels(dist_mat, threshold, 1)

    # Now rename the labels, get the fractions and pack attractors.
    grouped_labels, keys = group_attractors(att_v, 0.001)
    grouped_labels .+= maximum(keys)
    @show grouped_labels
    @show keys
    for (k,group_lab)  in enumerate(grouped_labels)
        if group_lab != -1
            att_num = keys[k]
            labels[labels .== att_num] .= grouped_labels[k]
        end
    end

    fractions_curves = get_fractions(labels, n, spp)
    return fractions_curves
end

function get_fractions(labels, n, spp)
    # finally we collect/group stuff into their dictionaries
    fractions_curves = Vector{Dict{Int, Float64}}(undef, n)
   for i in 1:n
        current_labels = view(labels, ((i - 1)*spp + 1):i*spp)
        current_ids = unique(current_labels)
        # getting fractions is easy; use API function that takes in arrays
        fractions_curves[i] = basins_fractions(current_labels)
    end
    return fractions_curves
end


function label_fractions_across_parameter(labels, n, spp, feature, info_extraction)
    # finally we collect/group stuff into their dictionaries
    fractions_curves = Vector{Dict{Int, Float64}}(undef, n)
    dummy_info = info_extraction(feature)
    attractors_info = Vector{Dict{Int, typeof(dummy_info)}}(undef, n)
    for i in 1:n
        # Here we know which indices correspond to which parameter value
        # because they are sequentially increased every `spp`
        # (steps per parameter)
        current_labels = view(labels, ((i - 1)*spp + 1):i*spp)
        current_ids = unique(current_labels)
        # getting fractions is easy; use API function that takes in arrays
        fractions_curves[i] = basins_fractions(current_labels)
        attractors_info[i] = Dict(
            id => info_extraction(
                view(current_labels, findall(isequal(id), current_labels))
            ) for id in current_ids
        )
    end
    return fractions_curves, attractors_info
end
