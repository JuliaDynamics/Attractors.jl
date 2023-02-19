export RecurrencesContinuation
import ProgressMeter
using Random: MersenneTwister

# The recurrences based method is rather flexible because it works
# in two independent steps: it first finds attractors and then matches them.
struct RecurrencesContinuation{A, M, S, E} <: BasinsFractionContinuation
    mapper::A
    method::M
    threshold::Float64
    seeds_from_attractor::S
    info_extraction::E
end

"""
    RecurrencesContinuation(mapper::AttractorsViaRecurrences; kwargs...)
A method for [`basins_fractions_continuation`](@ref). It performs a continuation using a
mapper from [`AttractorsViaRecurrences`](@ref). The method uses two different 
continuation method: 

- Attractors are grouped across a parameor a single parameter slice. The continuation
is performed by matching the labelled attractors from one parameter slice to the next 
using the special distance `method`. The threshold   
- The attractors are grouped over the full parameter range. Said differently, all 
attractors of all initial conditions across all parameter values are put into the same 
"pool" and then grouped using the DBSCAN algorithm. After the grouping is finished the 
feature label fractions are distributedi to each parameter value they came from.
The optimal radius parameter for the DBSCAN algorithm is set though the threshold 
keyword argument. 

It uses seeding of previous attractors to find new ones, which is the main performance
bottleneck. The method uses [`match_attractor_ids!`](@ref) to match attractors
as the system parameter is increased.

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
function RecurrencesContinuation(
        mapper::AttractorsViaRecurrences; method = Centroid(),
        threshold = Inf, seeds_from_attractor = _default_seeding_process,
        info_extraction = identity
    )
    return RecurrencesContinuation(
        mapper, method, threshold, seeds_from_attractor, info_extraction
    )
end

function _default_seeding_process(attractor::AbstractDataset; rng = MersenneTwister(1))
    max_possible_seeds = 10
    seeds = round(Int, log(10, length(attractor)))
    seeds = clamp(seeds, 1, max_possible_seeds)
    return (rand(rng, attractor.data) for _ in 1:seeds)
end

"""
    basins_fractions_continuation(continuation::RecurrencesContinuation, prange, pidx, ics; kwargs...) → fractions_curves::Dict, attractors_info::Dict

Performs the continuation using a `mapper` from [`AttractorsViaRecurrences`](@ref) that 
maps an initial condition to an attractor. The structure `continuation` is a instance of [`RecurrencesContinuation`](@ref) that contains the information of the dynamical system. 

`prange` is the range of parameters for the continuation. 
`pidx` is the number or name of the parameter in the dynamical system.    
`ics` is a grid of initial conditions or a dedicated function that generates initial 
conditions in the state space.  

RecurrencesContinuation Keyword Arguments
- `show_progress = true`: print information on the current computation. 
- `samples_per_parameter = 100`: number of initial conditions to process per parameter. 
- `group_method = :grouping`: selects the method to perform the continuation. `:grouping` 
is meant to group the features accross the parameter range while `:matching` will match the clusters of attractors from one parameter slice to the next. 
"""
function basins_fractions_continuation(
        continuation::RecurrencesContinuation,
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

    # Gather labels, fractions and attractors doing the seeding process for each parameter.
    sav_labs = (group_method == :grouping)
    sav_labs && (labels = Vector{Int}(undef, n*spp))
    fractions_curves = Vector{Dict{Int, Float64}}(undef, n)
    attractors_info = Vector{Dict}(undef, n)
    prev_atts = Dict()
    for (i,p) in enumerate(prange)
        current_atts, fs, lab = get_attractors_and_fractions(mapper, continuation,
            ics, pidx, p, prev_atts, spp)
        fractions_curves[i] = fs
        attractors_info[i] =  get_info(current_atts)
        sav_labs && (labels[((i - 1)*spp + 1):i*spp] .= lab)
        overwrite_dict!(prev_atts, current_atts)
        ProgressMeter.next!(progress; showvalues = [("previous parameter", p),])
    end

    if group_method == :matching
        # Do the matching from one parameter to the next.
        match_attractors_forward!(attractors_info, fractions_curves, method, threshold)
    elseif group_method == :grouping
        # Group over the all range of parameters
        group_attractors!(attractors_info, labels,
            fractions_curves, n, spp, method, threshold)
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
    fs, lab, att = basins_fractions(mapper, ics; show_progress = false, N = N, additional_fs = labels,
        return_all_info = true)
    return att, fs, lab
end

# This function matches the attractors from one parameter slice to the next. 
# increasing the value of the parameter. Other matching are possible. For example 
# backward. 
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
end

# This function groups the attractors using the DBSCAN algorithm. The optimal 
# radius for the algorithm is set by threshold.
function group_attractors!(attractors, labels, fractions_curves, n,
    spp, method,  threshold)
    att = merge(attractors...)
    # Do the clustering with custom threshold
    att_keys, grouped_labels = clustering(att, method, threshold)

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
function clustering(att::Dict, method, threshold)
    distances = datasets_sets_distances(att, att, method)
    dist_mat = [distances[i][j] for i in keys(distances), j in keys(distances)]
    att_keys = collect(keys(distances))
    grouped_labels = _cluster_distances_into_labels(dist_mat, threshold, 1)
    return att_keys, grouped_labels
end
