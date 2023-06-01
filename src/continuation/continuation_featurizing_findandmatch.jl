export FeaturizingFindAndMatch
import ProgressMeter
using Random: MersenneTwister

struct FeaturizingFindAndMatch{A, M, R<:Real, S, E} <: AttractorsBasinsContinuation
    mapper::A
    distance::M
    threshold::R
    seeds_from_attractor::S
    info_extraction::E
end

"""
Very similar to the recurrences version, only difference being that the seeding from attractors is different. 
"""
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

"""
Continuation here is very similar to the one done with recurrences. The difference is only
in how the ics from previous attractors are seeded to new parameters. In this case, we get ics sampled the previous attractors and pass them to 
basins_fractions, which extracts features from them and pushes them together with the other features. 
This could be generalized somehow so that one function could deal with both of the mappers, reducing this code duplication.
"""
function continuation(
        fam::FeaturizingFindAndMatch,
        prange, pidx, ics;
        samples_per_parameter = 100, show_progress = true,
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
    fs, _ = basins_fractions(mapper, ics; show_progress = false, N = samples_per_parameter)
    # At each parmaeter `p`, a dictionary mapping attractor ID to fraction is created.
    fractions_curves = [fs]
    # Furthermore some info about the attractors is stored and returned
    prev_attractors = deepcopy(mapper.attractors)
    get_info = attractors -> Dict(
        k => fam.info_extraction(att) for (k, att) in attractors
    )
    info = get_info(prev_attractors)
    attractors_info = [info]
    ProgressMeter.next!(progress; showvalues = [("previous parameter", prange[1]),])
    # Continue loop over all remaining parameters
    for p in prange[2:end]
        set_parameter!(mapper.ds, pidx, p)
        reset!(mapper)
        
        # Collect ics from previous attractors to pass as additional ics to basins fractions (seeding)
        # to ensure that the clustering will identify them as clusters, we need to guarantee that there
        # are at least `min_neighbors` entries
        additional_ics = Dataset(vcat(map(att-> 
            fam.seeds_from_attractor(att, fam.mapper.group_config.min_neighbors),
            values(prev_attractors))...)) #dataset with ics seeded from previous attractors
        
        # Now perform basin fractions estimation as normal, utilizing found attractors
        fs, _ = basins_fractions(mapper, ics;
            show_progress = false, N = samples_per_parameter, additional_ics
        )
        current_attractors = mapper.attractors
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

function reset!(mapper::AttractorsViaFeaturizing)
    empty!(mapper.attractors)
end