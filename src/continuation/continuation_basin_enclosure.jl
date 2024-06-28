struct BasinEnclosure <: AttractorsBasinsContinuation
    mapper
    coflowing_threshold
end

function continuation(basen::BasinEnclosure, prange, pidx, ics;
        samples_per_parameter = 100, show_progress = true,
    )
    progress = ProgressMeter.Progress(length(prange);
        desc = "Continuing attractors and basins:", enabled=show_progress
    )
    mapper = basen.mapper
    reset_mapper!(mapper)
    # first parameter is run in isolation, as it has no prior to seed from
    set_parameter!(referenced_dynamical_system(mapper), pidx, prange[1])
    fs = basins_fractions(mapper, ics; show_progress = false, N = samples_per_parameter)
    # At each parmaeter `p`, a dictionary mapping attractor ID to fraction is created.
    fractions_cont = [fs]
    # The attractors are also stored (and are the primary output)
    prev_attractors = deepcopy(extract_attractors(mapper))
    attractors_cont = [prev_attractors]
    ProgressMeter.next!(progress; showvalues = [("previous parameter", prange[1]),])
    # Continue loop over all remaining parameters
    for p in prange[2:end]
        # reset mapper
        set_parameter!(referenced_dynamical_system(mapper), pidx, p)
        reset_mapper!(mapper)
        # and obtain the previous attractors
        IDs₋ = sort(keys(prev_attractors))
        initconds_from_prev_atts = [prev_attractors[k][end] for k in IDs₋]

        # Optionally here we can do seeding from previously found attractors,
        # but let's skip it for now!
        # <seeding...>

        # The rest of the code is different depending on whether
        # the mapper supports the `mapper(u) -> id` functionality or not.
        # This we take care via multiple dispatch.
        # Here we obtain current fractions and attractors,
        # as well as the "enclosed map" which maps OLD IDs to NEWLY FOUND IDs
        fs, current_attractors, enclosed_map = basin_enclosure_next_step(
            mapper::AttractorsViaRecurrences, ics, N, IDs₋, initconds_from_prev_atts
        )
        # this enclosed map has to be cleaned, e.g., if we have multiple
        # attractors from the previous parameter converging to the same in current parameter.
        # In the end we get the replacement map, mapping CURRENT IDs,
        # to IDs in the previous parameter. I.e., "matching".
        rmap = resolve_coflowing(enclosed_map, args...) # @Kalel you have to do this
        # The rest of the code utilizes existing infrastructure based
        # on the replacement map

        # what we are left with is the replacment map. It maps keys of CURRENT
        # attractors to keys of the OLD attractors! So we match with the existing function:
        swap_dict_keys!(current_attractors, rmap)
        swap_dict_keys!(fs, rmap)
        # and store the result
        push!(fractions_cont, fs)
        push!(attractors_cont, current_attractors)
        # and replace old with new:
        overwrite_dict!(prev_attractors, current_attractors)
        ProgressMeter.next!(progress; showvalues = [("previous parameter", p),])
    end
    return fractions_cont, attractors_cont
end

# for mappers that support `mapper(u) -> ID`
# TODO: Have more correct dispatch here. Featurizing also works with histogram or
# nearest feature.
function basin_enclosure_next_step(mapper::AttractorsViaRecurrences, ics, N, IDs₋, initconds_from_prev_atts)
    # first we get attractors and basins as normal
    fs = basins_fractions(mapper, ics;
        additional_fs = seeded_fs, show_progress = false, N
    )

    # We extract the attractors
    current_attractors = deepcopy(extract_attractors(mapper))
    # and also map the initial conditions of the previous attractors to the new ones
    enclosed_map = Dict(k => mapper(initconds_from_prev_atts[k]) for k in IDs₋)
    return fs, current_attractors, enclosed_map
end

# for mappers that do not support `mapper(u) -> ID`
function basin_enclosure_next_step(mapper::AttractorsViaFeaturizing, ics, N, IDs₋, initconds_from_prev_atts)
    # We augment the list of initial conditions to run through the mapper with the extra ones
    M = length(IDs₋)
    if ics isa Function
        ics_collected = [copy(ics()) for _ in 1:N]
    else
        ics = copy(ics)
    end
    append!(ics_collected, initconds_from_prev_atts)

    fs, labels = basins_fractions(mapper, ics_collected; show_progress = false)
    # from here we do two modifications; first we obtain the labels of the special extra initial conditions
    extra_labels = labels[N+1:N+M]
    # which we transform to the enclosure map
    enclosed_map = Dict(k => extra_labels[i] for (i, k) in enumerate(IDs₋))

    # we also need to filter/adjust the basin fractions to the pure ones without these extra i.c.
    labels = labels[1:N]
    fs = basins_fractions(labels) # utility function

    # lastly, we extract the attractors
    current_attractors = deepcopy(extract_attractors(mapper))
    return fs, current_attractors, enclosed_map
end