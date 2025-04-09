function global_continuation(
    acam::AttractorSeedContinueMatch{<:StabilityMeasuresAccumulator}, pcurve, ics;
        samples_per_parameter = 100, show_progress = true,
    )
    N = samples_per_parameter
    progress = ProgressMeter.Progress(
        length(pcurve); desc = "Continuing attractors and stability:", enabled=show_progress
    )
    mapper = acam.mapper
    reset_mapper!(mapper)

    set_parameters!(referenced_dynamical_system(mapper), pcurve[1])
    if ics isa Function
        fs = basins_fractions(mapper, ics; show_progress = false, N = samples_per_parameter)
    else # we ignore labels in this continuation algorithm
        fs, = basins_fractions(mapper, ics; show_progress = false)
    end

    # this is the first difference with the standard continuation: we have to call finalize.
    measures = finalize_accumulator(mapper)
    measures_cont = [measures]

    prev_attractors = deepcopy(extract_attractors(mapper))
    attractors_cont = [deepcopy(prev_attractors)] # we need the copy
    ProgressMeter.next!(progress)

    # Continue loop over all remaining parameters
    for p in @view(pcurve[2:end])
        set_parameters!(referenced_dynamical_system(mapper), p)
        reset_mapper!(mapper)
        if typeof(mapper.mapper) <: AttractorsViaRecurrences
            fs = if allows_mapper_u0(mapper)
                seed_attractors_to_fractions_individual(mapper,
                    prev_attractors, ics, N, acam.seeding
                )
            else
                seed_attractors_to_fractions_grouped(mapper,
                    prev_attractors, ics, N, acam.seeding
                )
            end
            current_attractors = deepcopy(extract_attractors(mapper))
        else
            error("Unsupported mapper type: $(typeof(mapper.mapper))")
        end
        push!(attractors_cont, current_attractors)
        overwrite_dict!(prev_attractors, current_attractors)
        ProgressMeter.next!(progress)

        # this is the second difference in global continuation
        measures = finalize_accumulator(mapper)
        push!(measures_cont, measures)
    end

    rmaps = match_sequentially!(
        attractors_cont, acam.matcher; pcurve, ds = referenced_dynamical_system(mapper)
    )

    # This is the third difference in global continuation
    for i in 2:length(pcurve)
        for dict in values(measures_cont[i])
            swap_dict_keys!(dict, rmaps[i-1])
        end
    end
    transposed = Dict{String, Vector{Dict{Int64, Float64}}}()
    for measure in measures_cont[1]
        measure_name = measure[1]
        transposed[measure_name] = Vector{Dict{Int64, Float64}}()
    end

    for measures in measures_cont
        for (measure_name, measure_dict) in measures
            push!(transposed[measure_name], measure_dict)
        end
    end

    return transposed, attractors_cont
end



# make sure to allow the possiblity that the proximity options can also be
# vectors of same length as `pcurve`; Same for the distributions
function stability_measures_along_continuation(
    ds::DynamicalSystem, 
    attractors_cont, 
    pcurve,
    ics; 
    ε = nothing, 
    weighting_distribution = EverywhereUniform(), 
    finite_time = 1.0,
    samples_per_parameter = 1000, 
    metric=Euclidean(), 
    proximity_mapper_options=NamedTuple(), 
    show_progress=true
)
    N = samples_per_parameter
    progress = ProgressMeter.Progress(
        length(pcurve); desc = "Continuing attractors and stability:", enabled=show_progress
    )
    measures_cont = []
    for (i, p) in enumerate(pcurve)
        ε_ = ε isa AbstractVector ? ε[i] : ε # if its a vector, get i-th entry
        weighting_distribution = weighting_distribution isa AbstractVector ?
                                weighting_distribution[i] : weighting_distribution
        finite_time = finite_time isa AbstractVector ? finite_time[i] : finite_time
        set_parameters!(ds, p)
        attractors = attractors_cont[i]
        accumulator = StabilityMeasuresAccumulator(
            AttractorsViaProximity(ds, attractors; ε = ε_, proximity_mapper_options...);
            weighting_distribution=weighting_distribution, finite_time=finite_time,
            metric=metric
        )
        N = ics isa Function ? N : length(ics)
        for i ∈ 1:N
            ic = _get_ic(ics, i)
            id = accumulator(ic) # accumulate stability measures for given i.c.
        end
        push!(measures_cont, finalize_accumulator(accumulator))
        ProgressMeter.next!(progress)
    end

    transposed = Dict{String, Vector{Dict{Int64, Float64}}}()
    for measure in measures_cont[1]
        measure_name = measure[1]
        transposed[measure_name] = Vector{Dict{Int64, Float64}}()
    end

    for measures in measures_cont
        for (measure_name, measure_dict) in measures
            push!(transposed[measure_name], measure_dict)
        end
    end
    return transposed
end
