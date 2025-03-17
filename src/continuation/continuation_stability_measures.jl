# TODO: Can we make this function not duplicate the source code
# of the original `ASCM` implementation?

function global_continuation(
    acam::AttractorSeedContinueMatch{<:StabilityMeasuresAccumulator}, pcurve, ics; distributions = p->accumulator.d,
        samples_per_parameter = 100, show_progress = true,
    )
    N = samples_per_parameter
    progress = ProgressMeter.Progress(length(pcurve); desc = "Continuing attractors and stability:", enabled=show_progress)
    mapper = acam.mapper
    reset_mapper!(mapper)

    set_parameters!(referenced_dynamical_system(accumulator), pcurve[1])
    if ics isa Function
        fs = basins_fractions(accumulator, ics; show_progress = false, N = samples_per_parameter)
    else # we ignore labels in this continuation algorithm
        fs, = basins_fractions(accumulator, ics; show_progress = false)
    end

    # this is the first difference with the standard continuation: we have to call finalize.
    measures = finalize_accumulator(mapper)
    measures_cont = [measures]

    prev_attractors = deepcopy(extract_attractors(accumulator))
    attractors_cont = [deepcopy(prev_attractors)] # we need the copy
    ProgressMeter.next!(progress)

    # Continue loop over all remaining parameters
    ### For @Datseris: Is this what you had in mind when we emailed about 
    ### the global continuation algorithm with AttractorsViaProximity?
    ### Currently, line 44 has an error because we cant pass ics to the AttractorsViaRecurrences.
    for p in @view(pcurve[2:end])
        set_parameters!(referenced_dynamical_system(mapper), p)
        reset_mapper!(mapper)
        mapper.d = distributions(p)
        if typeof(mapper.mapper) <: AttractorsViaRecurrences
            fs = if allows_mapper_u0(mapper)
                seed_attractors_to_fractions_individual(mapper, prev_attractors, ics, N, seeding)
            else
                seed_attractors_to_fractions_grouped(mapper, prev_attractors, ics, N, seeding)
            end
            current_attractors = deepcopy(extract_attractors(mapper))
        elseif typeof(mapper.mapper) <: AttractorsViaProximity
            dummy_mapper = AttractorsViaRecurrences(referenced_dynamical_system(mapper), ics)
            fs = if allows_mapper_u0(dummy_mapper)
                seed_attractors_to_fractions_individual(dummy_mapper, prev_attractors, ics, N, seeding)
            else
                seed_attractors_to_fractions_grouped(dummy_mapper, prev_attractors, ics, N, seeding)
            end
            current_attractors = deepcopy(extract_attractors(dummy_mapper))
            mapper.mapper.attractors = current_attractors
            A = ics_from_grid(grid)
            for u0 in A
                id = mapper(u0)
            end
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

    rmaps = match_sequentially!(attractors_cont, acam.matcher; pcurve, ds = referenced_dynamical_system(mapper))

    # This is the second difference: we change the output to something more pleasant to handle
    transposed_measures_cont = Dict{String, Vector{Dict{Int64, Float64}}}()
    for measure in measures_cont[1]
        measure_name = measure[1]
        transposed_measures_cont[measure_name] = Vector{Dict{Int64, Float64}}()
    end
    for measures in measures_cont
        for (measure_name, measure_dict) in measures
            push!(transposed_measures_cont[measure_name], measure_dict)
        end
    end

    # And lastly, we also match sequentially all dictionaries here
    for measure_cont in values(measures_cont)
        match_sequentially!(measure_cont, rmaps)
    end

    return measures_cont, attractors_cont
end

# TODO: Put this in some other function
#=
# Calculate distance to bifurcation
distance_to_bifurcation = Array{Dict{Int64, Float64}}([])
char_ret_time_cont = measures_cont["characteristic_return_time"]
attractors_cont_keys = unique(vcat([collect(keys(attractors_cont[i])) for i in 1:length(attractors_cont)]...))
prange = [pdict[1][2] for pdict in pcurve]
for i in 1:length(prange)
    push!(distance_to_bifurcation, Dict{Int64, Float64}())
    for key in attractors_cont_keys
        ps_unstable = [prange[j] for j in 1:length(prange) if !haskey(char_ret_time_cont[j], key) || (char_ret_time_cont[j][key]==Inf64 || isnan(char_ret_time_cont[j][key]))]
        if haskey(char_ret_time_cont[i], key) && length(ps_unstable)>0 && (char_ret_time_cont[i][key]!=Inf64 && !isnan(char_ret_time_cont[i][key]))
            distance_to_bifurcation[end][key] = minimum([abs(ps_unstable[k] - prange[i]) for k in 1:length(ps_unstable)])
        else
            distance_to_bifurcation[end][key] = NaN
        end
    end
end
merge!(measures_cont, Dict("distance_to_bifurcation" => distance_to_bifurcation))
=#