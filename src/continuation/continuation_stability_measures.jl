# global continuation algorithm tracking attractors and stability measures
function global_continuation(accumulator::StabilityMeasuresAccumulator, matcher, pcurve, ics; samples_per_parameter = 100, show_progress = true, seeding = A->[])#_default_seeding ### AM I had to remove the ::IDMatcher after matcher for it to precompile but I dont know why
    N = samples_per_parameter
    progress = ProgressMeter.Progress(length(pcurve); desc = "Continuing attractors and basins:", enabled=show_progress)
    mapper = accumulator
    reset_mapper!(mapper)

    set_parameters!(referenced_dynamical_system(mapper), pcurve[1])
    if ics isa Function
        fs = basins_fractions(mapper, ics; show_progress = false, N = samples_per_parameter)
    else # we ignore labels in this continuation algorithm
        fs, = basins_fractions(mapper, ics; show_progress = false)
    end

    # This is the third difference in global continuation
    measures = finalize(accumulator)
    measures_cont = [measures]

    prev_attractors = deepcopy(extract_attractors(mapper))
    attractors_cont = [deepcopy(prev_attractors)] # we need the copy
    ProgressMeter.next!(progress)

    # Continue loop over all remaining parameters
    for p in @view(pcurve[2:end])
        set_parameters!(referenced_dynamical_system(mapper), p)
        reset_mapper!(mapper)
        fs = if allows_mapper_u0(mapper)
            seed_attractors_to_fractions_individual(mapper, prev_attractors, ics, N, seeding)
        else
            seed_attractors_to_fractions_grouped(mapper, prev_attractors, ics, N, seeding)
        end
        current_attractors = deepcopy(extract_attractors(mapper))
        push!(attractors_cont, current_attractors)
        overwrite_dict!(prev_attractors, current_attractors)
        ProgressMeter.next!(progress)

        # this is the second difference in global continuation
        measures = finalize(accumulator)
        push!(measures_cont, measures)
    end

    rmaps = match_sequentially!(attractors_cont, matcher; pcurve, ds = referenced_dynamical_system(mapper))

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
    measures_cont = transposed

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
    return measures_cont, attractors_cont
end
