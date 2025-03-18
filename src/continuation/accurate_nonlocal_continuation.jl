
# make sure to allow the possiblity that the proximity options can also be
# vectors of same length as `pcurve`; Same for the distributions
function find_my_name(ds::DynamicalSystem, attractors_cont, pcurve, distributions, εs, proximity_mapper_options...)

    results = [] #initialize properly
    for (i, p) in pcurve
        ε = εs <: AbstractVector ? εs[i] : εs # if its a vector, get i-th entry
        d = distributions <: AbstractVector ? distributions[i] : distributions # if its a vector, get i-th entry
        set_parameters!(ds, p)
        attractors = attractors_cont[i]
        mapper = NonlocalStabilityAccumulator(
            AttractorsViaProximity(ds, attractors, ε, proximity_mapper_options...)
        )
        fs = basins_fractions(mapper, ics)
        results[i] = finalizue_accumulator(mapper)
    end

    # transpoe results like in the new global continuation
    measures_cont = transpose
    return measures_cont
end
