using Test, Attractors
using Attractors.DynamicalSystemsBase
using Random

@testset "Henon map" begin
    ds = Systems.henon(; b = -0.9, a = 1.4)
    prange = range(0.6, 1.1; length = 20)
    pidx = 1; spp = 10000
    xg = yg = range(-2,2, length = 1000)
    grid = (xg,yg)
    mapper = Attractors.AttractorsViaRecurrences(ds, grid; sparse = true)
    sampler, = statespace_sampler(Random.MersenneTwister(1234); min_bounds = minimum.(grid), max_bounds = maximum.(grid))

    ## CLUSTERING CONTINUATION 
    dts_dis(x,y)= dataset_distance(x,y, Hausdorff())
    gc = GroupViaClustering(;clust_distance_metric = dts_dis, optimal_radius_method = 1., min_neighbors = 1)
    continuation = ClusteringAttractorsContinuation(mapper, par_weight = 1., group_config = gc)
    fs, att = basins_fractions_continuation(
        continuation,  prange, pidx, sampler; 
        show_progress = false, samples_per_parameter = spp)
    
    ## RECURENCE CONTINUATION
    sampler, = statespace_sampler(Random.MersenneTwister(1234); min_bounds = minimum.(grid), max_bounds = maximum.(grid))
    continuation = RecurrencesSeedingContinuation(mapper;
        threshold = 0.99, method = dts_dis
    )
    fs2, att2 = basins_fractions_continuation(
        continuation, prange, pidx, sampler;
        show_progress = false, samples_per_parameter = spp
    )
    # Compare fractions 
    for k in 1:length(prange)
        v1 = collect(values(fs[k]))
        v2 = collect(values(fs2[k]))
        for d in  ((sort(v1) .- sort(v2)) .< 0.001)
            @test d
        end
    end

    # compare number of labels detected
    k1 = Int[]; k2 = Int[];
    for k in 1:length(prange)
        push!(k1, collect(keys(fs[k]))...)
        push!(k2, collect(keys(fs2[k]))...)
    end
    @test sort(unique(k1)) == sort(unique(k2))
end

