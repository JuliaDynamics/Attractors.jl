using Test, Attractors
using Attractors.DynamicalSystemsBase
using Random

# @testset "Clustering Attractors: Henon map" begin
    ds = Systems.henon(; b = -0.9, a = 1.4)
    prange = range(0.6, 1.1; length = 10)
    pidx = 1; spp = 1000
    xg = yg = range(-2,2, length = 1000)
    grid = (xg,yg)
    mapper = Attractors.AttractorsViaRecurrences(ds, grid; sparse = true)
    sampler, = statespace_sampler(Random.MersenneTwister(1234); min_bounds = minimum.(grid), max_bounds = maximum.(grid))

    ## CLUSTERING CONTINUATION 
    dts_dis(x,y)= dataset_distance(x,y, Hausdorff())
    gc = GroupViaClustering(;clust_distance_metric = dts_dis, optimal_radius_method = 1., min_neighbors = 1)
    continuation = ClusteringAttractorsContinuation(mapper; par_weight = 1., group_config = gc)
    fs, att, fsj = basins_fractions_continuation(
        continuation,  prange, pidx, sampler; 
        show_progress = true, samples_per_parameter = spp)
    
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
        @show v1, v2
        # @test length(v1) == length(v2)
        # for d in  ((sort(v1) .- sort(v2)) .< 0.001)
        #     @test d
        # end
    end

    # compare number of labels detected
    k1 = Int[]; k2 = Int[];
    for k in 1:length(prange)
        push!(k1, collect(keys(fs[k]))...)
        push!(k2, collect(keys(fs2[k]))...)
    end
    @test sort(unique(k1)) == sort(unique(k2))
end

# @testset "Clustering Attractors: magnetic pendulum" begin
    d, α, ω = 0.3, 0.2, 0.5
    ds = Systems.magnetic_pendulum(; d, α, ω)
    xg = yg = range(-3, 3, length = 101); grid = (xg, yg)
    ds = projected_integrator(ds, 1:2, [0.0, 0.0])
    mapper = AttractorsViaRecurrences(ds, (xg, yg); Δt = 1.0)
    rr = range(1, 0; length = 21)
    ps = [[1, 1, γ] for γ in rr]
    pidx = :γs
    spp = 100

    # RECURENCE CONTINUATION
    sampler, = statespace_sampler(Random.MersenneTwister(1234); min_bounds = minimum.(grid), max_bounds = maximum.(grid))
    continuation = RecurrencesSeedingContinuation(mapper; threshold = Inf)
    fs2, att2 = basins_fractions_continuation(
        continuation, ps, pidx, sampler; show_progress = true, samples_per_parameter = spp
    )
    ## CLUSTERING CONTINUATION 
    sampler, = statespace_sampler(Random.MersenneTwister(1234); min_bounds = minimum.(grid), max_bounds = maximum.(grid))
    dts_dis(x,y)= dataset_distance(x,y, Hausdorff())
    gc = GroupViaClustering(;clust_distance_metric = dts_dis, optimal_radius_method = .2, min_neighbors = 1)
    continuation = ClusteringAttractorsContinuation(mapper; par_weight = 1., group_config = gc)
    fs, att = basins_fractions_continuation(
        continuation,  ps, pidx, sampler; 
        show_progress = true, samples_per_parameter = spp)
    
     
    # Compare fractions 
    for k in 1:length(ps)
        v1 = collect(values(fs[k]))
        v2 = collect(values(fs2[k]))
        @show v1, v2
        for d in  ((sort(v1) .- sort(v2)) .< 0.001)
            @test d
        end
    end

    # compare number of labels detected
    k1 = Int[]; k2 = Int[];
    for k in 1:length(ps)
        push!(k1, collect(keys(fs[k]))...)
        push!(k2, collect(keys(fs2[k]))...)
    end
    @test sort(unique(k1)) == sort(unique(k2))
# end
