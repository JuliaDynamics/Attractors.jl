using Test, Attractors
using Attractors.DynamicalSystemsBase, Attractors.DelayEmbeddings
using Random


@testset "Henon map" begin

    ds = Systems.henon(; b = -0.9, a = 1.4)
    psorig = range(0.6, 1.1; length = 10)
    pidx = 1
    sampler, = statespace_sampler(Random.MersenneTwister(1234);
        min_bounds = [-2,-2], max_bounds = [2,2]
    )

    ps = psorig
    # Feature based on period. 
    function featurizer(a, t) 
        if abs(a[end,1]) > 100 || isnan(a[end,1])
            return [100]
        end
        tol = 1e-5
        if abs(a[end-1,1] - a[end,1]) < tol
            # period 1
            return [1] 
        elseif abs(a[end-3,1] - a[end,1]) < tol
            # period 3
            return [3]
        else
            return [100]
        end
    end
    clusterspecs = Attractors.ClusteringConfig()
    mapper = Attractors.AttractorsViaFeaturizing(ds, featurizer, clusterspecs; T = 500)
    continuation = ClusteringAcrossParametersContinuation(mapper)
    fractions_curves, attractors_info = Attractors.basins_fractions_continuation(
    continuation, psorig, pidx, sampler;
    show_progress = true, samples_per_parameter = 1000, par_weight = 1., ϵ_optimal = 1.)

    for (i, p) in enumerate(psorig)

        fs = fractions_curves[i]
        if   p < 0.9 
            k = sort!(collect(keys(fs)))
            @test length(k) == 2
        elseif p > 1. 
            k = sort!(collect(keys(fs)))
            @test length(k) == 3
        end    
        @test sum(values(fs)) ≈ 1
    end


    # Test mmap 
    fractions_curves, attractors_info = Attractors.basins_fractions_continuation(
    continuation, psorig, pidx, sampler;
    show_progress = true, samples_per_parameter = 1000, par_weight = 1., ϵ_optimal = 1., mmap_limit = 1000)

    for (i, p) in enumerate(psorig)
        fs = fractions_curves[i]
        if   p < 0.9 
            k = sort!(collect(keys(fs)))
            @test length(k) == 2
        elseif p > 1. 
            k = sort!(collect(keys(fs)))
            @test length(k) == 3
        end    
        @test sum(values(fs)) ≈ 1
    end

end


@testset "Dummy bistable map" begin

    function dumb_map(dz, z, p, n)
        x,y = z
        r = p[1]

        if r < 0.5
            dz[1] = dz[2] = 0.
        else
            if x > 0 
                dz[1] = r
                dz[2] = r
            else
                dz[1] = -r
                dz[2] = -r
            end
        end
        return
    end

    r = 3.833
    ds = DiscreteDynamicalSystem(dumb_map, [0., 0.], [r])

    sampler, = statespace_sampler(Random.MersenneTwister(1234);
        min_bounds = [-3., -3.], max_bounds = [3., 3.])

    rrange = range(0., 2; length = 20)
    ridx = 1

    featurizer(a, t) = a[end,:] 
    clusterspecs = Attractors.ClusteringConfig()
    # clusterspecs.optimal_radius_method = "silhouettes"
    mapper = Attractors.AttractorsViaFeaturizing(ds, featurizer, clusterspecs; T = 20)
    continuation = ClusteringAcrossParametersContinuation(mapper)
    fractions_curves, attractors_info = Attractors.basins_fractions_continuation(
    continuation, rrange, ridx, sampler;
    show_progress = true, samples_per_parameter = 100, par_weight = 1., ϵ_optimal = 0.1)


    for (i, r) in enumerate(rrange)

        fs = fractions_curves[i]
        if r < 0.5
            k = sort!(collect(keys(fs)))
            @test length(k) == 1
        else 
            k = sort!(collect(keys(fs)))
            @test length(k) == 2
            v = values(fs) 
            for f in v
                @test (0.4 < f < 0.6)
            end
        end    
        @test sum(values(fs)) ≈ 1
    end
end
