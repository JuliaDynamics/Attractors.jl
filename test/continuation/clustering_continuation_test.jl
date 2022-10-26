using Test, Attractors
using Attractors.DynamicalSystemsBase, Attractors.DelayEmbeddings
using ChaosTools
using Random

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

    sampler, = ChaosTools.statespace_sampler(Random.MersenneTwister(1234);
        min_bounds = [-3., -3.], max_bounds = [3., 3.])

    rrange = range(0., 2; length = 20)
    ridx = 1

    # the features need some noise, if they are too well defined then the attractors are not 
    # clustered together.
    featurizer(a, t) = a[end,:] .+ rand(2)*0.1
    clusterspecs = Attractors.ClusteringConfig()
    clusterspecs.optimal_radius_method = "silhouettes"
    mapper = Attractors.AttractorsViaFeaturizing(ds, featurizer, clusterspecs; T = 20)
    continuation = ClusteringAcrossParametersContinuation(mapper)
    fractions_curves, attractors_info = Attractors.basins_fractions_continuation(
    continuation, rrange, ridx, sampler;
    show_progress = true, samples_per_parameter = 100, par_weight = 0.1)


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

