DO_EXTENSIVE_TESTS = get(ENV, "ATTRACTORS_EXTENSIVE_TESTS", "false") == "true"

using Test, Attractors
using Random

# TODO: Add a Histogram test with the competition dynamics model

@testset "Dummy bistable map" begin

    function dumb_map(dz, z, p, n)
        x, y = z
        r = p[1]
        if r < 0.5
            dz[1] = dz[2] = 0.0
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
        min_bounds = [-3.0, -3.0], max_bounds = [3.0, 3.0])

    rrange = range(0, 2; length = 20)
    ridx = 1

    featurizer(a, t) = a[end]
    clusterspecs = Attractors.GroupViaClustering(optimal_radius_method = "silhouettes", max_used_features = 200)
    mapper = Attractors.AttractorsViaFeaturizing(ds, featurizer, clusterspecs; T = 20, threaded = true)
    continuation = GroupAcrossParameterContinuation(mapper; par_weight = 1.)
    fractions_curves, attractors_info = Attractors.continuation(
    continuation, rrange, ridx, sampler; show_progress = false)


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
                @test (0.3 < f < 0.7)
            end
        end
        @test sum(values(fs)) ≈ 1
    end
end

if DO_EXTENSIVE_TESTS

    @testset "Henon period doubling" begin

        # Notice special parameter values:
        ds = Systems.henon(; b = -0.9, a = 1.4)
        ps = range(0.6, 1.1; length = 10)
        pidx = 1
        sampler, = statespace_sampler(Random.MersenneTwister(1234);
            min_bounds = [-2,-2], max_bounds = [2,2]
        )

        # Feature based on period.
        function featurizer(a, t)
            if any(isnan, a[end]) || abs(a[end,1]) > 100
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
        clusterspecs = Attractors.GroupViaClustering(optimal_radius_method = 1.)
        mapper = Attractors.AttractorsViaFeaturizing(ds, featurizer, clusterspecs; T = 500, threaded = true)
        continuation = GroupAcrossParameterContinuation(mapper; par_weight = 1.0)
        fractions_curves, attractors_info = Attractors.continuation(
            continuation, ps, pidx, sampler;
            samples_per_parameter = 1000, show_progress = false
        )

        for (i, p) in enumerate(ps)
            fs = fractions_curves[i]
            if p < 0.9
                k = sort!(collect(keys(fs)))
                @test length(k) == 2
            elseif p > 1.0 # (coexistence of period 1 and 3)
                k = sort!(collect(keys(fs)))
                @test length(k) == 3
            end
            @test sum(values(fs)) ≈ 1
        end

    end
end
