# Testing the generalized `AttractorSeedContinueMatch`.
# We don't test recurrences combo because it is already tested
# in another file
using Test, Attractors
using Random

@testset "analytic bistable: featurizing" begin
    # This is a fake bistable map that has two equilibrium points
    # for r > 0.5. It has analytically resolved fractions for any box.
    function dumb_map(dz, z, p, n)
        x, y = z
        r, q = p

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

    r = 1.0; q = 0.5
    ds = DeterministicIteratedMap(dumb_map, [0., 0.], [r, q])
    yg = xg = range(-10., 10, length = 100)
    grid = (xg,yg)
    mapper1 = AttractorsViaRecurrences(ds, grid; sparse = true, show_progress = false)

    sampler, = statespace_sampler(grid, 1234)

    rrange = range(0, 2; length = 20)
    ridx = 1

    function featurizer(A, t)
        return A[end]
    end

    # in all honesty, we don't have to test 2 grouping configs,
    # as the algorithm is agnostic to the grouping. But oh well!
    group1 = GroupViaClustering(optimal_radius_method = 0.1)
    mapper2 = AttractorsViaFeaturizing(ds, featurizer, group1; Ttr = 2, T = 2)

    group2 = GroupViaPairwiseComparison(threshold = 0.1)
    mapper3 = AttractorsViaFeaturizing(ds, featurizer, group2; Ttr = 2, T = 2)

    mappers = [mapper1, mapper2, mapper3, mapper1]

    @testset "case: $(i)" for (i, mapper) in enumerate(mappers)
        algo = AttractorSeedContinueMatch(mapper)

        if i < 4
            fractions_cont, a = global_continuation(
                algo, rrange, ridx, sampler;
                show_progress = false, samples_per_parameter = 1000
            )
        else # test parameter curve version
            pcurve = [[1 => r, 2 => 1.1] for r in rrange]
            fractions_cont, a = global_continuation(
                algo, pcurve, sampler;
                show_progress = false, samples_per_parameter = 1000
            )
        end


        for (i, r) in enumerate(rrange)

            fs = fractions_cont[i]
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

end