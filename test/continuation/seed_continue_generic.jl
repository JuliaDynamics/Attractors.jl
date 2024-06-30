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

    r = 1.0
    ds = DeterministicIteratedMap(dumb_map, [0., 0.], [r])
    yg = xg = range(-10., 10, length = 100)
    grid = (xg,yg)
    mapper = AttractorsViaRecurrences(ds, grid; sparse = true, show_progress = false)

    sampler, = statespace_sampler(grid, 1234)

    rrange = range(0, 2; length = 20)
    ridx = 1

    function featurizer(A, t)
        return A[end]
    end

    # in all honesty, we don't have to test 2 grouping configs,
    # as the algorithm is agnostic to the grouping. But oh well!
    group1 = GroupViaClustering(optimal_radius_method = 0.1)
    group2 = GroupViaPairwiseComparison(threshold = 0.1)

    @testset "grouping: $(nameof(typeof(group)))" for group in (group1, group2)
        mapper = AttractorsViaFeaturizing(ds, featurizer, group; Ttr = 2, T = 2)
        algo = AttractorSeedContinueMatch(mapper)
        fractions_cont, a = global_continuation(
            algo, rrange, ridx, sampler;
            show_progress = false, samples_per_parameter = 1000
        )

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