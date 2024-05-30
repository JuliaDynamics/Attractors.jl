# This file also tests `aggregate_attractor_fractions`!
using Test, Attractors
using Random

@testset "analytic bistable map" begin
    # This is a fake bistable map that has two equilibrium points
    # for r > 0.5. It has predictable fractions.
    function dumb_map(dz, z, p, n)
        x,y = z
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
    featurizer(A,t) = A[end]
    grouping_config = GroupViaPairwiseComparison(; threshold=0.2)
    mapper = AttractorsViaFeaturizing(ds, featurizer, grouping_config)

    yg = xg = range(-10., 10, length = 100)
    grid = (xg,yg)
    sampler, = statespace_sampler(grid, 1234)
    samples_per_parameter = 1000
    ics = Dataset([deepcopy(sampler()) for _ in 1:samples_per_parameter])

    rrange = range(0, 2; length = 20)
    ridx = 1
    fsc = FeaturizingFindAndMatch(mapper; threshold = 0.3)
    fractions_curves, a = continuation(
        fsc, rrange, ridx, ics;
        show_progress = false
    )

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
