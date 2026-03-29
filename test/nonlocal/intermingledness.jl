using Test, Attractors
using Random
using Distances: WeightedEuclidean

@testset "intermingledness def" begin
    # generate features that are not intermingled in 2D space
    # but are intermingled in their x coordinate
    randompoint() = (rand([-1, 1]) + 0.05randn(), 0.05randn())
    features = [randompoint() for _ in 1:1000]
    features = [SVector(x, y + (x > 0 ? rand([-0.5, 0.5]) : 0)) for (x, y) in features]
    # sort by x dimension so that key 1 is the left cluster always
    sort!(features)
    features = StateSpaceSet(features)

    gconfig = GroupViaPairwiseComparison(; rescale_features = false, threshold = 0.5)

    labels = group_features(features, gconfig)
    replace!(labels, 1 => 5)

    # scatter(features; color = labels)

    mingle = intermingledness(features, labels)
    @test sort!(collect(keys(mingle))) == [2, 3, 5]
    @test all(<(0.1), values(mingle))
    @test mingle[5] < mingle[2]
    @test mingle[2] ≈ mingle[3] atol = 1e-2 rtol = 1e-2

    # now multidim
    distances = [
        Euclidean(),
        WeightedEuclidean([1, 0]) # distance of x dimension only
    ]

    mingles = intermingledness(features, labels, distances)

    @test mingles[1] == mingle

    xmingle = mingles[2]
    @test xmingle[5] < 0.1
    @test xmingle[2] ≈ xmingle[3] atol = 1e-1 rtol = 1e-1
    @test xmingle[2] ≈ 1 atol = 1e-1 rtol = 1e-1

end

@testset "intermingledness accumulator" begin
    # we'll use a modified dumb map that is intermingled when projected to second dimension
    function dumb_map(z, p, n)
        x, y = z
        r = p[1]
        if r < 0.5
            return SVector(0.0, 0.0)
        else
            if x ≥ 0
                return SVector(r, 0)
            else
                return SVector(-r, 0)
            end
        end
    end
    dynamics = DiscreteDynamicalSystem(dumb_map, [0.0, 1.0], [0.9])
    grid = (-1:0.1:1, -1:0.1:1)
    mapper = AttractorsViaRecurrences(dynamics, grid; sparse = false, Ttr = 2)
    distances = [
        Euclidean(),
        WeightedEuclidean([0, 1]) # distance of x dimension only
    ]
    accumulator = StabilityMeasuresAccumulator(mapper; idistances = distances)

    A = ics_from_grid(grid)
    for u0 in A
        id = accumulator(u0) # run this to accumulate measures
    end
    results = finalize_accumulator(accumulator)

    @test haskey(results, "intermingledness1")
    @test haskey(results, "intermingledness2")

    i1 = results["intermingledness1"]
    i2 = results["intermingledness2"]

    @test i1[1] ≈ i1[2] atol = 1e-1
    @test i1[1] ≈ 0.66 atol = 1e-1
    @test i2[1] ≈ i2[2] atol = 1e-14
    @test i2[1] ≈ 1 atol = 1e-14
end
