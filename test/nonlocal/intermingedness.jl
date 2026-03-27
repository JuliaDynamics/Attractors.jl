using Test, Attractors
using Random
using Distances: WeightedEuclidean

@testset "intermingledness def" begin
    # generate features that are not intermingled in 2D space
    # but are intermingled in their x coordinate
    randompoint() = (rand([-1, 1]) + 0.05randn(), 0.05randn())
    features = [randompoint() for _ in 1:1000]
    features = [SVector(x, y + (x > 0 ? rand([-0.5, 0.5]) : 0)) for (x, y) in features]
    # sort by x dimension so that key 1 is the left cluster
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

end



