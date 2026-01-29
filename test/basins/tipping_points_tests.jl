using Attractors
using Test

@testset "basins_fractions(::Array)" begin
    b₋ = rand(Random.MersenneTwister(1234), Int16.(1:3), (10, 10, 10))
    fs = basins_fractions(b₋)
    @test keytype(fs) == Int16
    @test all(v -> 0.31 < v < 0.35, values(fs))
    # Also test analytically resolved just to be sure
    ba = [1 2; 2 1]
    fs = basins_fractions(ba)
    @test fs[1] == 0.5
    @test fs[2] == 0.5
end

@testset "Tipping probabilities (overlap)" begin
    basins_before = ones(20, 20)
    basins_before[1:10, 1:10] .= -1
    basins_after = ones(20, 20)
    basins_after[:, 11:20] .= 2
    P = tipping_probabilities(basins_before, basins_after)
    @test size(P) == (2, 2)
    @test P[1, 1] ≈ 1 / 3
    @test P[2, 1] == 1
    @test P[1, 2] ≈ 2 / 3
    @test P[2, 2] == 0
end
