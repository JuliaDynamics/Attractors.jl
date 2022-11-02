using Attractors
using Test

@testset "basins_fractions(::Array)" begin
    b₋ = rand(Random.MersenneTwister(1234), Int16.(1:3), (10, 10, 10))
    fs = basins_fractions(b₋)
    @test keytype(fs) == Int16
    @test all(v -> 0.31 < v < 0.35, values(fs))
    # Also test analytically resolved juuuuust to be sure
    ba = [1 2; 2 1]
    fs = basins_fractions(ba)
    @test fs[1] == 0.5
    @test fs[2] == 0.5
end

@testset "Tipping points" begin
    xg = yg = range(-4, 4, length = 100)
    ds = Systems.magnetic_pendulum(d=0.2, α=0.2, ω=0.8, N=3)
    pinteg = projected_integrator(ds, 1:2, [0,0]; diffeq=(reltol = 1e-9,))
    mapper = AttractorsViaRecurrences(pinteg, (xg, yg); sparse = false, Δt = 1.0)
    basins_before, = basins_of_attraction(mapper; show_progress = false)
    pinteg.integ.p.γs = [1, 1, 0.05] # Change magnet strength
    mapper = AttractorsViaRecurrences(pinteg, (xg, yg); sparse = false, Δt = 1.0)
    basins_after, = basins_of_attraction(mapper; show_progress = false)
    P = tipping_probabilities(basins_before, basins_after)
    @test size(P) == (3,2)
    @test all(0.495 .≤ P[1, :] .≤ 0.505) # must be ≈ 50%
    for i in 1:3; @test sum(P[i, :]) ≈ 1; end
end
