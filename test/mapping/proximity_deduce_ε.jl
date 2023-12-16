using Attractors
using Test

@testset "Proximity deduce ε" begin
    henon_rule(x, p, n) = SVector{2}(1.0 - p[1]*x[1]^2 + x[2], p[2]*x[1])
    henon() = DeterministicIteratedMap(henon_rule, zeros(2), [1.4, 0.3])
    ds = henon()
    @testset "single attractor, no ε" begin
        attractors = Dict(1 => trajectory(ds, 10000, [0.0, 0.0]; Δt = 1, Ttr=100)[1])
        mapper = AttractorsViaProximity(ds, attractors)
        @test trunc(mapper.ε, digits = 2) ≈ 0.18 # approximate size of attractor here
    end
    @testset "two attractors, analytically known ε" begin
        attractors = Dict(
            1 => StateSpaceSet([0 1.0]; warn = false),
            2 => StateSpaceSet([0 2.0]; warn = false)
        )
        mapper = AttractorsViaProximity(ds, attractors)
        @test mapper.ε == 0.5
    end
    @testset "one attractor, single point (invalid)" begin
        attractors = Dict(
            1 => StateSpaceSet([0 1.0]; warn = false),
        )
        @test_throws ArgumentError AttractorsViaProximity(ds, attractors)
    end
end

@testset "Fix #61" begin
    cubicmap(u, p, n) = SVector{1}(p[1]*u[1] - u[1]^3)
    ds = DeterministicIteratedMap(cubicmap, [1.0], [2.0])
    fp = [sqrt(2)]
    attrs = Dict(1 => StateSpaceSet([fp]), 2 => StateSpaceSet([-fp]))
    mapper = AttractorsViaProximity(ds, attrs)
    label = mapper([2.0])
    @test label == -1
end