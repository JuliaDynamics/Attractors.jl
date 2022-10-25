using Attractors
using Test
using OrdinaryDiffEq:Vern9
using Random
using Statistics

@testset "Non-sparse version" begin
    function test_compatibility_sparse_nonsparse(ds, grid; diffeq, kwargs...)
        sampler, = statespace_sampler(Random.MersenneTwister(1234);
            min_bounds = minimum.(grid), max_bounds = maximum.(grid)
        )
            ics = Dataset([sampler() for i in 1:1000])
            mapper = AttractorsViaRecurrences(ds, grid; sparse=true, diffeq, show_progress = false, kwargs...)
            fs_sparse, labels_sparse, approx_atts_sparse = basins_fractions(mapper, ics; show_progress = false)

            mapper = AttractorsViaRecurrences(ds, grid; sparse=false, diffeq, show_progress = false, kwargs...)
            fs_non, labels_non, approx_atts_non = basins_fractions(mapper, ics; show_progress = false)

            @test fs_sparse == fs_non
            @test labels_sparse == labels_non
            @test approx_atts_sparse == approx_atts_non
        nothing
    end
    @testset "Lorenz84" begin
        F = 6.886; G = 1.347; a = 0.255; b = 4.0
        ds = Systems.lorenz84(; F, G, a, b)
        u0s = [
            1 => [2.0, 1, 0], # periodic
            2 => [-2.0, 1, 0], # chaotic
            3 => [0, 1.5, 1.0], # fixed point
        ]
        diffeq = (alg = Vern9(), reltol = 1e-9, abstol = 1e-9)
        M = 200; z = 3; xg = yg = zg = range(-z, z; length = M); grid = (xg, yg, zg)
        test_compatibility_sparse_nonsparse(ds, grid; diffeq)
    end
end