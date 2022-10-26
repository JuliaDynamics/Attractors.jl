using Attractors
using Test
using Random
using Statistics

@testset "Compatibility sparse and nonsparse" begin
    function test_compatibility_sparse_nonsparse(ds, grid; diffeq=NamedTuple(), kwargs...)
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

    @testset "Henon map: discrete & divergence" begin
        ds = Systems.henon(zeros(2); a = 1.4, b = 0.3)
        u0 = [0.0, 0.0]
        xg = yg = range(-2.0, 2.0; length=100)
        grid = (xg, yg)
        test_compatibility_sparse_nonsparse(ds, grid)
    end

end