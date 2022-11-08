DO_EXTENSIVE_TESTS = parse(Bool, get(ENV, "ATTRACTORS_EXTENSIVE_TESTS", false))

if DO_EXTENSIVE_TESTS
# The functionality tested here has been resolved and is only added as a test
# for future security. It has no need to be tested in every commit.

using Attractors
using Test
using OrdinaryDiffEq: Vern9
using Random

@testset "Density of attractors" begin
    function test_morrislecar(; stop_at_Δt, store_once_per_cell)
        df = Systems.morris_lecar(; I = 0.115)
        diffeq = (reltol = 1e-9,  alg = Vern9())
        xg = yg = range(-1,1,length = 2000)
        mapper = AttractorsViaRecurrences(df, (xg, yg);
                mx_chk_loc_att = 1000,
                mx_chk_fnd_att = 2,
                stop_at_Δt,
                diffeq,
                store_once_per_cell,
                sparse = true,
                Δt = 0.15,
                Ttr = 10,
                )

        sampler, = Attractors.statespace_sampler(Random.MersenneTwister(1);
            min_bounds = [-0.5, 0], max_bounds = [0.5, 1])
        ics = Dataset([sampler() for i in 1:1000])

        fs, labels, atts = basins_fractions(mapper, ics; show_progress=false)
        num_att = length(atts)
        return num_att
    end

    @testset "Dense limit cycles" begin #More points are stored, and they sample the limit cycle more evenly
        num_att = test_morrislecar(; stop_at_Δt=true, store_once_per_cell=false)
        @test num_att == 1
    end

    @testset "Sparse limit cycles" begin #Finds two attractors instead of one!
        num_att = test_morrislecar(; stop_at_Δt=false, store_once_per_cell=true)
        @test num_att > 1
    end
end

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

end # extensive tests clause