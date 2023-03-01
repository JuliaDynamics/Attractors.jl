DO_EXTENSIVE_TESTS = get(ENV, "ATTRACTORS_EXTENSIVE_TESTS", "false") == "true"

if DO_EXTENSIVE_TESTS
# The functionality tested here has been resolved and is only added as a test
# for future security. It has no need to be tested in every commit.

using Attractors
using Test
using OrdinaryDiffEq: Vern9
using Random

@testset "point saving limit cycles" begin

    @inbounds function morris_lecar_rule(u, p, t)
        I, V3, V1, V2, V4, VCa, VL, VK, gCa, gK, gL, τ = p
        V, N = u
        M = 0.5*(1 + tanh((V-V1)/V2))
        G = 0.5*(1 + tanh((V-V3)/V4))
        du1 = -gCa*M*(V - VCa) -gK*N*(V - VK) -gL*(V-VL) + I
        du2 = 1/τ*(-N + G)
        return SVector{2}(du1, du2)
    end

    function test_morrislecar(; store_once_per_cell, stop_at_Δt, Δt = 0.15)
        u0=[0.1, 0.1];
        p = (I = 0.115, V3 = 0.1, V1 = -0.00, V2 = 0.15, V4 = 0.1,
        VCa = 1 ,  VL = -0.5, VK = -0.7, gCa = 1.2, gK = 2,
        gL = 0.5, τ = 3)
        # p = [I, V3, V1, V2, V4, VCa, VL, VK, gCa, gK, gL, τ]
        diffeq = (reltol = 1e-9, alg = Vern9(), adaptive = !stop_at_Δt, dt = Δt)

        df = CoupledODEs(morris_lecar_rule, u0, p; diffeq)

        xg = yg = range(-1,1,length = 2000)
        mapper = AttractorsViaRecurrences(df, (xg, yg);
            mx_chk_loc_att = 1000,
            mx_chk_fnd_att = 2,
            store_once_per_cell,
            sparse = true,
            Δt,
            Ttr = 10,
        )

        sampler, = Attractors.statespace_sampler(Random.MersenneTwister(1);
            min_bounds = [-0.5, 0], max_bounds = [0.5, 1])
        ics = StateSpaceSet([sampler() for i in 1:1000])

        fs, atts, labels = basins_fractions(mapper, ics; show_progress=false)
        num_att = length(atts)
        return num_att
    end

    @testset "Sparse limit cycles" begin # when using a
        num_att = test_morrislecar(; store_once_per_cell=true, stop_at_Δt = false)
        @test num_att > 1
        num_att = test_morrislecar(; store_once_per_cell=false, stop_at_Δt = false)
        @test num_att > 1
        num_att = test_morrislecar(; store_once_per_cell=false, stop_at_Δt = true)
        @test num_att == 1
        num_att = test_morrislecar(; store_once_per_cell=true, stop_at_Δt = true)
        @test num_att == 1
    end
end

@testset "Compatibility sparse and nonsparse" begin
    function test_compatibility_sparse_nonsparse(ds, grid; kwargs...)
        sampler, = statespace_sampler(Random.MersenneTwister(1234);
            min_bounds = minimum.(grid), max_bounds = maximum.(grid)
        )
            ics = StateSpaceSet([sampler() for i in 1:1000])

            mapper = AttractorsViaRecurrences(ds, grid; sparse=true, show_progress = false, kwargs...)
            fs_sparse, approx_atts_sparse, labels_sparse = basins_fractions(mapper, ics; show_progress = false)

            mapper = AttractorsViaRecurrences(ds, grid; sparse=false, show_progress = false, kwargs...)
            fs_non, approx_atts_non, labels_non = basins_fractions(mapper, ics; show_progress = false)

            @test fs_sparse == fs_non
            @test labels_sparse == labels_non
            @test approx_atts_sparse == approx_atts_non
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