using Attractors
using Test
using OrdinaryDiffEq:Vern9
using Random
using Statistics

function morris_lecar!(du, u, p, t)
    V1 = -0.00; V2 = 0.15; V4 = 0.1; VCa = 1 ;  VL = -0.5; VK = -0.7; gCa = 1.2; gK = 2; gL = 0.5; τ = 3;
    I, V3 = p
    V, N = u
    M(x) = 0.5*(1 + tanh((x-V1)/V2))
    G(x) = 0.5*(1 + tanh((x-V3)/V4))
    du[1] = -gCa*M(V)*(V - VCa) -gK*N*(V - VK) -gL*(V-VL) + I
    du[2] = 1/τ*(-N + G(V))
end

@testset "Density of attractors" begin
    function test_morrislecar(; stop_at_dt, store_once_per_cell)
        I = 0.115; V3 = 0.1; p = [I, V3];
        df = ContinuousDynamicalSystem(morris_lecar!, [0.0, 0.0], p)
        diffeq = (reltol = 1e-9,  alg = Vern9())
        xg = yg = range(-1,1,length = 20000)
        mapper = AttractorsViaRecurrences(df, (xg, yg);
                mx_chk_fnd_att = 5000,
                mx_chk_loc_att = 1000,
                stop_at_dt,
                store_once_per_cell,
                sparse = true,
                Δt = 0.1,
                Ttr = 10,
                )

        sampler, = Attractors.statespace_sampler(Random.MersenneTwister(1);
            min_bounds = [-0.5, -1], max_bounds = [0.5, 1])
        ics = Dataset([sampler() for i in 1:1000])

        fs, labels, atts = basins_fractions(mapper, ics; show_progress=false)
        num_att = length(atts)
        return num_att
    end

    @testset "Dense limit cycles" begin #More points are stored, and they sample the limit cycle more evenly
        num_att = test_morrislecar(; stop_at_dt=true, store_once_per_cell=false)
        @test num_att == 1
    end

    @testset "Sparse limit cycles" begin #Finds two attractors instead of one!
        num_att = test_morrislecar(; stop_at_dt=false, store_once_per_cell=true)
        @test num_att == 2
    end
end