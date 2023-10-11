# This file performs more extensive tests specifically for `AttractorsViaRecurrences`
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
        diffeq = (reltol = 1e-9, alg = Vern9())

        df = CoupledODEs(morris_lecar_rule, u0, p; diffeq)

        xg = yg = range(-1,1,length = 2000)
        mapper = AttractorsViaRecurrences(df, (xg, yg);
            mx_chk_loc_att = 1000,
            mx_chk_fnd_att = 2,
            store_once_per_cell,
            sparse = true,
            Δt,
            Ttr = 10,
            force_non_adaptive = stop_at_Δt,
        )

        sampler, = statespace_sampler(HRectangle([-0.5, 0], [0.5, 1]), 155)
        ics = StateSpaceSet([copy(sampler()) for i in 1:1000])

        fs, labels = basins_fractions(mapper, ics; show_progress=false)
        num_att = length(fs)
        return num_att
    end

    @testset "Sparse limit cycles" begin
        # how many times we store per cell doesn't change anything
        # however breaking the orbit commensurability does!
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
            sampler, = statespace_sampler(grid, 1244)
            ics = StateSpaceSet([copy(sampler()) for i in 1:1000])

            mapper = AttractorsViaRecurrences(ds, grid; sparse=true, show_progress = false, kwargs...)
            fs_sparse, labels_sparse = basins_fractions(mapper, ics; show_progress = false)
            approx_atts_sparse = extract_attractors(mapper)
            mapper = AttractorsViaRecurrences(ds, grid; sparse=false, show_progress = false, kwargs...)
            fs_non, labels_non = basins_fractions(mapper, ics; show_progress = false)
            approx_atts_non = extract_attractors(mapper)

            @test fs_sparse == fs_non
            @test labels_sparse == labels_non
            @test approx_atts_sparse == approx_atts_non
    end

    @testset "Henon map: discrete & divergence" begin
        henon_rule(x, p, n) = SVector{2}(1.0 - p[1]*x[1]^2 + x[2], p[2]*x[1])
        ds = DeterministicIteratedMap(henon_rule, zeros(2), [1.4, 0.3])
        xg = yg = range(-2.0, 2.0; length=100)
        grid = (xg, yg)
        test_compatibility_sparse_nonsparse(ds, grid)
    end
end

@testset "Escape to -1 test" begin
# This is for testing if the chk safety keyword is working
# as intended. The output should be only -1.

function dissipative_standard_map_rule(u, p, n)
    x, y = u
    ν, f₀ = p
    s = x + y
    xn = mod2pi(s)
    yn = (1 - ν)*y + f₀*(sin(s))
    return SVector(xn, yn)
end

p0 = (ν = 0.02, f0 = 4.0)
u0 = [0.1, 0.1]
ds = DeterministicIteratedMap(dissipative_standard_map_rule, u0, p0)
density = 10
xg = range(0, 2π; length = density+1)[1:end-1]
ymax = 2
yg = range(-ymax, ymax; length = density)
grid = (xg, yg)

mapper_kwargs = (
    mx_chk_safety = 10,
    sparse = false, # we want to compute full basins
)

mapper = AttractorsViaRecurrences(ds, grid; mapper_kwargs...)
basins, attractors = basins_of_attraction(mapper)
ids = sort!(unique(basins))
@test ids... == -1

end

@testset "continuing lost state" begin
    # see discussion in https://github.com/JuliaDynamics/Attractors.jl/pull/103#issuecomment-1754617229
    henon_rule(x, p, n) = SVector{2}(1.0 - p[1]*x[1]^2 + x[2], p[2]*x[1])
    henon() = DeterministicIteratedMap(henon_rule, zeros(2), [1.4, 0.3])
    ds = henon()
    xg = yg = range(-1.0, 1.0; length=100)
    grid = (xg, yg)
    mapper = AttractorsViaRecurrences(ds, grid)
    id = mapper([0. ,0.])
    @test id == 1 # we resume going into the attractor outside the grid
    id2 = mapper([10. ,10.])
    @test id2 == -1 # actual divergence to infinity is still detected
end

end # extensive tests clause
