using Attractors
using Test

function predator_prey_ds()
    function predator_prey_fastslow(u, p, t)
        α, γ, ϵ, ν, h, K, m = p
        N, P = u
        du1 = α*N*(1 - N/K) - γ*N*P / (N+h)
        du2 = ϵ*(ν*γ*N*P/(N+h) - m*P)
        return SVector(du1, du2)
    end
    γ = 2.5
    h = 1
    ν = 0.5
    m = 0.4
    ϵ = 1.0
    α = 0.8
    K = 15
    u0 = rand(2)
    p0 = [α, γ, ϵ, ν, h, K, m]
    ds = CoupledODEs(predator_prey_fastslow, u0, p0)
    return ds
end

@testset "Different grids" begin
    @testset "Irregular grid" begin
        # when pow = 1 we have regular grid that makes
        # the unstable fixed point be an attractor due to trajectory slow-down
        # The second initial condition is directly on the limit cycle
        u0 = [9.0, 3.5]
        u0 = [0.03, 0.5]
        As = []
        for pow in (1, 8)
            ds = predator_prey_ds()
            d = 21
            xg = range(0, 18.0^(1/pow); length = d).^pow
            yg = range(0, 18.0; length = d)
            mapper = AttractorsViaRecurrences(ds, (xg, yg);
                Dt = 0.05, sparse = true, force_non_adaptive = true,
                mx_chk_fnd_att = 20, mx_chk_loc_att = 100,
                store_once_per_cell = true,
                mx_chk_safety = 1000,
            )
            id = mapper(u0)
            A = extract_attractors(mapper)[1]
            push!(As, A)
            if pow == 1
                @test length(A) == 1 # fixed point
            elseif pow == 2
                @test length(A) > 1 # limit cycle
            end
        end
    end

    @testset "Subdivision based grid" begin
        u0 = [0.03, 0.5]
        ds = predator_prey_ds()
        xg = yg = range(0, 18; length = 21)
        grid = (xg, yg)
        gridsd = subdivision_based_grid(ds, (xg, yg))
        for (i, g) in enumerate((grid, gridsd))
            # we use the same mapper as the irregular grid (important!)
            mapper = AttractorsViaRecurrences(ds, g;
                Dt = 0.05, sparse = true, force_non_adaptive = true,
                mx_chk_fnd_att = 20, mx_chk_loc_att = 100,
                store_once_per_cell = true,
                mx_chk_safety = 1000,
            )
            id = mapper(u0)
            A = extract_attractors(mapper)[1]
            if i == 1
                @test length(A) == 1
            else
                @test length(A) > 1
            end
        end
    end
end

@testset "basin cell index" begin
    ds = predator_prey_ds()
    xg = yg = range(0, 5, length = 6)

    grid = subdivision_based_grid(ds, (xg, yg); maxlevel = 3)
    lvl_array = grid.lvl_array

    @test (Attractors.basin_cell_index((2, 1.1) , grid) == CartesianIndex((17,9)))
    @test (Attractors.basin_cell_index((3.6, 2.1) , grid) == CartesianIndex((29,17)))
    @test (Attractors.basin_cell_index((4.9, 4.9) , grid) == CartesianIndex((33,33)))
    @test (Attractors.basin_cell_index((0.5, 0.5) , grid) == CartesianIndex((5,5)))
    @test (Attractors.basin_cell_index((3.4, 1.6) , grid) == CartesianIndex((27,13)))
    @test (Attractors.basin_cell_index((4.7, 2.2) , grid) == CartesianIndex((37,17)))
end

@testset "automatic Δt" begin
    ds = predator_prey_ds()

    xg = yg = range(0, 18, length = 30)
    grid0 = (xg, yg)
    xg = yg = range(0, 18.0^(1/2); length = 20).^2
    grid1 = (xg, yg)
    grid2 = SubdivisionBasedGrid(grid0, rand([0, 1, 2], (30, 30)))

    using Attractors: RegularGrid, IrregularGrid

    Dt0 = automatic_Δt_basins(ds, RegularGrid(grid0))
    Dt1 = automatic_Δt_basins(ds, IrregularGrid(grid1))
    Dt2 = automatic_Δt_basins(ds, grid2)

    @test Dt0 > 0
    @test Dt1 > 0
    @test Dt2 > 0
end