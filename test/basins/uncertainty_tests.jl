using Attractors
using Test

function grebogi_map()

    function grebogi_map_rule(u, p, n)
        θ = u[1]; x = u[2]
        a, b, J₀ = p
        dθ = θ + a * sin(2 * θ) - b * sin(4 * θ) - x * sin(θ)
        dθ = mod(dθ, 2π) # to avoid problems with attractor at θ=π
        dx = -J₀ * cos(θ)
        return SVector{2}(dθ, dx)
    end
    p = (a = 1.32, b = 0.9, J₀ = 0.3)
    return DeterministicIteratedMap(grebogi_map_rule, rand(2), p)
end


@testset "original paper" begin
    ds = grebogi_map()
    θg = range(0, 2π, length = 251)
    xg = range(-0.5, 0.5, length = 251)
    mapper = AttractorsViaRecurrences(ds, (θg, xg); sparse = false)
    bsn, att = basins_of_attraction(mapper; show_progress = false)
    e, f, α = uncertainty_exponent(bsn; range_ε = 3:15)
    # In the paper the value is roughly 0.2
    @test (0.2 ≤ α ≤ 0.3)
end

@testset "Newton map" begin

    function newton_map(dz, z, p, n)
        f(x) = x^p[1] - 1
        df(x) = p[1] * x^(p[1] - 1)
        z1 = z[1] + im * z[2]
        dz1 = f(z1) / df(z1)
        z1 = z1 - dz1
        dz[1] = real(z1)
        dz[2] = imag(z1)
        return
    end

    ds = DiscreteDynamicalSystem(newton_map, [0.1, 0.2], [3])
    xg = yg = range(-1.0, 1.0, length = 300)
    mapper = AttractorsViaRecurrences(ds, (xg, yg); sparse = false)
    bsn, att = basins_of_attraction(mapper; show_progress = false)
    e, f, α = uncertainty_exponent(bsn; range_ε = 5:30)

    # Value (published) from the box-counting dimension is 1.42. α ≃ 0.6
    @test (0.55 ≤ α ≤ 0.65)

end

@testset "Basin entropy and Fractal test" begin
    ds = grebogi_map()
    θg = range(0, 2π, length = 300)
    xg = range(-0.5, 0.5, length = 300)
    mapper = AttractorsViaRecurrences(ds, (θg, xg); sparse = false)
    basin, attractors = basins_of_attraction(mapper; show_progress = false)
    Sb, Sbb = basin_entropy(basin, 6)
    @test 0.4 ≤ Sb ≤ 0.42
    @test 0.6 ≤ Sbb ≤ 0.61

    test_res, Sbb = basins_fractal_test(basin; ε = 5)
    @test test_res == :fractal

    henon_rule(x, p, n) = SVector{2}(1.0 - p[1] * x[1]^2 + x[2], p[2] * x[1])
    henon() = DeterministicIteratedMap(henon_rule, zeros(2), [1.4, 0.3])
    ds = henon()
    xg = yg = range(-2.0, 2.0, length = 300)
    mapper = AttractorsViaRecurrences(ds, (xg, yg); sparse = false)
    basin, attractors = basins_of_attraction(mapper; show_progress = false)
    test_res, Sbb = basins_fractal_test(basin; ε = 5)
    @test test_res == :smooth
end

@testset "Basin entropy API" begin
    basins = rand(Int, 50, 50)
    @test_throws ArgumentError basin_entropy(basins, 7)
end
