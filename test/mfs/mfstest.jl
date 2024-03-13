using Attractors
using Test
using LinearAlgebra

@testset "Newton 2d" begin

    function newton_map(z, p, n)
        z1 = z[1] + im*z[2]
        dz1 = newton_f(z1, p[1])/newton_df(z1, p[1])
        z1 = z1 - dz1
        return SVector(real(z1), imag(z1))
    end
    newton_f(x, p) = x^p - 1
    newton_df(x, p)= p*x^(p-1)

    ds = DiscreteDynamicalSystem(newton_map, [0.1, 0.2], [3.0])
    xg = yg = range(-1.5, 1.5; length = 400)

    newton = AttractorsViaRecurrences(ds, (xg, yg);
        sparse = false, consecutive_lost_steps = 1000
    )

    attractors = [[1.0, 0.0], [-0.5, 0.8660254037844386], [-0.5, -0.8660254037844386]]
    algo_r = Attractors.MFSBruteForce(sphere_decrease_factor = 0.96)

    randomised = Dict([atr => minimal_fatal_shock(newton, atr, [(-1.5, 1.5)], algo_r) for atr in attractors])

    algo_bb = Attractors.MFSBlackBoxOptim()
    blackbox = Dict([atr => (minimal_fatal_shock(newton, atr, [(-1.5, 1.5)], algo_bb)) for atr in attractors])

    random_seed = [[rand([-1,1])*rand()/2, rand([-1,1])*rand()/2] for _ in 1:20]
    randomised_r = Dict([atr => (minimal_fatal_shock(newton, atr, [(-1.5, 1.5)], algo_r)) for atr in random_seed])
    blackbox_r = Dict([atr => (minimal_fatal_shock(newton, atr, [(-1.5, 1.5)], algo_bb)) for atr in random_seed])

    test = true
    for i in (keys(randomised))
        if norm(randomised[i]) >= 0.65 || norm(randomised[i]) <= 0.60 || newton(randomised[i] + i) == newton(i)
            test = false
        end
    end
    @test test

    test = true
    for i in (keys(blackbox))
        if norm(blackbox[i]) >= 0.629 || norm(blackbox[i]) <= 0.62009 || newton(blackbox[i] + i) == newton(i)
            test = false
        end
    end
    @test test

    test = true
    for i in (keys(randomised_r))
        if norm(randomised_r[i]) >= 0.5 || newton(randomised_r[i] + i) == newton(i)
            test = false
        end
    end
    @test test

    test = true
    for i in (keys(blackbox_r))
        if norm(blackbox_r[i]) >= 0.5 || newton(blackbox_r[i] + i) == newton(i)
            test = false
        end
    end
    @test test

    @testset "target id" begin
        atrdict = Dict(i => StateSpaceSet([a]) for (i, a) in enumerate(attractors))
        mapper = AttractorsViaProximity(ds, atrdict, 0.01)
        algo_bb = Attractors.MFSBlackBoxOptim()
        # multiple target attractors
        mfs_1 = [
            i => minimal_fatal_shock(mapper, a, (-1.5, 1.5), algo_bb; target_id = setdiff(1:3, [i]))
            for (i, a) in enumerate(attractors)
        ]
        # specific target attractors
        mfs_2 = [
            i => minimal_fatal_shock(mapper, a, (-1.5, 1.5), algo_bb; target_id = setdiff(1:3, [i])[1])
            for (i, a) in enumerate(attractors)
        ]
        # due to the symmetries of the system all mfs should be identical in magnitude
        for mfs in (mfs_1, mfs_2)
            for f in mfs
                @test norm(f[2]) ≈ 0.622 rtol = 1e-2
            end
        end
    end
end



###############################################
#           Magnetic 2D                       #
###############################################
@testset "Magnetic 2D" begin

struct MagneticPendulum
    magnets::Vector{SVector{2, Float64}}
end
mutable struct MagneticPendulumParams2
    γs::Vector{Float64}
    d::Float64
    α::Float64
    ω::Float64
end

function (m::MagneticPendulum)(u, p, t)
    x, y, vx, vy = u
    γs::Vector{Float64}, d::Float64, α::Float64, ω::Float64 = p.γs, p.d, p.α, p.ω
    dx, dy = vx, vy
    dvx, dvy = @. -ω^2*(x, y) - α*(vx, vy)
    for (i, ma) in enumerate(m.magnets)
        δx, δy = (x - ma[1]), (y - ma[2])
        D = sqrt(δx^2 + δy^2 + d^2)
        dvx -= γs[i]*(x - ma[1])/D^3
        dvy -= γs[i]*(y - ma[2])/D^3
    end
    return SVector(dx, dy, dvx, dvy)
end

function magnetic_pendulum(u = [sincos(0.12553*2π)..., 0, 0];
    γ = 1.0, d = 0.3, α = 0.2, ω = 0.5, N = 3, γs = fill(γ, N))
    m = MagneticPendulum([SVector(cos(2π*i/N), sin(2π*i/N)) for i in 1:N])
    p = MagneticPendulumParams2(γs, d, α, ω)
    return CoupledODEs(m, u, p)
end

ds = magnetic_pendulum(d=0.2, α=0.2, ω=0.8, N=3)

psys = ProjectedDynamicalSystem(ds, [1, 2], [0.0, 0.0])

attractors_m = Dict(i => StateSpaceSet([dynamic_rule(ds).magnets[i]]) for i in 1:3)

mapper_m = AttractorsViaProximity(psys, attractors_m)

xg = yg = range(-4, 4; length = 201)
grid = (xg, yg)
algo_r = Attractors.MFSBruteForce(sphere_decrease_factor = 0.99)
algo_bb = Attractors.MFSBlackBoxOptim()

attractor3 = vec((collect(values(attractors_m)))[3])
attractor2 = vec((collect(values(attractors_m)))[2])
attractor1 = vec((collect(values(attractors_m)))[1])

randomised_r = Dict([atr => norm(minimal_fatal_shock(mapper_m, atr, [(-4, 4), (-4, 4)], algo_r))
                                         for atr in [attractor1[1], attractor2[1], attractor3[1]]])

blackbox_r = Dict([atr => norm(minimal_fatal_shock(mapper_m, atr, [(-4, 4), (-4, 4)], algo_bb))
                                            for atr in [attractor1[1], attractor2[1], attractor3[1]]])

@test map(x -> (x <= 0.4) && (x) > 0.39, values(randomised_r)) |> all
@test map(x -> (x <= 0.395) && (x) > 0.39, values(blackbox_r)) |> all

end




###############################################
#           Thomas 3D                       #
###############################################
@testset "3D symmetry" begin
    using LinearAlgebra: norm
    function thomas_rule(u, p, t)
        x,y,z = u
        b = p[1]
        xdot = sin(y) - b*x
        ydot = sin(z) - b*y
        zdot = sin(x) - b*z
        return SVector{3}(xdot, ydot, zdot)
    end

    thomas_cyclical(u0 = [1.0, 0, 0]; b = 0.2) = CoupledODEs(thomas_rule, u0, [b])

    ds = thomas_cyclical(b = 0.1665)
    xg = yg = zg = range(-6.0, 6.0; length = 251)
    mapper_3d = AttractorsViaRecurrences(ds, (xg, yg, zg))

    ux = SVector(1.5, 0, 0)
    uy = SVector(0, 1.5, 0)

    algo_bb = Attractors.MFSBlackBoxOptim(max_steps = 50000)

    ux_res = minimal_fatal_shock(mapper_3d, ux, (-6.0,6.0), algo_bb)
    uy_res = minimal_fatal_shock(mapper_3d, uy, (-6.0,6.0), algo_bb)

    @test norm(ux_res) - norm(uy_res) < 0.0001
end
