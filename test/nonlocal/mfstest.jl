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
    algo_r = Attractors.MCSBruteForce(sphere_decrease_factor = 0.96)

    randomised = Dict([atr => minimal_critical_shock(newton, atr, (-1.5, 1.5), algo_r) for atr in attractors])

    algo_bb = Attractors.MCSBlackBoxOptim()
    blackbox = Dict([atr => (minimal_critical_shock(newton, atr, (-1.5, 1.5), algo_bb)) for atr in attractors])

    random_seed = [[rand([-1,1])*rand()/2, rand([-1,1])*rand()/2] for _ in 1:20]
    randomised_r = Dict([atr => (minimal_critical_shock(newton, atr, (-1.5, 1.5), algo_r)) for atr in random_seed])
    blackbox_r = Dict([atr => (minimal_critical_shock(newton, atr, (-1.5, 1.5), algo_bb)) for atr in random_seed])

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
        algo_bb = Attractors.MCSBlackBoxOptim()
        # multiple target attractors
        mfs_1 = [
            i => minimal_critical_shock(mapper, a, (-1.5, 1.5), algo_bb; target_id = setdiff(1:3, [i]))
            for (i, a) in enumerate(attractors)
        ]
        # specific target attractors
        mfs_2 = [
            i => minimal_critical_shock(mapper, a, (-1.5, 1.5), algo_bb; target_id = setdiff(1:3, [i])[1])
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

attractors_dict = Dict(i => StateSpaceSet([dynamic_rule(ds).magnets[i]]) for i in 1:3)
attractors = collect(values(attractors_dict))

# CRUCIAL: Because we are using the projected magnetic pendulum, we need to make sure
# that we have transient time (because it passes above a magnet but with nonzero velocity)
mapper_m = AttractorsViaProximity(psys, attractors_dict; Ttr = 100)

searchrange =  [(-4, 4), (-4, 4)]

algo_r = MCSBruteForce(sphere_decrease_factor = 0.99, seed = 42, sphere_iterations = 1000, initial_iterations = 1000)

randomised_r = [
    norm(minimal_critical_shock(mapper_m, A[1], searchrange, algo_r))
    for A in attractors
]

# Due to symmetry, all MCS should be the same, and by plotting the basins
# you can also see how large they should be.
@test all(x -> isapprox(x, 0.395; atol = 1e-2), randomised_r)

algo_bb = MCSBlackBoxOptim()
randomised_bb = [
    norm(minimal_critical_shock(mapper_m, A[1], searchrange, algo_bb))
    for A in attractors
]

# the black box case is drastically more accurate
@test all(x -> isapprox(x, 0.3923; atol = 1e-4), randomised_bb)

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

    algo_bb = Attractors.MCSBlackBoxOptim(max_steps = 50000)

    ux_res = minimal_critical_shock(mapper_3d, ux, (-6.0, 6.0), algo_bb)
    uy_res = minimal_critical_shock(mapper_3d, uy, (-6.0, 6.0), algo_bb)

    # Due to the symmetry of the system, the shocks have to be the same
    @test norm(ux_res) - norm(uy_res) < 0.0001
end
