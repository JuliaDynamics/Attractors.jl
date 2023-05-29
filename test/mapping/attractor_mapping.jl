# This is an extensive test file. Its goal is to test all combination of mappers
# with all combinations of dynamical systems. However, this can take
# several years to complete. So, the environment parameter
# `ATTRACTORS_EXTENSIVE_TESTS` controls whether the tests should be done extensively or not.
# If not, a small, but representative subset of mappers and dynamical systems is used.

DO_EXTENSIVE_TESTS = get(ENV, "ATTRACTORS_EXTENSIVE_TESTS", "false") == "true"

using Test
using Attractors
using LinearAlgebra
using OrdinaryDiffEq: Vern9
using Random
using Statistics
using PredefinedDynamicalSystems: Systems

# Define generic testing framework
function test_basins(ds, u0s, grid, expected_fs_raw, featurizer;
        rerr = 1e-3, ferr = 1e-3, aerr = 1e-15, ε = nothing, max_distance = Inf,
        kwargs... # kwargs are propagated to recurrences
    )
    # u0s is Vector{Pair}
    known_attractors = Dict(
        k => trajectory(ds, 10000, v; Δt = 1, Ttr=100)[1] for (k,v) in u0s if k ≠ -1
    )
    sampler, = statespace_sampler(Random.MersenneTwister(1234);
        min_bounds = minimum.(grid), max_bounds = maximum.(grid)
    )
    ics = StateSpaceSet([sampler() for i in 1:1000])
    expected_fs = sort!(collect(values(expected_fs_raw)))
    known_ids = collect(u[1] for u in u0s)
    reduced_grid = map(g -> range(minimum(g), maximum(g); length = 10), grid)

    # reusable testing function
    function test_basins_fractions(mapper;
            err = 1e-3, known=false, single_u_mapping = true,
            known_ids = known_ids, expected_fs = expected_fs,
        )
        if single_u_mapping
            for (k, u0) in u0s
                @test k == mapper(u0)
            end
        end
        # Generic test
        fs = basins_fractions(mapper, sampler; show_progress = false)
        approx_atts_sampler = extract_attractors(mapper)
        for k in keys(fs)
            @test 0 ≤ fs[k] ≤ 1
        end
        @test sum(values(fs)) ≈ 1 atol=1e-14

        # Precise test with known initial conditions
        fs, labels = basins_fractions(mapper, ics; show_progress = false)
        approx_atts = extract_attractors(mapper)
        found_fs = sort(collect(values(fs)))
        if length(found_fs) > length(expected_fs)
            # drop -1 key if it corresponds to just unidentified points
            found_fs = found_fs[2:end]
        end
        @test length(found_fs) == length(expected_fs) #number of attractors
        errors = abs.(expected_fs .- found_fs)
        for er in errors
            @test er .≤ err
        end
        if known # also test whether the attractor index is correct
            for k in known_ids
                @test abs(fs[k] - expected_fs_raw[k]) ≤ err
            end
        end
        
        @test length(approx_atts_sampler) == length(approx_atts)
        @test all(approx_atts_sampler[i] == approx_atts[i] for i in eachindex(approx_atts_sampler) )
        
        # `basins_of_attraction` tests
        basins, approx_atts = basins_of_attraction(mapper, reduced_grid; show_progress=false)
        @test length(size(basins)) == length(grid)
        if known
            bids = sort!(unique(basins))
            @test all(x -> x ∈ known_ids, bids)
        end
    end

    @testset "Recurrences" begin
        mapper = AttractorsViaRecurrences(ds, grid; show_progress = false, kwargs...)
        test_basins_fractions(mapper; err = rerr)
    end

    @testset "Featurizing, clustering" begin
        optimal_radius_method = "silhouettes_optim"
        config = GroupViaClustering(; num_attempts_radius=20, optimal_radius_method)
        mapper = AttractorsViaFeaturizing(ds, featurizer, config; Ttr = 500)
        test_basins_fractions(mapper;
            err = ferr, single_u_mapping = false, known_ids = [-1, 1, 2, 3]
        )
    end

    @testset "Featurizing, nearest feature" begin
        # First generate the templates
        function features_from_u(u)
            A, t = trajectory(ds, 100, u; Ttr = 500, Δt = 1)
            featurizer(A, t)
        end
        t = [features_from_u(x[2]) for x in u0s]
        templates = Dict([u0[1] for u0 ∈ u0s] .=> t) # keeps labels of u0s

        config = GroupViaNearestFeature(templates; max_distance)
        mapper = AttractorsViaFeaturizing(ds, featurizer, config; Ttr=500)
        test_basins_fractions(mapper; err = ferr, single_u_mapping = false)
    end

    if DO_EXTENSIVE_TESTS
        # Proximity method is the simplest and not crucial to test due to the limited
        # practical use of not being able to find attractors
        @testset "Proximity" begin
            mapper = AttractorsViaProximity(ds, known_attractors, ε; Ttr = 100)
            test_basins_fractions(mapper; known = true, err = aerr)
        end
    end
end

# Actual tests
@testset "Henon map: discrete & divergence" begin
    u0s = [1 => [0.0, 0.0], -1 => [0.0, 2.0]] # template ics
    henon_rule(x, p, n) = SVector{2}(1.0 - p[1]*x[1]^2 + x[2], p[2]*x[1])
    henon() = DeterministicIteratedMap(henon_rule, zeros(2), [1.4, 0.3])
    ds = henon()

    xg = yg = range(-2.0, 2.0; length=100)
    grid = (xg, yg)
    expected_fs_raw = Dict(1 => 0.451, -1 => 0.549)
    function featurizer(A, t)
        # Notice that unsupervised clustering cannot support "divergence to infinity",
        # which it identifies as another attractor (in fact, the first one).
        x = SVector(mean(A[:, 1]), mean(A[:, 2]))
        return any(isinf, x) ? SVector(200.0, 200.0) : x
    end
    test_basins(ds, u0s, grid, expected_fs_raw, featurizer;
    max_distance = 20, ε = 1e-3)
end

# Okay, all of these aren't fundamentally new tests.
if DO_EXTENSIVE_TESTS

    @testset "Lorenz-84 system: interlaced close-by" begin
        F = 6.886; G = 1.347; a = 0.255; b = 4.0
        function lorenz84_rule(u, p, t)
            F, G, a, b = p
            x, y, z = u
            dx = -y^2 -z^2 -a*x + a*F
            dy = x*y - y - b*x*z + G
            dz = b*x*y + x*z - z
            return SVector{3}(dx, dy, dz)
        end
        diffeq = (alg = Vern9(), reltol = 1e-9, abstol = 1e-9)
        ds = CoupledODEs(lorenz84_rule, fill(0.1, 3), [F, G, a, b]; diffeq)

        u0s = [
            1 => [2.0, 1, 0], # periodic
            2 => [-2.0, 1, 0], # chaotic
            3 => [0, 1.5, 1.0], # fixed point
        ]
        M = 200; z = 3
        xg = yg = zg = range(-z, z; length = M)
        grid = (xg, yg, zg)
        expected_fs_raw = Dict(2 => 0.165, 3 => 0.642, 1 => 0.193)

        using ComplexityMeasures

        function featurizer(A, t)
            # `g` is the number of boxes needed to cover the set
            probs = probabilities(ValueHistogram(0.1), A)
            g = exp(entropy(Renyi(; q = 0), probs))
            return SVector(g, minimum(A[:,1]))
        end

        test_basins(ds, u0s, grid, expected_fs_raw, featurizer;
        ε = 0.01, ferr=1e-2, Δt = 0.2, mx_chk_att = 20)
    end

    @testset "Duffing oscillator: stroboscopic map" begin
        @inbounds function duffing_rule(x, p, t)
            ω, f, d, β = p
            dx1 = x[2]
            dx2 = f*cos(ω*t) - β*x[1] - x[1]^3 - d * x[2]
            return SVector(dx1, dx2)
        end
        ds = CoupledODEs(duffing_rule, [0.1, 0.25], [1.0, 0.2, 0.15, -1])

        xg = yg = range(-2.2, 2.2; length=200)
        grid = (xg, yg)
        diffeq = (alg = Vern9(), reltol = 1e-9, abstol = 1e-9)
        T = 2π/1.0
        ds = stroboscopicmap(ds, T; diffeq)
        u0s = [
            1 => [-0.8, 0],
            2 => [1.8, 0],
        ]
        expected_fs_raw = Dict(2 => 0.511, 1 => 0.489)
        function featurizer(A, t)
            return SVector(A[end][1], A[end][2])
        end

        test_basins(ds, u0s, grid, expected_fs_raw, featurizer;
        ε = 0.01, ferr=1e-2, rerr = 1e-2, aerr = 5e-3)
    end


    @testset "Magnetic pendulum: projected system" begin
        # TODO: replace this!
        ds = Systems.magnetic_pendulum(γ=1, d=0.2, α=0.2, ω=0.8, N=3)
        xg = range(-2,2,length = 201)
        yg = range(-2,2,length = 201)
        grid = (xg, yg)
        diffeq = (alg = Vern9(), reltol = 1e-9, abstol = 1e-9)
        ds = projected_integrator(ds, 1:2, [0.0, 0.0]; diffeq)
        u0s = [
            1 => [-0.5, 0.857],
            2 => [-0.5, -0.857],
            3 => [1.  , 0.],
        ]
        expected_fs_raw = Dict(2 => 0.318, 3 => 0.347, 1 => 0.335)

        function featurizer(A, t)
            return SVector(A[end][1], A[end][2])
        end

        test_basins(ds, u0s, grid, expected_fs_raw, featurizer; ε = 0.2, Δt = 1.0, ferr=1e-2)
    end


    @testset "Thomas cyclical: Poincaré map" begin
        function thomas_rule(u, p, t)
            x,y,z = u
            b = p[1]
            xdot = sin(y) - b*x
            ydot = sin(z) - b*y
            zdot = sin(x) - b*z
            return SVector{3}(xdot, ydot, zdot)
        end
        ds = CoupledODEs(thomas_rule, [1.0, 0, 0], [0.1665])

        xg = yg = range(-6.0, 6.0; length = 100) # important, don't use 101 here, because
        # the dynamical system has some fixed points ON the hyperplane.
        grid = (xg, yg)
        pmap = poincaremap(ds, (3, 0.0), 1e6;
            rootkw = (xrtol = 1e-8, atol = 1e-8), diffeq=(reltol=1e-9,)
        )
        u0s = [
            1 => [1.83899, -4.15575, 0],
            2 => [1.69823, -0.0167188, 0],
            3 => [-4.08547,  -2.26516, 0],
        ]
        expected_fs_raw = Dict(2 => 0.29, 3 => 0.237, 1 => 0.473)
        function thomas_featurizer(A, t)
            x, y = columns(A)
            return SVector(minimum(x), minimum(y))
        end

        test_basins(pmap, u0s, grid, expected_fs_raw, thomas_featurizer; ε = 1.0, ferr=1e-2)
    end
end