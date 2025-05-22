# This is an extensive test file. Its goal is to test all combination of mappers
# with all combinations of dynamical systems. However, this can take
# several years to complete. So, the environment parameter
# `ATTRACTORS_EXTENSIVE_TESTS` controls whether the tests should be done extensively or not.
# If not, a small, but representative subset of mappers and dynamical systems is used.

DO_EXTENSIVE_TESTS = get(ENV, "ATTRACTORS_EXTENSIVE_TESTS", "false") == "true"
using Test
using Attractors
using LinearAlgebra
using OrdinaryDiffEqVerner: Vern9
using Random
using Statistics

# Define generic testing framework
function test_basins(ds, u0s, grid, expected_fs_raw, featurizer;
        rerr = 1e-3, ferr = 1e-3, aerr = 1e-15, ε = nothing, max_distance = Inf,
        proximity_test = true, pairwise_comparison_matrix_test = false, featurizer_matrix = nothing,
        threshold_pairwise = 1,
        kwargs... # kwargs are propagated to recurrences
    )
    # u0s is Vector{Pair}
    sampler, = statespace_sampler(grid, 1234)
    ics = StateSpaceSet([copy(sampler()) for i in 1:1000])

    expected_fs = sort!(collect(values(expected_fs_raw)))
    known_ids = collect(u[1] for u in u0s)

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
        fs = basins_fractions(mapper, sampler; N = 100, show_progress = false)
        for k in keys(fs)
            @test 0 ≤ fs[k] ≤ 1
        end
        @test sum(values(fs)) ≈ 1 atol=1e-14

        # Precise test with known initial conditions
        fs, labels = basins_fractions(mapper, ics; show_progress = false)
        # @show nameof(typeof(mapper))
        # @show fs
        approx_atts = extract_attractors(mapper)
        found_fs = sort(collect(values(fs)))
        # @show found_fs
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
    end

    @testset "Recurrences" begin
        mapper = AttractorsViaRecurrences(ds, grid; kwargs...)
        test_basins_fractions(mapper; err = rerr)
    end


    @testset "Basins of Attraction" begin 
        # The reason of this test is to check whether 
        # basins_of_attraction does not modify the internal 
        # state of the mapper changing. (regression bug)
        function map_3_states(z, p, n)
            if z[1] < -1
                return SVector(-2.0, 0.) 
            elseif -1 ≤ z[1] ≤ 1
                return SVector(0.0, 0.)
            else 
                return SVector(2.0, 0.)
            end
        end
        ds = DiscreteDynamicalSystem(map_3_states, zeros(2)) 
        xg = yg = range(-3, 3; length = 20) 
        grid = (xg, yg)
        mapper = AttractorsViaRecurrences(ds, grid; sparse = false)  
        basins, atts = basins_of_attraction(mapper; show_progress = false)
        ics = [ [x,1.0] for x in range(-3,3,length = 20)]
        fractions, labels = basins_fractions(mapper, ics)

        @test 0 ∉ keys(fractions)
    end



    @testset "Featurizing, clustering" begin
        optimal_radius_method = "silhouettes_optim"
        config = GroupViaClustering(; num_attempts_radius=20, optimal_radius_method)
        mapper = AttractorsViaFeaturizing(ds, featurizer, config; Ttr = 500)
        test_basins_fractions(mapper;
            err = ferr, single_u_mapping = false, known_ids = [-1, 1, 2, 3]
        )
    end

    @testset "Featurizing, pairwise comparison" begin
        config = GroupViaPairwiseComparison(; threshold=threshold_pairwise,
        metric=Euclidean(), rescale_features=false)
        mapper = AttractorsViaFeaturizing(ds, featurizer, config; Ttr = 500)
        test_basins_fractions(mapper;
            err = ferr, single_u_mapping = false, known_ids = [-1, 1, 2, 3]
        )
    end

    if pairwise_comparison_matrix_test
        @testset "Featurizing, pairwise comparison, matrix features" begin
            function metric_hausdorff(A,B)
                set_distance(A, B, Hausdorff())
            end
            config = GroupViaPairwiseComparison(; threshold=threshold_pairwise,
            metric=metric_hausdorff, rescale_features=false)
            mapper = AttractorsViaFeaturizing(ds, featurizer_matrix, config; Ttr = 500)
            test_basins_fractions(mapper;
                err = ferr, single_u_mapping = false, known_ids = [-1, 1, 2, 3]
            )
        end
    end

    @testset "Featurizing, nearest feature" begin
        # First generate the templates
        function features_from_u(u)
            A, t = trajectory(ds, 100, u; Ttr = 1000, Δt = 1)
            featurizer(A, t)
        end
        # t = [features_from_u(x[2]) for x in u0s]
        # templates = Dict([u0[1] for u0 ∈ u0s] .=> t) # keeps labels of u0s
        templates = Dict(k => features_from_u(u) for (k, u) in u0s)

        config = GroupViaNearestFeature(templates; max_distance)
        mapper = AttractorsViaFeaturizing(ds, featurizer, config; Ttr=500)
        # test the functionality mapper(u0) -> label
        @test isinteger(mapper(current_state(ds))) == true
        test_basins_fractions(mapper; err = ferr, single_u_mapping = false)
    end

    if DO_EXTENSIVE_TESTS && proximity_test
        # Proximity method is the simplest and not crucial to test due to the limited
        # practical use of not being able to find attractors
        @testset "Proximity" begin
            known_attractors = Dict(
                k => trajectory(ds, 1000, v; Δt = 1, Ttr=100)[1] for (k,v) in u0s if k ≠ -1
            )
            mapper = AttractorsViaProximity(ds, known_attractors, ε; Ttr = 100, consecutive_lost_steps = 1000)
            test_basins_fractions(mapper; known = true, err = aerr)
        end
    end
end

# %% Actual tests
@testset "Analytic dummy map" begin
    function dumb_map(z, p, n)
        x, y = z
        r = p[1]
        if r < 0.5
            return SVector(0.0, 0.0)
        else
            if x ≥ 0
                return SVector(r, r)
            else
                return SVector(-r, -r)
            end
        end
    end

    r = 1.0
    ds = DeterministicIteratedMap(dumb_map, [0., 0.], [r])
    u0s = [1 => [r, r], 2 => [-r, -r]] # template ics

    xg = yg = range(-2.0, 2.0; length=100)
    grid = (xg, yg)
    expected_fs_raw = Dict(1 => 0.5, 1 => 0.5)
    featurizer(A, t) = SVector(A[1][1])
    test_basins(ds, u0s, grid, expected_fs_raw, featurizer;
    max_distance = 20, ε = 1e-1, proximity_test = false, threshold_pairwise=1,
    rerr = 1e-1, ferr = 1e-1, aerr = 1e-15)
end

@testset "Henon map: discrete & divergence" begin
    u0s = [1 => [0.0, 0.0], -1 => [0.0, 2.0]] # template ics
    henon_rule(x, p, n) = SVector{2}(1.0 - p[1]*x[1]^2 + x[2], p[2]*x[1])
    ds = DeterministicIteratedMap(henon_rule, zeros(2), [1.4, 0.3])
    xg = yg = range(-2.0, 2.0; length=100)
    grid = (xg, yg)
    expected_fs_raw = Dict(-1 => 0.575, 1 => 0.425)
    function featurizer(A, t)
        # Notice that unsupervised clustering cannot support "divergence to infinity",
        # which it identifies as another attractor (in fact, the first one).
        x = SVector(mean(A[:, 1]), mean(A[:, 2]))
        return any(isinf, x) ? SVector(200.0, 200.0) : x
    end
    test_basins(ds, u0s, grid, expected_fs_raw, featurizer;
    max_distance = 20, ε = 1e-3, proximity_test = false, threshold_pairwise=1)
end

@testset "Magnetic pendulum: projected system" begin
    mutable struct MagneticPendulumParams
        γs::Vector{Float64}
        d::Float64
        α::Float64
        ω::Float64
        magnets::Vector{SVector{2, Float64}}
    end
    function magnetic_pendulum_rule(u, p, t)
        x, y, vx, vy = u
        γs::Vector{Float64}, d::Float64, α::Float64, ω::Float64 = p.γs, p.d, p.α, p.ω
        dx, dy = vx, vy
        dvx, dvy = @. -ω^2*(x, y) - α*(vx, vy)
        for (i, ma) in enumerate(p.magnets)
            δx, δy = (x - ma[1]), (y - ma[2])
            D = sqrt(δx^2 + δy^2 + d^2)
            dvx -= γs[i]*(x - ma[1])/D^3
            dvy -= γs[i]*(y - ma[2])/D^3
        end
        return SVector(dx, dy, dvx, dvy)
    end
    function magnetic_pendulum(u = [sincos(0.12553*2π)..., 0, 0];
        γ = 1.0, d = 0.3, α = 0.2, ω = 0.5, N = 3, γs = fill(γ, N), diffeq)
        m = [SVector(cos(2π*i/N), sin(2π*i/N)) for i in 1:N]
        p = MagneticPendulumParams(γs, d, α, ω, m)
        return CoupledODEs(magnetic_pendulum_rule, u, p; diffeq)
    end

    diffeq = (alg = Vern9(), reltol = 1e-9, abstol = 1e-9)
    ds = magnetic_pendulum(γ=1, d=0.2, α=0.2, ω=0.8, N=3; diffeq)
    xg = range(-2,2; length = 201)
    yg = range(-2,2; length = 201)
    grid = (xg, yg)
    ds = ProjectedDynamicalSystem(ds, 1:2, [0.0, 0.0])
    u0s = [
        1 => [-0.5, 0.857],
        2 => [-0.5, -0.857],
        3 => [1.  , 0.],
    ]
    expected_fs_raw = Dict(2 => 0.314, 3 => 0.309, 1 => 0.377)
    function featurizer(A, t)
        return SVector(A[end][1], A[end][2])
    end

    function featurizer_matrix(A, t)
        return A
    end

    test_basins(ds, u0s, grid, expected_fs_raw, featurizer; ε = 0.2, Δt = 1.0, ferr=1e-2, featurizer_matrix, pairwise_comparison_matrix_test=true, threshold_pairwise=1)
end

# Okay, all of these aren't fundamentally new tests.
if DO_EXTENSIVE_TESTS

    @testset "Lorenz-84 system: interlaced close-by" begin # warning, super expensive test
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
        expected_fs_raw = Dict(2 => 0.18, 3 => 0.645, 1 => 0.175)

        using ComplexityMeasures

        function featurizer(A, t)
            # `g` is the number of boxes needed to cover the set
            probs = probabilities(ValueHistogram(0.1), A)
            g = exp(entropy(Renyi(; q = 0), probs))
            return SVector(g, minimum(A[:,1]))
        end

        test_basins(ds, u0s, grid, expected_fs_raw, featurizer;
        ε = 0.01, ferr=1e-2, Δt = 0.2, consecutive_attractor_steps = 5, Ttr = 100, threshold_pairwise=100) #threshold is very high because features haven't really converged yet here
    end

    @testset "Duffing oscillator: stroboscopic map" begin
        @inbounds function duffing_rule(x, p, t)
            ω, f, d, β = p
            dx1 = x[2]
            dx2 = f*cos(ω*t) - β*x[1] - x[1]^3 - d * x[2]
            return SVector(dx1, dx2)
        end
        diffeq = (alg = Vern9(), reltol = 1e-9, abstol = 1e-9)
        ds = CoupledODEs(duffing_rule, [0.1, 0.25], [1.0, 0.2, 0.15, -1]; diffeq)

        xg = yg = range(-2.2, 2.2; length=200)
        grid = (xg, yg)
        T = 2π/1.0
        ds = StroboscopicMap(ds, T)
        u0s = [
            1 => [-2.2, -2.2],
            2 => [2.2, 2.2],
        ]
        expected_fs_raw = Dict(2 => 0.488, 1 => 0.512)
        function featurizer(A, t)
            return SVector(A[end][1], A[end][2])
        end

        test_basins(ds, u0s, grid, expected_fs_raw, featurizer;
        ε = 1.0, ferr=1e-2, rerr = 1e-2, aerr = 5e-3, threshold_pairwise=1)
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
        ds = CoupledODEs(thomas_rule, [1.0, 0, 0], [0.1665]; diffeq=(reltol=1e-9,))

        xg = yg = range(-6.0, 6.0; length = 100) # important, don't use 101 here, because
        # the dynamical system has some fixed points ON the hyperplane.
        grid = (xg, yg)
        pmap = PoincareMap(ds, (3, 0.0);
            Tmax = 1e6,
            rootkw = (xrtol = 1e-8, atol = 1e-8),
        )
        u0s = [
            1 => [1.83899, -4.15575, 0],
            2 => [1.69823, -0.0167188, 0],
            3 => [-4.08547,  -2.26516, 0],
        ]
        expected_fs_raw = Dict(2 => 0.28, 3 => 0.266, 1 => 0.454)
        function thomas_featurizer(A, t)
            x, y = columns(A)
            return SVector(minimum(x), minimum(y))
        end

        test_basins(pmap, u0s, grid, expected_fs_raw, thomas_featurizer; ε = 1.0, ferr=1e-2)
    end

end
