using Attractors
using Test
using OrdinaryDiffEqVerner
using LinearAlgebra

@testset "Saddle point of cubic map" begin
    cubicmap(u, p, n) = SVector{1}(p[1]*u[1] - u[1]^3)
    ds = DeterministicIteratedMap(cubicmap, [1.0], [2.0])
    attrs = Dict(1 => StateSpaceSet([1.0]), 2 => StateSpaceSet([-1.0]))
    saddle = edgetracking(ds, attrs; Δt=1, abstol=1e-8).edge[end]
    @test saddle[1] < 1e-5
end

@testset "Saddle point of FitzHugh-Nagumo system" begin
    fhn(u, p, t) = SVector{2}([10*(u[1] - u[1]^3 - u[2]), -3*u[2] + u[1]])
    ds = CoupledODEs(fhn, ones(2), diffeq=(;alg = Vern9()))
    fp = [sqrt(2/3), sqrt(2/27)]
    attrs = Dict([1 => StateSpaceSet([fp]), 2 => StateSpaceSet([-fp])])
    bisect_thresh, diverge_thresh, maxiter, abstol = 1e-8, 1e-7, 100, 1e-9
    edge = edgetracking(ds, attrs; u1=[-1.0, 0.2], u2=[1.0, 0.2],
        bisect_thresh, diverge_thresh, maxiter, abstol, Ttr = 100).edge
    @test sqrt(sum(abs, (edge[end]-zeros(2)).^2)) < 1e-5
end

@testset "Thomas' rule" begin
    # Chaotic dynamical system
    function thomas_rule(u, p, t)
        x,y,z = u
        b = p[1]
        xdot = sin(y) - b*x
        ydot = sin(z) - b*y
        zdot = sin(x) - b*z
        return SVector{3}(xdot, ydot, zdot)
    end
    ds = CoupledODEs(thomas_rule, [1.0, 0, 0], [0.16]; diffeq=(reltol=1e-12,))

    # Find attractors on a 3D grid
    xg = yg = yz = range(-6.0, 6.0; length = 101)
    grid = (xg, yg, yz)
    mapper = AttractorsViaRecurrences(ds, grid; consecutive_recurrences = 1000)
    sampler, = statespace_sampler(grid)
    basins_fractions(mapper, sampler)
    attractors = extract_attractors(mapper)

    # Run edgetracking between pairs of points lying on different attractors
    n_sample = 25
    pairs12, pairs13, pairs23 = [], [], []
    for i in 1:Int(sqrt(n_sample))
        _pairs12, _pairs13, _pairs23 = [], [], []
        for j in 1:Int(sqrt(n_sample))
            et12 = edgetracking(ds, attractors;
                    u1=attractors[1][i], u2=attractors[2][j], bisect_thresh=1e-4, Ttr = 100,
                    diverge_thresh=1e-3, maxiter=10000, abstol=1e-5, verbose=false)
            et13 = edgetracking(ds, attractors;
                    u1=attractors[1][i], u2=attractors[3][j], bisect_thresh=1e-4, Ttr = 100,
                    diverge_thresh=1e-3, maxiter=10000, abstol=1e-5, verbose=false)
            et23 = edgetracking(ds, attractors;
                    u1=attractors[2][i], u2=attractors[3][j], bisect_thresh=1e-4, Ttr = 100,
                    diverge_thresh=1e-3, maxiter=10000, abstol=1e-5, verbose=false)

            et12.success ? push!(_pairs12, et12.edge[end]) : nothing
            et13.success ? push!(_pairs13, et13.edge[end]) : nothing
            et23.success ? push!(_pairs23, et23.edge[end]) : nothing
        end
        push!(pairs12, _pairs12)
        push!(pairs13, _pairs13)
        push!(pairs23, _pairs23)
    end
    edgestates = reduce(vcat, [pairs12 pairs13 pairs23])

    # Verify that all found edge states have the same Euclidean norm of `norm_value`
    # (due to the symmetry of the system `ds`)
    norm_value = 4.06585
    norm_deviations = [norm(edgestates[i]) - norm_value for i in 1:length(edgestates)]
    @test maximum(norm_deviations) < 1e-3
end
