using Attractors
using Test
using OrdinaryDiffEq
using ProgressMeter

include("../../src/boundaries/edgetracking.jl")

@testset "Saddle point of cubic map" begin
    cubicmap(u, p, n) = SVector{1}(p[1]*u[1] - u[1]^3)    
    ds = DeterministicIteratedMap(cubicmap, [1.0], [2.0])
    attrs = Dict(1 => StateSpaceSet([1.0]), 2 => StateSpaceSet([-1.0]))
    saddle = edgetracking(ds, attrs; Δt=1, abstol=1e-8).edge[end]
    @test saddle[1] <= 1e-5
end

@testset "Saddle point of FitzHugh-Nagumo system" begin
    fhn(u, p, t) = SVector{2}([10*(u[1] - u[1]^3 - u[2]), -3*u[2] + u[1]])
    ds = CoupledODEs(fhn, ones(2), diffeq=(;alg = Vern9()))
    fp = [sqrt(2/3), sqrt(2/27)]
    attrs = Dict([1 => StateSpaceSet([fp]), 2 => StateSpaceSet([-fp])])
    bisect_thresh, diverge_thresh, maxiter, abstol = 1e-8, 1e-7, 100, 1e-8
    edge = edgetracking(ds, attrs; u1=[-1.0, 0.2], u2=[1.0, 0.2],
        bisect_thresh, diverge_thresh, maxiter, abstol).edge
    println(edge[end], edge[end-1])
    @test sqrt(sum(abs, (edge[end]-zeros(2)).^2)) <= 1e-6
end

@testset "Thomas' rule" begin
    
end

function thomas_rule(u, p, t)
    x,y,z = u
    b = p[1]
    xdot = sin(y) - b*x
    ydot = sin(z) - b*y
    zdot = sin(x) - b*z
    return SVector{3}(xdot, ydot, zdot)
end
ds = CoupledODEs(thomas_rule, [1.0, 0, 0], [0.16]; diffeq=(reltol=1e-12,))
xg = yg = zg = range(-6.0, 6.0; length = 101); grid = (xg, yg, zg)
mapper = AttractorsViaRecurrences(ds, grid; consecutive_recurrences = 1000)
sampler, = statespace_sampler(grid)
basins_fractions(mapper, sampler)
attractors = extract_attractors(mapper)