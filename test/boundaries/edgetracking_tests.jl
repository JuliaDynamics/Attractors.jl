using Attractors
using Test
using OrdinaryDiffEq

@testset "Saddle point of FitzHugh-Nagumo system" begin
    fhn(u, p, t) = SVector{2}([10*(u[1] - u[1]^3 - u[2]), -3*u[2] + u[1]])
    ds = CoupledODEs(fhn, ones(2), diffeq=(;alg = Vern9()))
    fp = [sqrt(2/3), sqrt(2/27)]
    attrs = Dict([1 => StateSpaceSet([fp]), 2 => StateSpaceSet([-fp])])
    bisect_thresh, diverge_thresh, maxiter, abstol = 1e-7, 1e-6, 100, 1e-8
    edge, track1, track2 = edgetracking(ds, attrs; u1=[-1.0, 0.2], u2=[1.0, 0.2],
        bisect_thresh, diverge_thresh, maxiter, abstol)
    @test sqrt(sum(abs, (edge[end]-zeros(2)).^2)) <= 1e-8
end