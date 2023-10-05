using Attractors
using DynamicalSystemsBase

fhn(u, p, t) = SVector{2}([10*(u[1] - u[1]^3 - u[2]), -3*u[2] + u[1]])
ds = CoupledODEs(fhn, ones(2))
fp = [sqrt(2/3), sqrt(2/27)]
attrs = Dict([1 => StateSpaceSet([fp]), 2 => StateSpaceSet([-fp])])

attrs[2]

 pds = ParallelDynamicalSystem(ds, [Matrix(attrs[1]), Matrix(attrs[2])])

Matrix(attrs[1])

ParallelDynamicalSystem(ds, [[0.,0], [1.,0]])

@testset "edgetracking (FitzHugh-Nagumo)" begin
    
    ds = CoupledODEs(fhn, ones(2), diffeq=(;alg = Vern9()))
    fp = [sqrt(2/3), sqrt(2/27)]
    attrs = Dict([1 => StateSpaceSet([fp]), 2 => StateSpaceSet([-fp])])
    abstol, maxiter, eps1, eps2 = 1e-9, 500, 1e-7, 1e-6
    M = edgetracking(ds, [-1.0, 0.3], [1.0, 0.3], attrs, abstol=abstol, maxiter=maxiter, eps1=eps1, eps2=eps2, output_level=0)
    @test sqrt(sum(abs, (M-zeros(2)).^2)) <= 1e-8
end