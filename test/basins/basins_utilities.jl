using Attractors
using Test

@testset "Basin types" begin
    # also test that type of IC can be different than attractors eltype
    grid = (range(-1, 1; length = 3), range(-1, 1; length = 3), )
    basins = rand(1:2, (3, 3))
    attractors = Dict(k => StateSpaceSet([rand(2) for _ in 1:5]) for k in 1:2)
    BoA = ArrayBasinsOfAttraction(basins, attractors, grid)
    @test BoA isa ArrayBasinsOfAttraction

    points = [rand(-1:1, (2, )) for _ in 1:10]
    basins = rand(1:2, 10)
    BoA = SampledBasinsOfAttraction(basins, attractors, points)
    @test BoA isa SampledBasinsOfAttraction
end

@testset "convergence_and_basins_of_attraction(mapper::BasinMap)" begin

    # Dynamical system
    test_ds(u, p, t) = SVector(0.5 * u[1] * (1 - u[1]^2))

    # Test continuous system
    ds = CoupledODEs(test_ds, [1.0])
    grid = (range(0, 1; length = 5),)
    mapper = BasinMapRecurrences(ds, grid)

    basins, attractors, convergence = convergence_and_basins_of_attraction(mapper, grid; show_progress = true)

    @test eltype(convergence) == typeof(current_time(ds))

    # Test iterated map
    ds = DeterministicIteratedMap(test_ds, [1.0])
    mapper = BasinMapRecurrences(ds, grid)
    basins, attractors, convergence = convergence_and_basins_of_attraction(mapper, grid; show_progress = true)

    @test eltype(convergence) == typeof(current_time(ds))

end
