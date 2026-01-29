using Attractors
using Test

@testset "convergence_and_basins_of_attraction(mapper::AttractorMapper)" begin

    # Dynamical system
    test_ds(u, p, t) = SVector(0.5 * u[1] * (1 - u[1]^2))

    # Test continuous system
    ds = CoupledODEs(test_ds, [1.0])
    grid = (range(0, 1; length = 5),)
    mapper = AttractorsViaRecurrences(ds, grid)

    basins, attractors, convergence = convergence_and_basins_of_attraction(mapper, grid; show_progress = true)

    @test eltype(convergence) == typeof(current_time(ds))

    # Test iterated map
    ds = DeterministicIteratedMap(test_ds, [1.0])
    mapper = AttractorsViaRecurrences(ds, grid)
    basins, attractors, convergence = convergence_and_basins_of_attraction(mapper, grid; show_progress = true)

    @test eltype(convergence) == typeof(current_time(ds))

end
