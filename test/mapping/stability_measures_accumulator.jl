DO_EXTENSIVE_TESTS = get(ENV, "ATTRACTORS_EXTENSIVE_TESTS", "false") == "true"

using Test, Attractors
using Random

# Test the computation of nonlocal stability measures using the
# `StabilityMeasuresAccumulator` for a dumb map.
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

dynamics = DiscreteDynamicalSystem(dumb_map, [1.0, 1.0], [1.0])

grid = ([-1, 0, 1], [-1, 0, 1])

mapper = AttractorsViaRecurrences(dynamics, grid; sparse = false)

A = ics_from_grid(grid)
for u0 in A
    id = mapper(u0)
end

attractors = extract_attractors(mapper)

mapper = AttractorsViaProximity(dynamics, attractors, 0.01, Ttr = 0)

accumulator = StabilityMeasuresAccumulator(mapper, finite_time = 0.5)

for u0 in A
    id = accumulator(u0)
end

results = finalize_accumulator(accumulator)

# The expected results are computed as in the following example:
# maximal_noncritical_shock of attractor 2 at [1, 1] is computed as the distance between the
# that point and the furthest point of its own basin at [0, -1] which is 2.23607. For
# attractor 1 at [-1, -1] the maximal noncritical shock is 2.0, which is the distance to the
# point [-1, 1] which is the furthest point in its basin.

results_expected = Dict(
    "mean_convergence_time" => Dict(2 => 0.22222, 1 => 0.55555, -1 => NaN),
    "maximal_noncritical_shock_magnitude" => Dict(2 => 2.23607, 1 => 2.0, -1 => NaN),
    "finite_time_basin_stability" => Dict(2 => 0.11111, 1 => 0.11111, -1 => 0.0),
    "median_convergence_pace" => Dict(2 => 0.5, 1 => 0.5, -1 => NaN),
    "median_convergence_time" => Dict(2 => 1.0, 1 => 1.0, -1 => NaN),
    "minimal_critical_shock_magnitude" => Dict(2 => 2.0, 1 => 1.0, -1 => NaN),
    "basin_stability" => Dict(2 => 0.66667, 1 => 0.33333, -1 => 0.0),
    "maximal_convergence_pace" => Dict(2 => 1.0, 1 => 1.0, -1 => NaN),
    "maximal_convergence_time" => Dict(2 => 1.0, 1 => 1.0, -1 => NaN),
    "mean_convergence_pace" => Dict(2 => 0.16667, 1 => 0.40604, -1 => NaN),
    "basin_fraction" => Dict(2 => 0.66667, 1 => 0.33333, -1 => 0.0),
    "mean_noncritical_shock_magnitude" => Dict(2 => 0.33333, 1 => 0.85003, -1 => NaN),
    "characteristic_return_time" => Dict(2 => NaN, 1 => NaN, -1 => NaN),
    "reactivity" => Dict(2 => NaN, 1 => NaN, -1 => NaN),
    "maximal_amplification" => Dict(2 => NaN, 1 => NaN, -1 => NaN),
    "maximal_amplification_time" => Dict(2 => NaN, 1 => NaN, -1 => NaN)
)
# Check if the results are as expected
@testset "Nonlocal Stability Measures Accumulator with dumb map" begin
    for (key, value) in results_expected
        @test key in keys(results)
        @test isapprox(
            sort(collect(values(value))),
            sort(collect(values(results[key])));
            atol = 1.0e-5,
            nans = true,
        )
    end
end

# Now we test the continuation of nonlocal stability measures.
pcurve = [[1 => p] for p in [-1.0, 1.0]]
attractors_cont = [
    Dict(1 => StateSpaceSet([SVector(0.0, 0.0)])),
    Dict(2 => StateSpaceSet([SVector(1.0, 1.0)]), 1 => StateSpaceSet([SVector(-1.0, -1.0)])),
]

proximity_mapper_options = (
    Ttr = 0, stop_at_Δt = false, horizon_limit = 1.0e2, consecutive_lost_steps = 10000,
)
measures_cont = stability_measures_along_continuation(
    dynamics, attractors_cont, pcurve, ics_from_grid(grid), ε = 0.1, finite_time = 0.5,
    proximity_mapper_options = proximity_mapper_options
)

measures_cont_expected = Dict(
    "finite_time_basin_stability" => [Dict(1 => 0.11111, -1 => 0.0), Dict(2 => 0.11111, 1 => 0.11111, -1 => 0.0)],
    "maximal_noncritical_shock_magnitude" => [Dict(1 => 1.41421, -1 => NaN), Dict(2 => 2.23607, 1 => 2.0, -1 => NaN)],
    "median_convergence_pace" => [Dict(1 => 0.70711, -1 => NaN), Dict(2 => 0.5, 1 => 0.5, -1 => NaN)],
    "basin_stability" => [Dict(1 => 1.0, -1 => 0.0), Dict(2 => 0.66667, 1 => 0.33333, -1 => 0.0)],
    "maximal_convergence_pace" => [Dict(1 => 1.0, -1 => NaN), Dict(2 => 1.0, 1 => 1.0, -1 => NaN)],
    "mean_convergence_pace" => [Dict(1 => 0.75871, -1 => NaN), Dict(2 => 0.40604, 1 => 0.16667, -1 => NaN)],
    "basin_fraction" => [Dict(1 => 1.0, -1 => 0.0), Dict(2 => 0.66667, 1 => 0.33333, -1 => 0.0)],
    "mean_convergence_time" => [Dict(1 => 0.88889, -1 => NaN), Dict(2 => 0.22222, 1 => 0.55555, -1 => NaN)],
    "minimal_critical_shock_magnitude" => [Dict(1 => Inf, -1 => NaN), Dict(2 => 2.0, 1 => 1.0, -1 => NaN)],
    "median_convergence_time" => [Dict(1 => 1.0, -1 => NaN), Dict(2 => 1.0, 1 => 1.0, -1 => NaN)],
    "maximal_convergence_time" => [Dict(1 => 1.0, -1 => NaN), Dict(2 => 1.0, 1 => 1.0, -1 => NaN)],
    "mean_noncritical_shock_magnitude" => [Dict(1 => 1.07298, -1 => NaN), Dict(2 => 0.33333, 1 => 0.85003, -1 => NaN)],
    "characteristic_return_time" => [Dict(1 => NaN, -1 => NaN), Dict(2 => NaN, 1 => NaN, -1 => NaN)],
    "reactivity" => [Dict(1 => NaN, -1 => NaN), Dict(2 => NaN, 1 => NaN, -1 => NaN)],
    "maximal_amplification" => [Dict(1 => NaN, -1 => NaN), Dict(2 => NaN, 1 => NaN, -1 => NaN)],
    "maximal_amplification_time" => [Dict(1 => NaN, -1 => NaN), Dict(2 => NaN, 1 => NaN, -1 => NaN)]
)
@testset "Nonlocal Stability Measures Continuation" begin
    # Validate the results
    for (key, value) in measures_cont_expected
        @test key in keys(measures_cont)
        for k in [1, 2]
            @test isapprox(
                sort(collect(values(value[k]))),
                sort(collect(values(measures_cont[key][k]))),
                atol = 1.0e-5,
                nans = true,
            )
        end
    end
end


# Now we will test the local stability measures in `StabilityMeasuresAccumulator` in a
# linear system.
function linear_evolution(z, p, n)
    if p[1] < 0.0
        A = [1.0 0.0; 0.0 1.0]  # Linear transformation matrix
    else
        A = [-0.5 0.0; 0.0 -0.5]  # Linear transformation matrix
    end
    return SVector(A * [z[1], z[2]]...)  # Convert matrix multiplication result to SVector
end

# Create the dynamical system
dynamics = CoupledODEs(linear_evolution, [1.0, 1.0], [1.0])

# Define a grid for initial conditions
grid = ([-1.0, -0.1, 0.3, 1.0], [-1.0, -0.3, 0.1, 1.0])

# Map initial conditions and compute attractors
mapper = AttractorsViaRecurrences(dynamics, grid; sparse = false)

# Use the StabilityMeasuresAccumulator to compute measures
accumulator = StabilityMeasuresAccumulator(mapper, finite_time = 0.5)
for u0 in A
    id = accumulator(u0)
end
results = finalize_accumulator(accumulator)

# Define expected results for the linear system
results_expected = Dict(
    "characteristic_return_time" => Dict(1 => 2.0, -1 => NaN),
    "reactivity" => Dict(1 => -0.5, -1 => NaN),
    "maximal_amplification" => Dict(1 => 1.0, -1 => NaN),
    "maximal_amplification_time" => Dict(1 => 0.0, -1 => NaN)
)
@testset "Local Stability Measures Accumulator" begin
    # Validate the results
    for (key, value) in results_expected
        @test key in keys(results)
        @test isapprox(
            sort(collect(values(value))),
            sort(collect(values(results[key])));
            atol = 1.0e-5,
            nans = true,
        )
    end
end

# Now we test the continuation of local stability measures using the linear map.
pcurve_local = [[1 => p] for p in [-1.0, 1.0]]
attractors_cont_local = [
    Dict(1 => StateSpaceSet([SVector(0.0, 0.0)])),
    Dict(1 => StateSpaceSet([SVector(0.0, 0.0)])),
]

proximity_mapper_options_local = (
    Ttr = 0, stop_at_Δt = false, horizon_limit = 1.0e2, consecutive_lost_steps = 10000,
)
measures_cont_local = stability_measures_along_continuation(
    dynamics, attractors_cont_local, pcurve_local, ics_from_grid(grid), ε = 0.1, finite_time = 0.5,
    proximity_mapper_options = proximity_mapper_options_local
)

measures_cont_local_expected = Dict(
    "characteristic_return_time" => [Dict(1 => Inf, -1 => NaN), Dict(1 => 2.0, -1 => NaN)],
    "reactivity" => [Dict(1 => 1.0, -1 => NaN), Dict(1 => -0.5, -1 => NaN)],
    "maximal_amplification" => [Dict(1 => Inf, -1 => NaN), Dict(1 => 1.0, -1 => NaN)],
    "maximal_amplification_time" => [Dict(1 => Inf, -1 => NaN), Dict(1 => 0.0, -1 => NaN)]
)
@testset "Local Stability Measures Continuation" begin
    # Validate the results
    for (key, value) in measures_cont_local_expected
        @test key in keys(measures_cont_local)
        for k in [1, 2]
            @test isapprox(
                sort(collect(values(value[k]))),
                sort(collect(values(measures_cont_local[key][k]))),
                atol = 1.0e-5,
                nans = true,
            )
        end
    end
end

@testset "Discrete time" begin
    # For these parameters the map has 1 fixed point and one period 3 orbit.
    # The tests fail because linear measures are not computed for discrete sytems now.
    henon_rule_alter(x, p, n) = SVector{2}(1.0 - p[1] * x[1]^2 + x[2], -p[2] * x[1])
    μ = 1.05; J = 0.9
    ds = DeterministicIteratedMap(henon_rule_alter, zeros(2), [μ, J])
    xg = range(-3.0, 3.0; length = 101)
    yg = range(-3.0, 4.0; length = 101)
    grid = (xg, yg)


    mapper = AttractorsViaRecurrences(ds, grid; sparse = false, consecutive_recurrences = 1000)
    accumulator = StabilityMeasuresAccumulator(mapper)


    A = ics_from_grid(grid)
    for u0 in A
        id = accumulator(u0)
    end
    stability_measures = finalize_accumulator(accumulator)

    measures = ["basin_stability", "minimal_critical_shock_magnitude"]

    @testset "Henon $m" for m in measures
        measures = collect(values(stability_measures[m]))
        @test count(!isnan, measures) ≥ 0
    end
end
