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

grid = ([-1, 0, 1], [-1, 0, 1],)

mapper = AttractorsViaRecurrences(dynamics, grid; sparse = false)

A = ics_from_grid(grid)
for u0 in A
    id = mapper(u0)
end

attractors = extract_attractors(mapper)

mapper = AttractorsViaProximity(dynamics, attractors, 0.01, Ttr=0)

accumulator = StabilityMeasuresAccumulator(mapper, finite_time=0.5)

for u0 in A
    id = accumulator(u0)
end

results = finalize_accumulator(accumulator)

# The expected results are:
# maximal_nonfatal_shock of attractor 2 at [1, 1] is computed as the distance between the
# that point and the furthest point of its own basin at [0, -1] which is 2.23607. For
# attractor 1 at [-1, -1] the maximal nonfatal shock is 2.0, which is the distance to the
# point [-1, 1] which is the furthest point in its basin.
# Similarly, the median convergence pace is the constant convergence time 1.0 divided by the
# distance to the point of median distance in the basin. For attractor 2, there are an even
# number of basin points, so the median is the mean of the two middle values. Here,
# 1.0 / 1.41421 = 0.707106 and 1.0 / 1.0, so 0.85355 for the points at [0,0] and [1, 0] for 
# example. 
results_expected = Dict(
    "mean_convergence_time"            => Dict(2=>1.0, 1=>1.0),
    "maximal_nonfatal_shock_magnitude" => Dict(2=>2.23607, 1=>2.0),
    "finite_time_basin_stability"      => Dict(2=>0.0, 1=>0.0),
    "median_convergence_pace"          => Dict(2=>0.85355, 1=>1.0),
    "median_convergence_time"          => Dict(2=>1.0, 1=>1.0),
    "minimal_fatal_shock_magnitude"    => Dict(2=>2.0, 1=>1.0),
    "basin_stability"                  => Dict(2=>0.66667, 1=>0.33333),
    "maximal_convergence_pace"         => Dict(2=>Inf, 1=>Inf),
    "maximal_convergence_time"         => Dict(2=>1.0, 1=>1.0),
    "mean_convergence_pace"            => Dict(2=>Inf, 1=>Inf),
    "basin_fraction"                   => Dict(2=>0.66667, 1=>0.33333)
)
# Check if the results are as expected
@testset "Nonlocal Stability Measures Accumulator with dumb map" begin
    for (key, value) in results_expected
        @test key in keys(results)
        @test value ≈ results[key] atol=1e-5
    end
end

# Now we test the continuation of nonlocal stability measures.
pcurve = [[1 => p] for p in [-1.0, 1.0]]
acam = AttractorSeedContinueMatch(accumulator, MatchBySSSetDistance(), seeding=A->[])
measures_cont, attractors_cont = global_continuation(acam, pcurve, ics_from_grid(grid))

proximity_mapper_options = (
    Ttr=0, stop_at_Δt = false, horizon_limit = 1e2, consecutive_lost_steps = 10000
)
measures_cont = stability_measures_along_continuation(
    dynamics, attractors_cont, pcurve, ics_from_grid(grid), ε=0.1, finite_time=0.5, 
    proximity_mapper_options = proximity_mapper_options
)

measures_cont_expected = Dict(
    "finite_time_basin_stability"      => [Dict(1=>0.0), Dict(2=>0.0, 1=>0.0)],
    "maximal_nonfatal_shock_magnitude" => [Dict(1=>Inf), Dict(2=>2.0, 1=>2.23607)],
    "median_convergence_pace"          => [Dict(1=>1.0), Dict(2=>1.0, 1=>0.85355)],
    "basin_stability"                  => [Dict(1=>1.0), Dict(2=>0.33333, 1=>0.66667)],
    "maximal_convergence_pace"         => [Dict(1=>Inf), Dict(2=>Inf, 1=>Inf)],
    "mean_convergence_pace"            => [Dict(1=>Inf), Dict(2=>Inf, 1=>Inf)],
    "basin_fraction"                   => [Dict(1=>1.0), Dict(2=>0.33333, 1=>0.66667)],
    "mean_convergence_time"            => [Dict(1=>1.0), Dict(2=>1.0, 1=>1.0)],
    "minimal_fatal_shock_magnitude"    => [Dict(1=>Inf), Dict(2=>1.0, 1=>2.0)],
    "median_convergence_time"          => [Dict(1=>1.0), Dict(2=>1.0, 1=>1.0)],
    "maximal_convergence_time"         => [Dict(1=>1.0), Dict(2=>1.0, 1=>1.0)],
)
@testset "Nonlocal Stability Measures Continuation" begin
    # Validate the results
    for (key, value) in measures_cont_expected
        @test key in keys(measures_cont)
        for k in [1, 2] 
            for (i, v) in value[k]
                @test i in keys(measures_cont[key][k])
                @test v ≈ measures_cont[key][k][i] atol=1e-5
            end
        end
    end
end


# Now we will test the local stability measures in `StabilityMeasuresAccumulator` in a 
# linear system.
function linear_evolution(z, p, n)
    A = [-0.5 0.0; 0.0 -0.5]  # Linear transformation matrix
    return SVector(A * [z[1], z[2]]...)  # Convert matrix multiplication result to SVector
end

# Create the dynamical system
dynamics = CoupledODEs(linear_evolution, [1.0, 1.0], [0.0])

# Define a grid for initial conditions
grid = ([-1.0, -0.1, 0.3, 1.0], [-1.0, -0.3, 0.1, 1.0])

# Map initial conditions and compute attractors
mapper = AttractorsViaRecurrences(dynamics, grid; sparse = false)

# Use the StabilityMeasuresAccumulator to compute measures
accumulator = StabilityMeasuresAccumulator(mapper, finite_time=0.5)
for u0 in A
    id = accumulator(u0)
end
results = finalize_accumulator(accumulator)

# Define expected results for the linear system
results_expected = Dict(
    "characteristic_return_time" => Dict(1 => 2.0),
    "reactivity"       => Dict(1 => -0.5),
    "maximal_amplification" => Dict(1 => 1.0),
    "maximal_amplification_time" => Dict(1 => 0.0)
)
@testset "Local Stability Measures Accumulator" begin
    # Validate the results
    for (key, value) in results_expected
        @test key in keys(results)
        @test value ≈ results[key] atol=1e-5
    end
end