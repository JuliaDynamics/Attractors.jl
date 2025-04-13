DO_EXTENSIVE_TESTS = get(ENV, "ATTRACTORS_EXTENSIVE_TESTS", "false") == "true"

using Test, Attractors
using Random

# use this dumb map, map initial conditiosn and then call finalize.
# you should know analytically the value for all nonlocal stability measures.
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
# TODO: @andreasmorr please justify where the values of nonfatal shock and median
# convergence pace are coming from
results_expected = Dict(
    "characteristic_return_time"       => Dict(2=>NaN, 1=>NaN),
    "maximal_amplification_time"       => Dict(2=>NaN, 1=>NaN),
    "mean_convergence_time"            => Dict(2=>1.0, 1=>1.0),
    "maximal_nonfatal_shock_magnitude" => Dict(2=>2.23607, 1=>2.0),
    "finite_time_basin_stability"      => Dict(2=>0.0, 1=>0.0),
    "median_convergence_pace"          => Dict(2=>0.853553, 1=>1.0),
    "reactivity"                       => Dict(2=>NaN, 1=>NaN),
    "median_convergence_time"          => Dict(2=>1.0, 1=>1.0),
    "minimal_fatal_shock_magnitude"    => Dict(2=>2.0, 1=>1.0),
    "basin_stability"                  => Dict(2=>0.666667, 1=>0.333333),
    "maximal_convergence_pace"         => Dict(2=>Inf, 1=>Inf),
    "maximal_convergence_time"         => Dict(2=>1.0, 1=>1.0),
    "mean_convergence_pace"            => Dict(2=>Inf, 1=>Inf),
    "maximal_amplification"            => Dict(2=>NaN, 1=>NaN),
    "basin_fraction"                   => Dict(2=>0.666667, 1=>0.333333)
)
# Check if the results are as expected
@testset "Stability Measures Accumulator with dumb map" begin
    for (key, value) in results_expected
        @test key in keys(results)
        @test value ≈ results[key] atol=1e-5
    end
end

# TODO: @andreasmorr please add one more test for the **LINEAR** measures by creating a
# LINEAR dynamical system and estimating the measures for its fixed point.

# TODO: @andreasmorr since the `stability_measures_along_continuation` function is
# exported and part of the public API it also needs to be tested here!