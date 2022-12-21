# Additional tests may be run in this test suite according to an environment variable
# `ATTRACTORS_EXTENSIVE_TESTS` which can be true or false.
# If false, a small, but representative subset of tests is used.

# ENV["ATTRACTORS_EXTENSIVE_TESTS"] = true or false (set externally)

using Test

defaultname(file) = uppercasefirst(replace(splitext(basename(file))[1], '_' => ' '))
testfile(file, testname=defaultname(file)) = @testset "$testname" begin; include(file); end

@testset "Attractors.jl" begin
    @testset "mapping" begin
        testfile("mapping/grouping.jl")
        testfile("mapping/recurrence.jl")
        testfile("mapping/proximity_deduce_ε.jl")
        testfile("mapping/attractor_mapping.jl")
    end

    @testset "basins analysis" begin
        testfile("basins/tipping_points_tests.jl")
        testfile("basins/uncertainty_tests.jl")
    end

    @testset "continuation" begin
        testfile("continuation/matching_attractors.jl")
        testfile("continuation/recurrences_continuation.jl")
        testfile("continuation/grouping_continuation.jl")
    end
end
