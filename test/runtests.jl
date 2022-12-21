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

    @testset "analysis" begin
        testfile("basins/tipping_points_tests.jl")
        testfile("basins/uncertainty_tests.jl")
    end

    @testset "continuation" begin
        testfile("continuation/matching_attractors_tests.jl")
        testfile("continuation/recurrences_continuation_tests.jl")
        testfile("continuation/clustering_continuation_tests.jl")
    end
end
