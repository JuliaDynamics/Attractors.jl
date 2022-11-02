using Test

defaultname(file) = uppercasefirst(replace(splitext(basename(file))[1], '_' => ' '))
testfile(file, testname=defaultname(file)) = @testset "$testname" begin; include(file); end

@testset "Attractors.jl" begin
    testfile("mapping/clustering_tests.jl")
    testfile("mapping/recurrence_tests.jl")
    testfile("mapping/proximity_deduce_ε_tests.jl")
    testfile("mapping/attractor_mapping_tests.jl")

    testfile("basins/tipping_points_tests.jl")
    testfile("basins/uncertainty_tests.jl")

    testfile("continuation/matching_attractors_tests.jl")
    testfile("continuation/recurrences_continuation_tests.jl")
end
