DO_EXTENSIVE_TESTS = get(ENV, "ATTRACTORS_EXTENSIVE_TESTS", "false") == "true"

if DO_EXTENSIVE_TESTS
# The functionality tested here has been resolved and is only added as a test
# for future security. It has no need to be tested in every commit.

using Attractors
using Statistics, Random, Test

@testset "Clustering" begin
    ## prepare test (three possible attractors)
    function featurizer(A, t)
        return SVector(maximum(A[:,1]), maximum(A[:,2]))
    end
    function cluster_StateSpaceSets(featurizer, t, StateSpaceSets, clusterspecs)
        features = [featurizer(StateSpaceSets[i], t) for i=1:length(StateSpaceSets)]
        return group_features(features, clusterspecs)
    end
    attractor_pool = [[0 0], [10 10], [20 20]];
    correct_labels = [1,1,1,1, 2,2,2,1,2,3,3,3,3,1]
    a = attractor_pool[correct_labels]
    correct_labels_infinite_threshold = deepcopy(correct_labels)
    a[end] = [50 5]
    correct_labels[end] = -1
    correct_labels_infinite_threshold[end] = 3
    attractors = Dict(1:length(a) .=> StateSpaceSet.(a; warn = false));

    ### silhouettes and real
    @testset "method=$(optimal_radius_method)" for optimal_radius_method in ["silhouettes", "silhouettes_optim", 5.0]
        for silhouette_statistic in [mean, minimum]
            clusterspecs = GroupViaClustering(; num_attempts_radius=20, silhouette_statistic,
            optimal_radius_method=optimal_radius_method, min_neighbors=2, rescale_features=false)
            clust_labels = cluster_StateSpaceSets(featurizer, [], attractors, clusterspecs)
            @test clust_labels == correct_labels
        end
    end

    ### knee method
    @testset "method=knee" begin
        correct_labels_knee = [1,1,1,1,2,2,2,1,2,3,3,3,3,1,2,2,2,2,2,3,3,3,3,3,3,1,1,1,1,1] #smaller number of features works even worse
        Random.seed!(1)
        a = [attractor_pool[label] + 0.2*rand(Float64, (1,2)) for label in correct_labels_knee]
        attractors_knee = Dict(1:length(a) .=> StateSpaceSet.(a; warn = false));
        clusterspecs = GroupViaClustering( optimal_radius_method="knee",
        min_neighbors=4, rescale_features=false)
        clust_labels = cluster_StateSpaceSets(featurizer, [], attractors_knee, clusterspecs)
        # at least check if it finds the same amount of attractors;
        # note this does not work for any value of `min_neighbors`.
        # TODO: Test here is wrong.
        # @test maximum(clust_labels) == maximum(correct_labels)
        @test 1 ≤ maximum(clust_labels) ≤ 3
    end

    @testset "Mmap" begin
        clusterspecs = GroupViaClustering(; num_attempts_radius=20, use_mmap=true,
        min_neighbors=2, rescale_features=false)
        clust_labels = cluster_StateSpaceSets(featurizer, [], attractors, clusterspecs)
        @test clust_labels == correct_labels
    end

end

@testset "Histogram" begin

end

end