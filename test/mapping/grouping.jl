using Attractors, Test

DO_EXTENSIVE_TESTS = get(ENV, "ATTRACTORS_EXTENSIVE_TESTS", "false") == "true"

@testset "0-dist features" begin
    clusterspecs = GroupViaClustering(rescale_features = false)
    features = [SVector(0.1, 0.2) for _ in 1:100]
    labels = group_features(features, clusterspecs)
    @test labels == fill(1, 100)
end

@testset "nearest feature" begin
    templates = Dict(
        10 => SVector(0.0, 0.0),
        20 => SVector(5.0, 5.0),
        30 => SVector(10.0, -2.0),
    )
    features = [
        SVector(0.1, -0.2),
        SVector(4.7, 5.2),
        SVector(9.5, -2.2),
    ]

    @testset "dictionary templates + group_features" begin
        cfg = GroupViaNearestFeature(templates)
        labels = group_features(features, cfg)
        @test labels == [10, 20, 30]
    end

    @testset "kdtree and brute-force agree" begin
        cfg_tree = GroupViaNearestFeature(templates; use_kdtree = true)
        cfg_brute = GroupViaNearestFeature(templates; use_kdtree = false)
        @test group_features(features, cfg_tree) == group_features(features, cfg_brute)
    end

    @testset "max_distance assigns -1 when too far" begin
        cfg = GroupViaNearestFeature(templates; max_distance = 0.6)
        labels = group_features([
            SVector(0.2, 0.2),
            SVector(40.0, 40.0),
        ], cfg)
        @test labels == [10, -1]
    end

    @testset "vector templates use positional labels" begin
        template_vec = [SVector(1.0, 1.0), SVector(-3.0, 2.0)]
        cfg = GroupViaNearestFeature(template_vec; use_kdtree = true)
        labels = group_features([
            SVector(1.1, 1.1),
            SVector(-2.9, 2.1),
        ], cfg)
        @test labels == [1, 2]
    end

end


if DO_EXTENSIVE_TESTS
    # The functionality tested here has been resolved and is only added as a test
    # for future security. It has no need to be tested in every commit.

    using Attractors
    using Statistics, Random, Test

    @testset "Clustering" begin
        ## prepare test (three possible attractors)
        function featurizer(A, t)
            return SVector(maximum(A[:, 1]), maximum(A[:, 2]))
        end
        function cluster_StateSpaceSets(featurizer, t, StateSpaceSets, clusterspecs)
            features = [featurizer(StateSpaceSets[i], t) for i in 1:length(StateSpaceSets)]
            return group_features(features, clusterspecs)
        end
        attractor_pool = [[0 0], [10 10], [20 20]]
        correct_labels = [1, 1, 1, 1, 2, 2, 2, 1, 2, 3, 3, 3, 3, 1]
        a = attractor_pool[correct_labels]
        correct_labels_infinite_threshold = deepcopy(correct_labels)
        a[end] = [50 5]
        correct_labels[end] = -1
        correct_labels_infinite_threshold[end] = 3
        attractors = Dict(1:length(a) .=> StateSpaceSet.(a; warn = false))

        ### silhouettes and real
        @testset "method=$(optimal_radius_method)" for optimal_radius_method in ["silhouettes", "silhouettes_optim", 5.0]
            for silhouette_statistic in [mean, minimum]
                clusterspecs = GroupViaClustering(;
                    num_attempts_radius = 20, silhouette_statistic,
                    optimal_radius_method = optimal_radius_method, min_neighbors = 2, rescale_features = false
                )
                clust_labels = cluster_StateSpaceSets(featurizer, [], attractors, clusterspecs)
                @test clust_labels == correct_labels
            end
        end

        ### knee method
        @testset "method=knee" begin
            correct_labels_knee = [1, 1, 1, 1, 2, 2, 2, 1, 2, 3, 3, 3, 3, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1] #smaller number of features works even worse
            Random.seed!(1)
            a = [attractor_pool[label] + 0.2 * rand(Float64, (1, 2)) for label in correct_labels_knee]
            attractors_knee = Dict(1:length(a) .=> StateSpaceSet.(a; warn = false))
            clusterspecs = GroupViaClustering(
                optimal_radius_method = "knee",
                min_neighbors = 4, rescale_features = false
            )
            clust_labels = cluster_StateSpaceSets(featurizer, [], attractors_knee, clusterspecs)
            # at least check if it finds the same amount of attractors;
            # note this does not work for any value of `min_neighbors`.
            # TODO: Test here is wrong.
            # @test maximum(clust_labels) == maximum(correct_labels)
            @test 1 ≤ maximum(clust_labels) ≤ 3
        end

        @testset "Mmap" begin
            clusterspecs = GroupViaClustering(;
                num_attempts_radius = 20, use_mmap = true,
                min_neighbors = 2, rescale_features = false
            )
            clust_labels = cluster_StateSpaceSets(featurizer, [], attractors, clusterspecs)
            @test clust_labels == correct_labels
        end

    end

end
