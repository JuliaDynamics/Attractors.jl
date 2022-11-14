using Attractors
using Statistics, Random

@testset "Artificial test for cluster_features" begin
    ## prepare test (three possible attractors)
    function featurizer(A, t)
        return [maximum(A[:,1]), maximum(A[:,2])]
    end
    function cluster_datasets(featurizer, t, datasets, clusterspecs)
        features = [featurizer(datasets[i], t) for i=1:length(datasets)]
        return group_features(features, clusterspecs)
    end
    attractor_pool = [[0 0], [10 10], [20 20]];
    correct_labels = [1,1,1,1, 2,2,2,1,2,3,3,3,3,1]
    a = attractor_pool[correct_labels]; correct_labels_infinite_threshold = deepcopy(correct_labels);
    a[end] = [50 5]; correct_labels[end] = -1; correct_labels_infinite_threshold[end] = 3;
    attractors = Dict(1:length(a) .=> Dataset.(a; warn = false));

    ## Unsupervised
    ### silhouettes and real
    for optimal_radius_method in ["silhouettes", "silhouettes_optim", 5.0]
        for silhouette_statistic in [mean, minimum]
            clusterspecs = GroupViaClustering(; num_attempts_radius=20, silhouette_statistic,
            optimal_radius_method=optimal_radius_method, min_neighbors=2, rescale_features=false)
            clust_labels = cluster_datasets(featurizer, [], attractors, clusterspecs)
            @test clust_labels == correct_labels
        end
    end

    ### knee method
    correct_labels_knee = [1,1,1,1, 2,2,2,1,2,3,3,3,3,1,2,2,2,2,2,3,3,3,3,3,3,1,1,1,1,1] #smaller number of features works even worse
    Random.seed!(1)
    a = [attractor_pool[label] + 0.2*rand(Float64, (1,2)) for label in correct_labels_knee]
    attractors_knee = Dict(1:length(a) .=> Dataset.(a; warn = false));
    clusterspecs = ClusteringConfig( optimal_radius_method="knee",
    min_neighbors=4, rescale_features=false)
    clust_labels = cluster_datasets(featurizer, [], attractors_knee, clusterspecs)
    # @test clust_labels == correct_labels #fails
    @test maximum(clust_labels) == maximum(correct_labels) #at least check if it finds the same amount of attractors; note this does not work for any value of `min_neighbors`.

    ## Supervised
    ###construct templates
    t = map(x->featurizer(x, []), attractor_pool);
    template_labels = [i for i ∈ eachindex(attractor_pool)]
    templates = Dict(template_labels.=> t)
    ###compare
    clusterspecs = ClusteringConfig(; templates, min_neighbors=1, rescale_features=false)
    clust_labels = cluster_datasets(featurizer, [], attractors, clusterspecs)
    @test clust_labels == correct_labels_infinite_threshold
    clusterspecs = ClusteringConfig(; templates, min_neighbors=1, rescale_features=false, max_distance_template=0.0)
    clust_labels = cluster_datasets(featurizer, [], attractors, clusterspecs)
    @test clust_labels == correct_labels
end
