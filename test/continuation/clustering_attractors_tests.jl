using Test, Attractors
using Attractors.DynamicalSystemsBase
using Random


    ds = Systems.henon(; b = -0.9, a = 1.4)
    psorig = range(0.6, 1.1; length = 20)
    pidx = 1
    xg = yg = range(-2,2, length = 1000)
    grid = (xg,yg)
    sampler, = statespace_sampler(Random.MersenneTwister(1234);
        min_bounds = minimum.(grid), max_bounds = maximum.(grid)
    )

    ps = psorig
    dts_dis(x,y)= dataset_distance(x,y, Hausdorff())
    clusterspecs = Attractors.GroupViaClustering(;clust_distance_metric = dts_dis, optimal_radius_method = 1.)
    mapper = Attractors.AttractorsViaRecurrences(ds, grid; sparse = true)
    continuation = ClusteringAttractorsContinuation(mapper, par_weight = 1., group_config = clusterspecs)
    fs, att = Attractors.basins_fractions_continuation(continuation,  psorig, pidx, sampler; samples_per_parameter = 10000)
    
