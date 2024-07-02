# This file also tests `aggregate_attractor_fractions`!

DO_EXTENSIVE_TESTS = get(ENV, "ATTRACTORS_EXTENSIVE_TESTS", "false") == "true"

using Test, Attractors
using Random

@testset "analytic bistable map" begin
    # This is a fake bistable map that has two equilibrium points
    # for r > 0.5. It has predictable fractions.
    function dumb_map(dz, z, p, n)
        x,y = z
        r = p[1]

        if r < 0.5
            dz[1] = dz[2] = 0.0
        else
            if x > 0
                dz[1] = r
                dz[2] = r
            else
                dz[1] = -r
                dz[2] = -r
            end
        end
        return
    end

    r = 1.0
    ds = DeterministicIteratedMap(dumb_map, [0., 0.], [r])
    yg = xg = range(-10., 10, length = 100)
    grid = (xg,yg)
    mapper = AttractorsViaRecurrences(ds, grid; sparse = true, show_progress = false)

    sampler, = statespace_sampler(grid, 1234)

    rrange = range(0, 2; length = 20)
    ridx = 1
    rsc = RecurrencesFindAndMatch(mapper; threshold = 0.3)
    fractions_cont, a = global_continuation(
        rsc, rrange, ridx, sampler;
        show_progress = false, samples_per_parameter = 1000
    )

    for (i, r) in enumerate(rrange)

        fs = fractions_cont[i]
        k = sort!(collect(keys(fs)))
        @test k == sort!(collect(keys(a[i]))) # no -1 ID, so keys must match
        if r < 0.5
            k = sort!(collect(keys(fs)))
            @test length(k) == 1
        else
            k = sort!(collect(keys(fs)))
            @test length(k) == 2
            v = values(fs)
            for f in v
                @test (0.4 < f < 0.6)
            end
        end
        @test sum(values(fs)) ≈ 1
    end

end


@testset "Multistable aggregating map" begin
    # This is a fake multistable map helps at testing the grouping
    # capabilities. We know what it does analytically!
    function dumb_map(dz, z, p, n)
        x,y = z
        r = p[1]

        θ = mod(angle(x + im*y),2pi)
        rr = abs(x+ im*y)

        # This map works as follows: for the parameter r< 0.5 there are
        # 8 attractors close to the origin (radius 0.1)
        # For r > 0.5 there are 8 attractors very close to the origin (radius 0.01
        # and 8 attractors far away and well separated.
        if r < 0.5
            rr = 0.1
        else
            if rr > 1
                rr = 3
            else
                rr = 0.01
            end
        end

        if 0 < θ ≤ π/4
           x = rr*cos(π/8); y = rr*sin(π/8)
        elseif π/4 < θ ≤ π/2
            x = rr*cos(3π/8); y = rr*sin(3π/8)
        elseif π/2 < θ ≤ 3π/4
            x = rr*cos(5π/8); y = rr*sin(5π/8)
        elseif 3π/4 < θ ≤ π
            x = rr*cos(7π/8); y = rr*sin(7π/8)
        elseif π < θ ≤ 5π/4
            x = rr*cos(9π/8); y = rr*sin(9π/8)
        elseif 5π/4 < θ ≤ 6π/4
            x = rr*cos(11π/8); y = rr*sin(11π/8)
        elseif 6π/4 < θ ≤ 7π/4
            x = rr*cos(13π/8); y = rr*sin(13π/8)
        elseif 7π/4 < θ ≤ 8π/4
            x = rr*cos(15π/8); y = rr*sin(15π/8)
        end
        dz[1]= x; dz[2] = y
        return
    end

    function test_fs(fractions_cont, rrange, frac_results)
        # For Grouping there should be one cluster for r < 0.5 and then 9 groups of attractors
        # For matching, all attractors are detected and matched
        # for r < 0.5 there are 4 attractors and for r > 0.5 there are 12
        for (i, r) in enumerate(rrange)
            fs = fractions_cont[i]
            # non-zero keys
            k = sort!([k for k in keys(fs) if fs[k] > 0])
            if r < 0.5
                @test length(k) == frac_results[1]
            else
                @test length(k) == frac_results[2]
            end
            @test sum(values(fs)) ≈ 1
        end
    end

    r = 0.3
    ds = DeterministicIteratedMap(dumb_map, [0.0, 0.0], [r])
    yg = xg = range(-10., 10, length = 100)
    grid = (xg, yg)

    sampler, = statespace_sampler(grid, 1234)

    rrange = range(0, 2; length = 21)
    ridx = 1
    # First, test the normal function of finding attractors
    mapper = AttractorsViaRecurrences(ds, grid; sparse = true, show_progress = false)
    rsc = RecurrencesFindAndMatch(mapper; threshold = 0.1)
    fractions_cont, attractors_cont = global_continuation(
        rsc, rrange, ridx, sampler;
        show_progress = false, samples_per_parameter = 1000,
    )
    test_fs(fractions_cont, rrange, [4, 12])

    # Then, test the aggregation of features via featurizing and histogram
    using Statistics
    featurizer = (x) -> mean(x)

    hconfig = GroupViaHistogram(
        FixedRectangularBinning(range(-4, 4; step = 0.005), 2)
    )

    aggr_fracs, aggr_info = aggregate_attractor_fractions(
        fractions_cont, attractors_cont, featurizer, hconfig
    )
    test_fs(aggr_fracs, rrange, [4, 12])

end

if DO_EXTENSIVE_TESTS

@testset "Henon map" begin
    henon_rule(x, p, n) = SVector{2}(1.0 - p[1]*x[1]^2 + x[2], p[2]*x[1])
    ds = DeterministicIteratedMap(henon_rule, zeros(2), [1.4, 0.3])
    ps = range(1.2, 1.25; length = 101)
    # In these parameters we go from a chaotic attractor to a period 7 orbit at a≈1.2265
    # (you can see this by launching our wonderful `interactive_orbitdiagram` app).
    # So we can use this to test different matching processes
    # (because "by distance" matches the two kind of attractors already)
    # Notice that the length=101 is rather sensitive and depending on it, some
    # much smaller periodic windows exist in the range.
    # (For 101, a period-14 window exists in the second parameter entry)
    acritical = 1.2265

    xg = yg = range(-2.5, 2.5, length = 500)
    pidx = 1
    sampler, = statespace_sampler(HRectangle([-2,-2], [2,2]), 1234)
    distance_function = function (A, B)
        # length of attractors within a factor of 2, then distance is ≤ 1
        return abs(log(2, length(A)) - log(2, length(B)))
    end
    # notice that without this special distance function, even with a
    # really small threshold like 0.2 we still get a "single" attractor
    # throughout the range. Now we get one with period 14, a chaotic,
    # and one with period 7 that spans the second half of the parameter range
    mapper = AttractorsViaRecurrences(ds, (xg, yg); sparse=false,
        consecutive_recurrences = 3000,
        attractor_locate_steps = 3000
    )
    rsc = RecurrencesFindAndMatch(mapper;
        threshold = 0.99, distance = distance_function
    )
    fractions_cont, attractors_cont = global_continuation(
        rsc, ps, pidx, sampler;
        show_progress = false, samples_per_parameter = 100
    )

    for (i, p) in enumerate(ps)
        fs = fractions_cont[i]
        attractors = attractors_cont[i]
        @test sum(values(fs)) ≈ 1
        # Test that keys are the same (-1 doesn't have attractor)
        k = sort!(collect(keys(fs)))
        -1 ∈ k && deleteat!(k, 1)
        attk = sort!(collect(keys(attractors)))
        @test k == attk
    end

    # unique keys
    ukeys = Attractors.unique_keys(attractors_cont)

    # We must have 3 attractors: initial chaotic, period 14 in the middle,
    # chaotic again, and period 7 at the end. ALl of these should be matched to each other.
    # Notice that due to the storage of "ghost" attractors in `RAFM`, the
    # first and second chaotic attractors are mapped to each other.
    @test ukeys == 1:3

    # # Animation of henon attractors
    # animate_attractors_continuation(ds, attractors_cont, fractions_cont, ps, pidx)
end

@testset "non-found attractors" begin
    # This is standard henon map
    henon_rule(x, p, n) = SVector{2}(1.0 - p[1]*x[1]^2 + x[2], p[2]*x[1])
    ds = DeterministicIteratedMap(henon_rule, zeros(2), [1.4, 0.3])
    ps = range(1.2, 1.25; length = 3)
    # This grid is chosen such that no attractors are in there!
    xg = yg = range(-25, -5; length = 500)
    pidx = 1
    sampler, = statespace_sampler(HRectangle([-2,-2], [2,2]), 1234)
    mapper = AttractorsViaRecurrences(ds, (xg, yg); sparse=false)
    rsc = RecurrencesFindAndMatch(mapper)

    fractions_cont, attractors_cont = global_continuation(
        rsc, ps, pidx, sampler;
        show_progress = false, samples_per_parameter = 100
    )
    @test all(i -> isempty(i), attractors_cont)
end

end

