using Test, Attractors

DO_EXTENSIVE_TESTS = get(ENV, "ATTRACTORS_EXTENSIVE_TESTS", "false") == "true"

@testset "matching utils" begin
    @testset "No double duplication" begin
        rmap = Dict(4 => 3, 3 => 2)
        a = Dict(4 => 1)
        swap_dict_keys!(a, rmap)
        @test haskey(a, 3)
        @test !haskey(a, 2)
    end
end

@testset "MatchBySSSetDistance" begin
    default = MatchBySSSetDistance()

@testset "analytic" begin
    a_befo = Dict(1 => [SVector(0.0, 0.0)], 2 => [SVector(1.0, 1.0)])
    a_befo = Dict(keys(a_befo) .=> StateSpaceSet.(values(a_befo)))
    @testset "infinite threshold" begin
        a_afte = Dict(2 => [SVector(0.0, 0.0)], 1 => [SVector(2.0, 2.0)])
        a_afte = Dict(keys(a_afte) .=> StateSpaceSet.(values(a_afte)))
        rmap = matching_map!(a_afte, a_befo, default)
        @test rmap == Dict(1 => 2, 2 => 1)
        @test a_afte[1] == a_befo[1] == StateSpaceSet([SVector(0.0, 0.0)])
        @test haskey(a_afte, 2)
        @test a_afte[2] == StateSpaceSet([SVector(2.0, 2.0)])
    end
    @testset "separating threshold" begin
        a_afte = Dict(2 => [SVector(0.0, 0.0)], 1 => [SVector(2.0, 2.0)])
        a_afte = Dict(keys(a_afte) .=> StateSpaceSet.(values(a_afte)))
        matcher = MatchBySSSetDistance(threshold = 0.1)
        rmap = matching_map!(a_afte, a_befo, matcher)
        @test rmap == Dict(1 => 3, 2 => 1)
        @test a_afte[1] == a_befo[1] == StateSpaceSet([SVector(0.0, 0.0)])
        @test !haskey(a_afte, 2)
        @test a_afte[3] == StateSpaceSet([SVector(2.0, 2.0)])
    end
end

@testset "synthetic multistable matching" begin
    attractors_cont_simple = Dict{Int64, SVector{1, Float64}}[Dict(1 => [0.0]), Dict(1 => [0.0]), Dict(2 => [2.0], 1 => [0.0]), Dict(2 => [0.0], 1 => [2.0]), Dict(2 => [2.0], 3 => [4.0], 1 => [-5.0]), Dict(2 => [4.0], 3 => [-5.0], 1 => [2.0]), Dict(4 => [6.0], 2 => [0.0], 3 => [2.0], 1 => [4.0]), Dict(4 => [4.0], 2 => [0.0], 3 => [2.0], 1 => [6.0]), Dict(5 => [8.0], 4 => [6.0], 2 => [0.0], 3 => [2.0], 1 => [4.0])]
    attractors_cont = [Dict(k=>StateSpaceSet(Vector(v)) for (k,v) in atts) for atts in attractors_cont_simple]

    mapped_atts = deepcopy(attractors_cont)
    rmaps = match_sequentially!(mapped_atts, default)

    fractions_cont = [Dict(1 => 1.0), Dict(1 => 1.0), Dict(2 => 0.8, 1 => 0.2), Dict(2 => 0.2, 1 => 0.8), Dict(2 => 0.2, 3 => 0.6, 1 => 0.2), Dict(2 => 0.6, 3 => 0.2, 1 => 0.2), Dict(4 => 0.4, 2 => 0.2, 3 => 0.2, 1 => 0.2), Dict(4 => 0.2, 2 => 0.2, 3 => 0.2, 1 => 0.4), Dict(5 => 0.2, 4 => 0.2, 2 => 0.2, 3 => 0.2, 1 => 0.2)]
    mapped_fracs = deepcopy(fractions_cont)
    match_sequentially!(mapped_fracs, rmaps)

    @test all(keys.(attractors_cont) .== keys.(mapped_atts)    )
    @test all(keys.(attractors_cont) .== keys.(mapped_fracs)    )
    @test all(Set.(values.(fractions_cont)) .== Set.(values.(mapped_fracs)))
end

@testset "global_continuation matching" begin
    # Make fake attractors with points that become more "separated" as "parameter"
    # is increased
    jrange = 0.1:0.1:1
    allatts = [Dict(1 => [SVector(0.0, 0.0)], 2 => [SVector(j, j)]) for j in jrange]
    allatts = [Dict(keys(d) .=> StateSpaceSet.(values(d))) for d in allatts]
    for i in eachindex(jrange)
        if isodd(i) && i ≠ 1
            # swap key of first attractor to from 1 to i
            allatts[i][i] = allatts[i][1]
            delete!(allatts[i], 1)
        end
    end
    # Test with distance not enough to increment
    match_sequentially!(allatts, MatchBySSSetDistance(threshold = 100.0)) # all odd keys become 1
    @test all(haskey(d, 1) for d in allatts)
    @test all(haskey(d, 2) for d in allatts)
    # Test with distance enough to increment
    allatts2 = deepcopy(allatts)
    match_sequentially!(allatts2, MatchBySSSetDistance(threshold = 0.1)) # all keys there were `2` get incremented
    @test all(haskey(d, 1) for d in allatts2)
    for i in 2:length(jrange)
        @test haskey(allatts2[i], i+1)
        @test !haskey(allatts2[i], 2)
    end
    @test haskey(allatts2[1], 2)
end

@testset "global_continuation matching advanced" begin
    jrange = 0.1:0.1:1
    allatts = [Dict(1 => [SVector(0.0, 0.0)], 2 => [SVector(1.0, 1.0)], 3 => [SVector((10j)^2, 0)]) for j in jrange]
    allatts = [Dict(keys(d) .=> StateSpaceSet.(values(d))) for d in allatts]
    # delete attractors every other parameter
    for i in eachindex(jrange)
        if iseven(i)
            delete!(allatts[i], 3)
        end
    end

    @testset "ignore vanished" begin
        @testset "no retract" begin
            atts = deepcopy(allatts)
            rmaps = match_sequentially!(atts, default; retract_keys = false)
            # After the first 3 key, all subsequent keys 3 become the next integer,
            # and since we started cutting away keys 3 from `i = 2` we have
            # 4 extra 3 keys to add.
            @test unique_keys(atts) == [1, 2, 3, 5, 7, 9, 11]
        end
        @testset "with retract" begin
            # with retraction it becomes nice
            atts = deepcopy(allatts)
            rmaps = match_sequentially!(atts, default; retract_keys = true)
            @test unique_keys(atts) == 1:7
        end
    end

    @testset "use vanished" begin
        @testset "Inf thresh" begin
            atts = deepcopy(allatts)
            match_sequentially!(atts, MatchBySSSetDistance(use_vanished = true))
            @test unique_keys(atts) == 1:3
            for i in eachindex(jrange)
                if iseven(i)
                    @test sort!(collect(keys(atts[i]))) == 1:2
                else
                    @test sort!(collect(keys(atts[i]))) == 1:3
                end
            end
        end
        @testset "finite thresh" begin
            # okay here we test the case that the threshold becomes too large
            threshold = 10.0 # at the 5th index, we cannot match anymore
            atts = deepcopy(allatts)
            match_sequentially!(atts, MatchBySSSetDistance(; use_vanished = true, threshold))
            @testset "i=$(i)" for i in eachindex(jrange)
                if iseven(i)
                    @test sort!(collect(keys(atts[i]))) == 1:2
                else
                    if i < 5
                        @test sort!(collect(keys(atts[i]))) == 1:3
                    else
                        @test sort!(collect(keys(atts[i]))) == [1, 2, i÷2 + 2]
                    end
                end
            end
        end
        @testset "no attractors for a step" begin
            atts = deepcopy(allatts)
            atts[5] = Dict()
            match_sequentially!(atts, MatchBySSSetDistance(use_vanished = true))
            @test unique_keys(atts) == 1:3
            for i in eachindex(jrange)
                if iseven(i)
                    @test sort!(collect(keys(atts[i]))) == 1:2
                else
                    if i == 5
                        @test isempty(atts[i])
                    else
                        @test sort!(collect(keys(atts[i]))) == 1:3
                    end
                end
            end
        end
    end
end

end # Matcher by distance tests


@testset "basin overlap" begin
    b1 = ones(Int, 10, 10)
    b1[:, 6:10] .= 2
    b2 = copy(b1)
    b2[:, 1:2] .= 2
    b2[:, 3:10] .= 3 # 3 has more overlap to 2, while 2 has more overlap to 1
    b3 = copy(b2)
    b3[:, 1:10] .= 3 # so now this would be teh same as the previous 3, which is will become 2

    m1 = MatchByBasinOverlap()

    rmap = matching_map(b2, b1, m1)
    @test rmap[3] == 2
    @test rmap[2] == 1

    bx = copy(b2)
    matching_map!(bx, b1, m1)
    @test sort(unique(bx)) == [1, 2]

    rmap = matching_map(b3, bx, m1)
    @test length(rmap) == 1
    @test rmap[3] == 2

    # test that there won't be matching with threshold
    m2 = MatchByBasinOverlap(1.99)
    rmap = matching_map(b2, b1, m2)
    # here 3 should go to 2, as it overlaps all of 2
    # while 2 cannot go to 1, as it doesn't 50% or more of 1.
    # The next available ID is not 4 though, it is 3, as this is
    # the next available integer given the previous keys (1, 2)
    @test length(rmap) == 2
    @test rmap[2] == 3
    @test rmap[3] == 2
end

@testset "BasinEncloure" begin

    @testset "synthetic multistable continuation" begin
        function dummy_multistable_equilibrium!(dx, x, p, n)
            r = p[1]
            num_max_atts = 5
            x_max_right = 10
            x_pos_atts = [x_max_right*(i-1)/num_max_atts for i=1:num_max_atts]
            if 3 <=r < 4 && (x[1] <= 2)
                x_pos_atts[1] = -5
            end
            num_atts = r < (num_max_atts + 1) ? floor(Int, r) : ceil(Int, 2*num_max_atts - r)
            x_atts = [x_pos_atts[i] for i=1:num_atts]
            att_of_x = findlast(xatt->xatt<=x[1], x_atts)
            if x[1] < 0 att_of_x = 1 end
            dx .= x_atts[att_of_x]

            return nothing
        end

        ds = DeterministicIteratedMap(dummy_multistable_equilibrium!, [0.], [1.0])
        featurizer(A,t) = A[end]
        grouping_config = GroupViaPairwiseComparison(; threshold=0.2)
        mapper = AttractorsViaFeaturizing(ds, featurizer, grouping_config)

        xg = range(0, 10, length = 100)
        grid = (xg,)
        sampler, = statespace_sampler(grid, 1234)
        samples_per_parameter = 1000
        ics = StateSpaceSet([deepcopy(sampler()) for _ in 1:samples_per_parameter])

        rrange = range(1, 9.5; step=0.5)
        ridx = 1

        mapper = AttractorsViaFeaturizing(ds, featurizer, grouping_config; T=10, Ttr=1)
        matcher  = MatchByBasinEnclosure(;ε=0.1)
        assc = AttractorSeedContinueMatch(mapper, matcher)
        fs_curves, atts_all = global_continuation(assc, rrange, ridx, ics; show_progress = true)

        atts_all_endpoint_solution = Dict{Int64, SVector{1, Float64}}[Dict(1 => [0.0]), Dict(1 => [0.0]), Dict(2 => [2.0], 1 => [0.0]), Dict(2 => [2.0], 1 => [0.0]), Dict(2 => [2.0], 3 => [4.0], 1 => [-5.0]), Dict(2 => [2.0], 3 => [4.0], 1 => [-5.0]), Dict(4 => [6.0], 2 => [2.0], 3 => [4.0], 1 => [0.0]), Dict(4 => [6.0], 2 => [2.0], 3 => [4.0], 1 => [0.0]), Dict(5 => [8.0], 4 => [6.0], 2 => [2.0], 3 => [4.0], 1 => [0.0]), Dict(5 => [8.0], 4 => [6.0], 2 => [2.0], 3 => [4.0], 1 => [0.0]), Dict(4 => [6.0], 2 => [2.0], 3 => [4.0], 1 => [0.0]), Dict(4 => [6.0], 2 => [2.0], 3 => [4.0], 1 => [0.0]), Dict(2 => [2.0], 3 => [4.0], 1 => [0.0]), Dict(2 => [2.0], 3 => [4.0], 1 => [0.0]), Dict(2 => [2.0], 1 => [0.0]), Dict(2 => [2.0], 1 => [0.0]), Dict(1 => [0.0]), Dict(1 => [0.0])]
        atts_all_endpoint = [Dict(k=>v[end] for (k,v) in atts) for atts in atts_all]

        @test atts_all_endpoint == atts_all_endpoint_solution

        fs_curves_solution = [Dict(1 => 1.0), Dict(1 => 1.0), Dict(2 => 0.8, 1 => 0.2), Dict(2 => 0.8, 1 => 0.2), Dict(2 => 0.2, 3 => 0.6, 1 => 0.2), Dict(2 => 0.2, 3 => 0.6, 1 => 0.2), Dict(4 => 0.4, 2 => 0.2, 3 => 0.2, 1 => 0.2), Dict(4 => 0.4, 2 => 0.2, 3 => 0.2, 1 => 0.2), Dict(5 => 0.2, 4 => 0.2, 2 => 0.2, 3 => 0.2, 1 => 0.2), Dict(5 => 0.2, 4 => 0.2, 2 => 0.2, 3 => 0.2, 1 => 0.2), Dict(4 => 0.4, 2 => 0.2, 3 => 0.2, 1 => 0.2), Dict(4 => 0.4, 2 => 0.2, 3 => 0.2, 1 => 0.2), Dict(2 => 0.2, 3 => 0.6, 1 => 0.2), Dict(2 => 0.2, 3 => 0.6, 1 => 0.2), Dict(2 => 0.8, 1 => 0.2), Dict(2 => 0.8, 1 => 0.2), Dict(1 => 1.0), Dict(1 => 1.0)]
        @test all(keys.(fs_curves_solution) .== keys.(fs_curves)    )
        for (fs_curve, fs_curve_solution) in zip(fs_curves, fs_curves_solution)
            for (k, fs) in fs_curve
                @test isapprox(fs, fs_curve_solution[k], atol=1e-1)
            end
        end
    end

end
