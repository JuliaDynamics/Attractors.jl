# This file also tests `aggregate_attractor_fractions`

# DO_EXTENSIVE_TESTS = get(ENV, "ATTRACTORS_EXTENSIVE_TESTS", "false") == "true"
DO_EXTENSIVE_TESTS = true
using Test, Attractors
using Attractors.DynamicalSystemsBase
using Random

# @testset "magnetic pendulum" begin
#     d, α, ω = 0.3, 0.2, 0.5
#     ds = Systems.magnetic_pendulum(; d, α, ω)
#     xg = yg = range(-3, 3; length = 101)
#     ds = projected_integrator(ds, 1:2, [0.0, 0.0])
#     mapper = AttractorsViaRecurrences(ds, (xg, yg); Δt = 1.0)
#     rr = range(1, 0; length = 101)
#     psorig = [[1, 1, γ] for γ in rr]
#     pidx = :γs
#     # important to make a sampler that respects the symmetry of the system
#     sampler, isinside = statespace_sampler(Xoshiro(1234); spheredims = 2, radius = 3.0)
#     for (j, ps) in enumerate((psorig, reverse(psorig)))
#         # test that both finding and removing attractor works
#         mapper = AttractorsViaRecurrences(ds, (xg, yg); sparse=false, Δt = 1.0)

#         continuation = RecurrencesSeedingContinuation(mapper; threshold = Inf)
#         # With this threshold all attractors are mapped to each other, they are within
#         # distance 1 in state space.
#         fractions_curves, attractors_info = basins_fractions_continuation(
#             continuation, ps, pidx, sampler; show_progress = false, samples_per_parameter = 1000
#         )

#         # Keys of the two attractors that always exist
#         twokeys = collect(keys(fractions_curves[(j == 2 ? 1 : 101)]))

#         @testset "symmetry respect" begin
#             # Initially fractions are all 0.33 but at the end only two of 0.5 remain
#             # because only two attractors remain (with equal magnetic pull)
#             startfracs, endfracs = j == 1 ? [0.33, 0.5] : [0.5, 0.33]
#             @test all(v -> isapprox(v, startfracs; atol = 1e-1), values(fractions_curves[1]))
#             @test all(v -> isapprox(v, endfracs; atol = 1e-1), values(fractions_curves[end]))
#         end

#         for (i, p) in enumerate(ps)
#             γ = p[3]
#             fs = fractions_curves[i]
#             attractors = attractors_info[i]
#             k = sort!(collect(keys(fs)))
#             @test maximum(k) ≤ 3
#             attk = sort!(collect(keys(attractors)))
#             @test k == attk
#             @test all(fk -> fk ∈ k, twokeys)

#             # It is arbitrary what id we get, because the third
#             # fixed point that vanishes could have any of the three ids
#             # But we can test for sure how many ids we have
#             # (depending on where we come from we find the attractor for longer)
#             if γ < 0.2
#                 @test length(k) == 2
#             elseif γ > 0.24
#                 @test length(k) == 3
#             else
#                 # There is a bit of varaibility of exactly when the transition
#                 # occurs, and also depends on randomness for when we get exactly 0
#                 # fraction for one of the attractors
#                 @test length(k) ∈ (2, 3)
#             end
#             @test sum(values(fs)) ≈ 1
#         end
#         # # Plot code for fractions
#         # using GLMakie
#         # x = [fs[finalkeys[1]] for fs in fractions_curves]
#         # y = [fs[finalkeys[2]] for fs in fractions_curves]
#         # z = zeros(length(x))
#         # fig = Figure(resolution = (400, 300))
#         # ax = Axis(fig[1,1])
#         # display(fig)
#         # γs = [p[3] for p in ps]
#         # band!(ax, γs, z, x; color = Cycled(1), label = "1")
#         # band!(ax, γs, x, x .+ y; color = Cycled(2), label  = "2")
#         # band!(ax, γs, x .+ y, 1; color = Cycled(3), label = "3")
#         # xlims!(ax, 0, 1)
#         # ylims!(ax, 0, 1)
#         # ax.ylabel = "fractions"
#         # ax.xlabel = "magnet strength"
#         # axislegend(ax)
#         # Makie.save("magnetic_fracs.png", fig; px_per_unit = 4)
#     end
# end

if DO_EXTENSIVE_TESTS

@testset "Henon map" begin
    ds = Systems.henon(; b = 0.3, a = 1.4)
    psorig = range(1.2, 1.25; length = 101)
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
    sampler, = statespace_sampler(Random.MersenneTwister(1234);
        min_bounds = [-2,-2], max_bounds = [2,2]
    )
    distance_function = function (A, B)
        # length of attractors within a factor of 2, then distance is ≤ 1
        return abs(log(2, length(A)) - log(2, length(B)))
    end
    # notice that without this special distance function, even with a
    # really small threshold like 0.2 we still get a "single" attractor
    # throughout the range. Now we get one with period 14, a chaotic,
    # and one with period 7 that spans the second half of the parameter range
    mapper = AttractorsViaRecurrences(ds, (xg, yg); sparse=false,
        mx_chk_fnd_att = 3000,
        mx_chk_loc_att = 3000
    )
    continuation = RecurrencesSeedingContinuation(mapper;
        threshold = 0.99, method = distance_function
    )
    ps = psorig
    fractions_curves, attractors_info = basins_fractions_continuation(
        continuation, ps, pidx, sampler;
        show_progress = false, samples_per_parameter = 100
    )

    for (i, p) in enumerate(ps)
        fs = fractions_curves[i]
        attractors = attractors_info[i]
        @test sum(values(fs)) ≈ 1
        # Test that keys are the same (-1 doesn't have attractor)
        k = sort!(collect(keys(fs)))
        -1 ∈ k && deleteat!(k, 1)
        attk = sort!(collect(keys(attractors)))
        @test k == attk
    end

    # unique keys
    ukeys = Attractors.unique_keys(attractors_info)
    # We must have 4 attractors: initial chaotic, period 14 in the middle,
    # chaotic again, and period 7 at the end. ALl of these should be matched to each other.
    # Since we retract keys, we have 1:4
    @test ukeys == 1:4

    # # Animation of henon attractors
    # using GLMakie
    # fig = Figure(); display(fig)
    # ax = Axis(fig[1,1]; limits = (-2,2,-1,1))
    # colors = Dict(k => Cycled(i) for (i, k) in enumerate(ukeys))
    # att_obs = Dict(k => Observable(Point2f[]) for k in ukeys)
    # for k in ukeys
    #     scatter!(ax, att_obs[k]; color = colors[k],
    #     label = "$k", markersize = 8)
    # end
    # axislegend(ax)
    # display(fig)
    # record(fig, "henon_test.mp4", eachindex(ps); framerate = 5) do i
    #     p = ps[i]
    #     ax.title = "p = $p"
    #     # fs = fractions_curves[i]
    #     attractors = attractors_info[i]
    #     set_parameter!(ds, pidx, p)
    #     for (k, att) in attractors
    #         tr = trajectory(ds, 1000, att[1]; Δt = 1)
    #         att_obs[k][] = vec(tr)
    #         notify(att_obs[k])
    #     end
    #     # also ensure that attractors that don't exist are cleared
    #     for k in setdiff(ukeys, collect(keys(attractors)))
    #         att_obs[k][] = Point2f[]; notify(att_obs[k])
    #     end
    # end
end

@testset "non-found attractors" begin
    # This is standard henon map
    ds = Systems.henon(; b = 0.3, a = 1.4)
    ps = range(1.2, 1.25; length = 3)
    # This grid is chosen such that no attractors are in there!
    xg = yg = range(-25, -5; length = 500)
    pidx = 1
    sampler, = statespace_sampler(Random.MersenneTwister(1234);
        min_bounds = [-2,-2], max_bounds = [2,2]
    )
    mapper = AttractorsViaRecurrences(ds, (xg, yg); sparse=false)
    continuation = RecurrencesSeedingContinuation(mapper)
    fractions_curves, attractors_info = basins_fractions_continuation(
        continuation, ps, pidx, sampler;
        show_progress = false, samples_per_parameter = 100
    )
    @test all(i -> isempty(i), attractors_info)
end


@testset "dumb bistable map" begin
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


    r = 1.
    ds = DiscreteDynamicalSystem(dumb_map, [0., 0.], [r])
    yg = xg = range(-10., 10, length = 100)
    grid = (xg,yg)
    mapper = Attractors.AttractorsViaRecurrences(ds, grid; sparse = true, show_progress = false)

    sampler, = statespace_sampler(Random.MersenneTwister(1234);
        min_bounds = minimum.(grid), max_bounds = maximum.(grid))

    rrange = range(0., 2; length = 20)
    ridx = 1
    continuation = Attractors.RecurrencesSeedingContinuation(mapper; threshold = 0.3)
    fractions_curves, a = Attractors.basins_fractions_continuation(
        continuation, rrange, ridx, sampler;
        show_progress = false, samples_per_parameter = 1000
    )

    # fractions_curves, a = Attractors.basins_fractions_continuation_group(
    #     continuation, rrange, ridx, sampler;
    #     show_progress = false, samples_per_parameter = 1000
    # )

    for (i, r) in enumerate(rrange)

        fs = fractions_curves[i]
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

end
