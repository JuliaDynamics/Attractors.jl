# Tests for `BayesianUpdateSampler`, the adaptive initial condition sampler of
# `src/continuation/sampler_api.jl`.

using Test, Attractors
using Attractors: n_boxes, resampling_required

# dummy dynamical system 
function dumb_multistable(u, p, n)
    x = u[1]
    w = p[1]
    xnext = if w ≤ 0
        x < 0 ? -1.0 : 1.0
    else
        x < -w ? -1.0 : (x > w ? 1.0 : 0.0)
    end
    return SVector(xnext, 0.0)
end

# Exact basin fractions over the region [-1, 1]².
function dumb_fractions(w)  
    if  w ≤ 0 
         Dict(1 => 0.5, 3 => 0.5) 
    else
        Dict(1 => (1 - w)/2, 2 => w, 3 => (1 - w)/2)
    end
end

const REGION = ((-1.0, 1.0), (-1.0, 1.0))

function recurrences_mapper(w)
    ds = DeterministicIteratedMap(dumb_multistable, [0.1, 0.1], [w, 0.0])
    grid = (range(-2.0, 2.0; length = 201), range(-2.0, 2.0; length = 201))
    return ds, BasinMapRecurrences(ds, grid; sparse = true, show_progress = false)
end

# The parameter sweeps `w` from the bistable regime across zero, where a third fixed
# point is born and the sampler does ask for re-sampling rounds. `history = true` keeps
# the per-parameter `etas`, so the alarms raised along the way can be read back after
# the fact.
const WS = [-0.5, -0.4, -0.3, 0.2, 0.4, 0.4]
const PCURVE = [[1 => w, 2 => 0.0] for w in WS]

_, BMAP = recurrences_mapper(first(WS))
SAMPLER = BayesianUpdateSampler(REGION, 4;
    sparse_n = 40, dense_n = 400, seed = 09092026, history = true)
FRACTIONS_CONT, ATTRACTORS_CONT = global_continuation(
    AttractorSeedContinueMatch(BMAP), PCURVE, SAMPLER; show_progress = false,
)

@testset "continuation test" begin
    @test length(FRACTIONS_CONT) == length(ATTRACTORS_CONT) == length(WS)
    # Make sure all attractors have been found
    for (atts, w) in zip(ATTRACTORS_CONT, WS)
        @test length(atts) == length(dumb_fractions(w))
    end
    # The sampler is left in a usable state: no flag pending, one prior per box, and the
    # next round is the routine sparse one.
    @test !resampling_required(SAMPLER)
    @test length(SAMPLER) == SAMPLER.sparse_n * n_boxes(SAMPLER)
    @test all(!isempty, SAMPLER.alphas)
end

@testset "fractions test" begin
    for (fs, w) in zip(FRACTIONS_CONT, WS)
        @test sum(values(fs)) ≈ 1
        # We have to sort values because we dont't know the
        # ordering of the labels.
        exact = sort!(collect(values(dumb_fractions(w))))
        got = sort!(collect(values(fs)))
        @test length(got) == length(exact)
        @test got ≈ exact atol = 0.02
    end
end

@testset "alarms test" begin
    # get history
    h = sampler_history(SAMPLER)
    @test length(h.etas) == length(h.alphas) == length(WS)
    alarms = [findall(<(0), e) for e in h.etas]

    # dense initialisation, no test is performed.
    @test all(iszero, h.etas[1])
 
    # η = 0 for w < 0: no alarm is possible, false or otherwise.
    @test all(isempty, alarms[2:3])

    # The fixed point is born at w = 0. Only the boxes overlapping its
    # basin |x| < 0.4 can see it.
    overlapping = [i for (i, b) in enumerate(SAMPLER.boxes)
                   if b.mins[1] < 0.4 && b.maxs[1] > -0.4]
    @test length(overlapping) == 8
    @test alarms[4] == overlapping

    # Last parameter is repeated and the alarms shouldn't 
    # trigger
    @test length(alarms[6]) ≤ 2
end
