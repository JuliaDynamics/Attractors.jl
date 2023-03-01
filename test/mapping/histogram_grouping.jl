using Test
using Attractors
import Random
using Statistics: mean

ds = Systems.henon(zeros(2); a = 1.4, b = 0.3)
xg = yg = range(-2.0, 2.0; length=100)
grid = (xg, yg)
expected_fs_raw = Dict(1 => 0.451, -1 => 0.549)
function featurizer(A, t)
    # Notice that unsupervised clustering cannot support "divergence to infinity",
    # which it identifies as another attractor (in fact, the first one).
    x = SVector(mean(A[:, 1]), mean(A[:, 2]))
    if any(x -> (isinf(x) || isnan(x) || x > 10), x)
        SVector(200.0, 200.0)
    else
        x
    end
end

binning = FixedRectangularBinning(range(-10, 240; step = 50), 2)
gconfig = GroupViaHistogram(binning)

mapper = AttractorsViaFeaturizing(ds, featurizer, gconfig)

xg = yg = range(-2.0, 2.0; length=100)
grid = (xg, yg)
expected_fs_raw = Dict(1 => 0.451, -1 => 0.549)

sampler, = statespace_sampler(Random.MersenneTwister(1234);
min_bounds = minimum.(grid), max_bounds = maximum.(grid)
)
ics = StateSpaceSet([sampler() for i in 1:1000])

fs, approx_atts, labels = basins_fractions(mapper, ics; show_progress = false)
@test length(keys(fs)) == 2
@test fs[1] ≈ 0.451 rtol = 1e-2
# the divergent points go to last bin, which is 25 = 5x5
@test fs[25] ≈ 0.549 rtol = 1e-2