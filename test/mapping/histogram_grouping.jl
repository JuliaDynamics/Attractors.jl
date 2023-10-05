using Test
using Attractors
import Random
using Statistics: mean

henon_rule(x, p, n) = SVector{2}(1.0 - p[1]*x[1]^2 + x[2], p[2]*x[1])
henon() = DeterministicIteratedMap(henon_rule, zeros(2), [1.4, 0.3])
ds = henon()
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

sampler, _ = statespace_sampler(grid, 12444)

ics = StateSpaceSet([copy(sampler()) for _ in 1:1000])
fs, = basins_fractions(mapper, ics; show_progress = false)
@test length(keys(fs)) == 2
@test fs[1] ≈ 0.45 rtol = 1e-1
# the divergent points go to last bin, which is 25 = 5x5
@test fs[25] ≈ 0.55 rtol = 1e-1