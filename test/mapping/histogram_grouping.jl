using Test
using Attractors

ds = Systems.henon(zeros(2); a = 1.4, b = 0.3)
xg = yg = range(-2.0, 2.0; length=100)
grid = (xg, yg)
expected_fs_raw = Dict(1 => 0.451, -1 => 0.549)
function featurizer(A, t)
    # Notice that unsupervised clustering cannot support "divergence to infinity",
    # which it identifies as another attractor (in fact, the first one).
    x = @SVector[mean(A[:, 1]), mean(A[:, 2])]
    return any(>(10), x) ? @SVector[200.0, 200.0] : x
end

binning = FixedRectangularBinning((-10, 290; step = 100), 2)