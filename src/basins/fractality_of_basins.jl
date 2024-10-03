export uncertainty_exponent, basins_fractal_dimension, basins_fractal_test, basin_entropy

"""
    basin_entropy(basins::Array, ε = 20) -> Sb, Sbb

Return the basin entropy [Daza2016](@cite) `Sb` and basin boundary entropy `Sbb`
of the given `basins` of attraction by considering `ε`-sized boxes along each dimension.

## Description

First, the n-dimensional input `basins`
is divided regularly into n-dimensional boxes of side `ε`.
If `ε` is an integer, the same size is used for all dimensions, otherwise `ε` can be
a tuple with the same size as the dimensions of `basins`.
Assuming that there are ``N`` `ε`-boxes that cover the `basins`, the basin entropy is estimated
as [Daza2016](@cite)

```math
S_b = \\tfrac{1}{N}\\sum_{i=1}^{N}\\sum_{j=1}^{m_i}-p_{ij}\\log(p_{ij})
```
where ``m`` is the number of unique IDs (integers of `basins`) in box ``i``
and ``p_{ij}`` is the relative frequency (probability) to obtain ID ``j``
in the ``i`` box (simply the count of IDs ``j`` divided by the total in the box).

`Sbb` is the boundary basin entropy.
This follows the same definition as ``S_b``, but now averaged over only
only boxes that contains at least two
different basins, that is, for the boxes on the boundaries.

The basin entropy is a measure of the uncertainty on the initial conditions of the basins.
It is maximum at the value `log(n_att)` being `n_att` the number of unique IDs in `basins`. In
this case the boundary is intermingled: for a given initial condition we can find
another initial condition that lead to another basin arbitrarily close. It provides also
a simple criterion for fractality: if the boundary basin entropy `Sbb` is above `log(2)`
then we have a fractal boundary. It doesn't mean that basins with values below cannot
have a fractal boundary, for a more precise test see [`basins_fractal_test`](@ref).
An important feature of the basin entropy is that it allows
comparisons between different basins using the same box size `ε`.
"""
function basin_entropy(basins::AbstractArray{<:Integer, D}, ε::Integer = 20) where {D}
    es = ntuple(i -> ε, Val(D))
    return basin_entropy(basins, es)
end

function basin_entropy(basins::AbstractArray{<:Integer, D}, es::NTuple{D, <: Integer}) where {D}
    Sb = 0.0; Nb = 0
    εranges = map((d, ε) -> 1:ε:d, size(basins), es)
    box_iterator = Iterators.product(εranges...)
    for box_start in box_iterator
        box_ranges = map((d, ε) -> d:(d+ε-1), box_start, es)
        box_values = view(basins, box_ranges...)
        uvals = unique(box_values)
        if length(uvals) > 1
            Nb += 1
            # we only need to estimate entropy for boxes with more than 1 val,
            # because in other cases the entropy is zero
            Sb = Sb + _box_entropy(box_values, uvals)
        end
    end
    return Sb/length(box_iterator), Sb/Nb
end

function _box_entropy(box_values, unique_vals = unique(box_values))
    h = 0.0
    for v in unique_vals
        p = count(x -> (x == v), box_values)/length(box_values)
        h += -p*log(p)
    end
    return h
end



"""
    basins_fractal_test(basins; ε = 20, Ntotal = 1000) -> test_res, Sbb

Perform an automated test to decide if the boundary of the basins has fractal structures
based on the method of Puy et al. [Puy2021](@cite).
Return `test_res` (`:fractal` or `:smooth`) and the mean basin boundary entropy.

## Keyword arguments

* `ε = 20`: size of the box to compute the basin boundary entropy.
* `Ntotal = 1000`: number of balls to test in the boundary for the computation of `Sbb`

## Description

The test "looks" at the basins with a magnifier of size `ε` at random.
If what we see in the magnifier looks like a smooth boundary (onn average) we decide that
the boundary is smooth. If it is not smooth we can say that at the scale `ε` we have
structures, i.e., it is fractal.

In practice the algorithm computes the boundary basin entropy `Sbb` [`basin_entropy`](@ref)
for `Ntotal`
random boxes of radius `ε`. If the computed value is equal to theoretical value of a smooth
boundary
(taking into account statistical errors and biases) then we decide that we have a smooth
boundary. Notice that the response `test_res` may depend on the chosen ball radius `ε`.
For larger size,
we may observe structures for smooth boundary and we obtain a *different* answer.

The output `test_res` is a symbol describing the nature of the basin and the output `Sbb` is
the estimated value of the boundary basin entropy with the sampling method.
"""
function basins_fractal_test(basins; ε = 20, Ntotal = 1000)
    dims = size(basins)
    # Sanity check.
    if minimum(dims)/ε < 50
        @warn "Maybe the size of the grid is not fine enough."
    end
    if Ntotal < 100
        error("Ntotal must be larger than 100 to gather enough statistics.")
    end

    v_pts = zeros(Float64, length(dims), prod(dims))
    I = CartesianIndices(basins)
    for (k,coord) in enumerate(I)
         v_pts[:, k] = [Tuple(coord)...]
    end
    tree = searchstructure(KDTree, v_pts, Euclidean())
    # Now get the values in the boxes.
    Nb = 1
    N_stat = zeros(Ntotal)
    while Nb < Ntotal
        p = [rand()*(sz-ε)+ε for sz in dims]
        idxs = isearch(tree, p, WithinRange(ε))
        box_values = basins[idxs]
        bx_ent = _box_entropy(box_values)
        if bx_ent > 0
            Nb = Nb + 1
            N_stat[Nb] = bx_ent
        end
    end

    Ŝbb = mean(N_stat)
    σ_sbb = std(N_stat)/sqrt(Nb)
    # Table of boundary basin entropy of a smooth boundary for dimension 1 to 5:
    Sbb_tab = [0.499999, 0.4395093, 0.39609176, 0.36319428, 0.33722572]
    if length(dims) ≤ 5
        Sbb_s = Sbb_tab[length(dims)]
    else
        Sbb_s = 0.898*length(dims)^-0.4995
    end
    # Systematic error approximation for the disk of radius ε
    δub = 0.224*ε^-1.006

    tst_res = :smooth
    if Ŝbb < (Sbb_s - σ_sbb) ||  Ŝbb > (σ_sbb + Sbb_s + δub)
        # println("Fractal boundary for size of box ε=", ε)
        tst_res = :fractal
    else
        # println("Smooth boundary for size of box ε=", ε)
        tst_res = :smooth
    end

    return tst_res, Ŝbb
end

# as suggested in https://github.com/JuliaStats/StatsBase.jl/issues/398#issuecomment-417875619
linreg(x, y) = hcat(fill!(similar(x), 1), x) \ y

"""
    basins_fractal_dimension(basins; kwargs...) -> V_ε, N_ε, d

Estimate the fractal dimension `d` of the boundary between basins of attraction using
a box-counting algorithm for the boxes that contain at least two different basin IDs.

## Keyword arguments

* `range_ε = 2:maximum(size(basins))÷20` is the range of sizes of the box to
  test (in pixels).

## Description

The output `N_ε` is a vector with the number of the balls of radius `ε` (in pixels)
that contain at least two initial conditions that lead to different attractors. `V_ε`
is a vector with the corresponding size of the balls. The output `d` is the estimation
of the box-counting dimension of the boundary by fitting a line in the `log.(N_ε)`
vs `log.(1/V_ε)` curve. However it is recommended to analyze the curve directly
for more accuracy.

It is the implementation of the popular algorithm of the estimation of the box-counting
dimension. The algorithm search for a covering the boundary with `N_ε` boxes of size
`ε` in pixels.
"""
function basins_fractal_dimension(basins::AbstractArray; range_ε = 3:maximum(size(basins))÷20)

    dims = size(basins)
    num_step = length(range_ε)
    N_u = zeros(Int, num_step) # number of uncertain box
    N = zeros(Int, num_step) # number of boxes
    V_ε = zeros(1, num_step) # resolution

    # Naive box counting estimator
    for (k,eps) in enumerate(range_ε)
        Nb, Nu = 0, 0
        # get indices of boxes
        bx_tuple = ntuple(i -> range(1, dims[i] - rem(dims[i],eps), step = eps), length(dims))
        box_indices = CartesianIndices(bx_tuple)
        for box in box_indices
            # compute the range of indices for the current box
            ind = CartesianIndices(ntuple(i -> range(box[i], box[i]+eps-1, step = 1), length(dims)))
            c = basins[ind]
            if length(unique(c))>1
                Nu = Nu + 1
            end
            Nb += 1
        end
        N_u[k] = Nu
        N[k] = Nb
        V_ε[k] = eps
    end
    N_ε = N_u
    # remove zeros in case there are any:
    ind = N_ε .> 0.0
    N_ε = N_ε[ind]
    V_ε = V_ε[ind]
    # get exponent via liner regression on `f_ε ~ ε^α`
    b, d = linreg(vec(-log10.(V_ε)), vec(log10.(N_ε)))
    return V_ε, N_ε, d
end

"""
    uncertainty_exponent(basins; kwargs...) -> ε, N_ε, α

Estimate the uncertainty exponent[Grebogi1983](@cite) of the basins of attraction. This exponent
is related to the final state sensitivity of the trajectories in the phase space.
An exponent close to `1` means basins with smooth boundaries whereas an exponent close
to `0` represent completely fractalized basins, also called riddled basins.

The output `N_ε` is a vector with the number of the balls of radius `ε` (in pixels)
that contain at least two initial conditions that lead to different attractors.
The output `α` is the estimation of the uncertainty exponent using the box-counting
dimension of the boundary by fitting a line in the `log.(N_ε)` vs `log.(1/ε)` curve.
However it is recommended to analyze the curve directly for more accuracy.

## Keyword arguments
* `range_ε = 2:maximum(size(basins))÷20` is the range of sizes of the ball to
  test (in pixels).

## Description

A phase space with a fractal boundary may cause a uncertainty on the final state of the
dynamical system for a given initial condition. A measure of this final state sensitivity
is the uncertainty exponent. The algorithm probes the basin of attraction with balls
of size `ε` at random. If there are a least two initial conditions that lead to different
attractors, a ball is tagged "uncertain". `f_ε` is the fraction of "uncertain balls" to the
total number of tries in the basin. In analogy to the fractal dimension, there is a scaling
law between, `f_ε ~ ε^α`. The number that characterizes this scaling is called the
uncertainty exponent `α`.

Notice that the uncertainty exponent and the box counting dimension of the boundary are
related. We have `Δ₀ = D - α` where `Δ₀` is the box counting dimension computed with
[`basins_fractal_dimension`](@ref) and `D` is the dimension of the phase space.
The algorithm first estimates the box counting dimension of the boundary and
returns the uncertainty exponent.
"""
function uncertainty_exponent(basins::AbstractArray; range_ε = 2:maximum(size(basins))÷20)
    V_ε, N_ε, d = basins_fractal_dimension(basins; range_ε)
    return V_ε, N_ε, length(size(basins)) - d
end
