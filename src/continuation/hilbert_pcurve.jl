export hilbert_pcurve
import BijectiveHilbert # must be 0.6.1 or later

"""
    hilbert_pcurve(input) → pcurve

Create a parameter curve to perform global continuation over by constructing a Hilbert
curve in a multiparameter space. The input specifies which parameters to cover,
and it must be dictionary mapping parameter indices to (min, max, density) tuples.
For example,
```julia
input = Dict(1 => (0, 5, 32), 2 => (2.2, 3.4, 2^5), 6 => (-2, 2, 2^7))
```
will create a 3-dimensional Hilbert curve over the 1st, 2nd, and 6th parameter.
The 1st parameter for example will cover the range `0` to `5` using `32` values.

!!! note "Density is in powers of 2"
    Due to the constraints in multidimensional Hilbert curve construction, the provided
    densities are always converted to the closest power of 2. It is recommended to
    provided powers of 2 to begin with, to avoid confusion.
"""
hilbert_pcurve(input::AbstractDict) = hilbert_pcurve(collect(input))
function hilbert_pcurve(input::Vector{<:Pair}) # need ordered input for hilbert
    # transform to hilbert appropriate
    pidxs = getindex.(input, 1)
    powers2 = [round(Int, log2(v[3])) for (k, v) in input]
    pranges = [range(v[1], v[2]; length = 2^powers2[i]) for (i, (k, v)) in enumerate(input)]
    # produce hilbert parameter curve
    hencoder = BijectiveHilbert.Compact(Int, powers2)
    intseq = 1:prod(2 .^ powers2)
    point_out = zeros(Int, length(input))
    # output however can be `Dict` so that parameter indices can be used
    hilbert_curve = map(intseq) do i
        BijectiveHilbert.decode_hilbert!(hencoder, point_out, i)
        pvals = getindex.(pranges, point_out)
        return Dict(k => v for (k, v) in zip(pidxs, pvals))
    end
    return hilbert_curve
end
