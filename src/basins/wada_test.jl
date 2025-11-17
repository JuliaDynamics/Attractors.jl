export test_wada_merge, haussdorff_distance



"""
    test_wada_merge(BoA::ArrayBasinsOfAttraction,r) -> p
    test_wada_merge(basins, r) -> p

Test if the 2D array `basins` has the [Wada
property](https://en.wikipedia.org/wiki/Lakes_of_Wada)
using the merging technique of [Daza2018](@cite).

The second function signature exists for backwards compatibility. 

## Description

The technique consists in computing the generalized basins
of each attractor. These new basins are formed with on of the basins
and the union of the other basins. A new boundary is defined by
these two objects. The algorithm then computes the distance between
each boundaries of these basins pairwise. If all the boundaries are
within some distance `r`, there is a unique boundary separating
the basins and we have the wada property.
The algorithm returns the maximum proportion of pixels of a boundary
with distance strictly greater than `r` from another boundary.

If `p == 0`,  we have the Wada property for this value of `r`.
If `p > 0`, the criteria to decide if
the basins are Wada is left to the user. Numerical inaccuracies
may be responsible for a small percentage of points with distance larger than `r`

"""
function test_wada_merge(basins,r)
    ids = unique(basins)
    if length(ids) < 3
        @error "There must be at least 3 attractors"
        return  Inf
    end
    M = Vector{BitMatrix}(undef,length(ids))
    for k in ids
       M[k] = merged_basins(basins,k)
    end
    v = Vector{Float64}()
    for k in 1:length(M)-1, j in k+1:length(M)
        push!(v,wada_fractions(M[k],M[j],r))
    end
    return maximum(v)
end

test_wada_merge(BoA::ArrayBasinsOfAttraction,r) = test_wada_merge(BoA.basins,r)

"""
    haussdorff_distance(M1, M2) -> hd

Compute the Hausdorff distance between two binary matrices of
the same size. First a distance matrix is computed using
`distance_matrix` for each matrix M1 and M2. The entries
of this matrix are the distance in L1 metric to the
closest 0 pixel in the initial matrix. The distance being
0 if it belongs to the basin with id 0.
"""
function haussdorff_distance(M1::BitMatrix, M2::BitMatrix)
    bd1 = distance_matrix(M1)
    bd2 = distance_matrix(M2)
    hd1 = maximum(bd1.*M2)
    hd2 = maximum(bd2.*M1)
    return max(hd1,hd2)
end

# wada_fractions computes the distance between
# boundaries using the hausdorff distance.
# The distance is computed using the Shonwilker
# algorithm. The result is a matrix with the distance
# of the pixels from one boundary to the other.
# The algorithm returns the number of pixels
# that are at a distance strictly more than r
# (with the hausdorff distance).
function wada_fractions(bas1::BitMatrix, bas2::BitMatrix, r::Int)
    bnd1 = get_boundary(bas1)
    bd1 = distance_matrix(bnd1)

    bnd2 = get_boundary(bas2)
    bd2 = distance_matrix(bnd2)

    c1 = count(bd1.*bnd2 .> r)
    c2 = count(bd2.*bnd1 .> r)
    return max(c1,c2)./length(bnd1)
end



# Return a matrix with two basins: the first is the basins
# from the attractor in `id`  and the other basin is formed
# with the union of all others basins.
# The function returns a BitMatrix such that the basins of
# `id` is mapped to 0 and the other basins to 1.
function merged_basins(basins, id)
    mrg_basins = fill!(BitMatrix(undef,size(basins)), false)
    ids = setdiff(unique(basins), id)
    for k in ids
        I = findall(basins .== k)
        mrg_basins[I] .= 1
    end
    return mrg_basins
end



# Generate all pairs of the number in ids without
# repetition
function generate_pairs(ids)
 p = Vector{Tuple{Int,Int}}()
 for (k,n1) in enumerate(ids)
     for (j,n2) in enumerate(ids[k+1:end])
        push!(p,(n1,n2))
     end
 end
 return p
end


# This function returns the boundary between two basins:
# First we compute the distance matrix of the basins and its
# complementary. The pixels at distance 1 from a basin are
# the closest to the boundary. We compute the union the
# pixels computed from one basins and its complementary.
# We get the boundary in L1 distance!
function get_boundary(basins::BitMatrix)
    bd1 = distance_matrix(basins)
    bdd1 = distance_matrix( .! basins)
    bnd = (bd1 .== 1) .| (bdd1 .== 1)
    return bnd
end


# Function for L1 metric
w(M) = min(M[1,2]+1, M[2,1]+1, M[2,2])
w2(M) = min(M[1,2]+1, M[2,1]+1, M[1,1])

# R. Shonkwilker, An image algorithm for computing the Hausdorff distance efficiently in linear time.
# https://doi.org/10.1016/0020-0190(89)90114-2
# This compute the Distance transform matrix.
# Given a matrix, we compute the distance from the
# basin with label 0. The result is a matrix whose
# entry is the distance to the closest 0 pixel in the
# L1 metric (Manhattan).
function distance_matrix(basins::BitMatrix)
   r,c = size(basins)
   basdist = ones(Int32,r+2,c+2)*(r+c+4)
   # Assign the maximum distance to the pixels not in the basin
   basdist[2:r+1,2:c+1] .= (1 .- basins) .*(r+c+4)
   # first pass right to left, up to bottom
   for j in 2:c+1, k in 2:r+1
       basdist[j,k] = w(view(basdist,j-1:j,k-1:k))
    end
   # second pass left to right, bottom up
   for j in c+1:-1:1, k in r+1:-1:1
       basdist[j,k] = w2(view(basdist,j:j+1,k:k+1))
   end
   # Remove the extra rows and cols necessary to the alg.
   return basdist[2:r+1,2:c+1]
end
