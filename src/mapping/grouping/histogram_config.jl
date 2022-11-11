# using Entropies: FixedRectangularBinning, probabilities_and_outcomes
# export GroupViaHistogram
# export FixedRectangularBinning

# """
#     GroupViaHistogram(binning, min_entries_per_bin = 1)

# Initialize a struct that contains instructions on how to group features in
# [`AttractorsViaFeaturizing`](@ref). `GroupViaHistogram` performs a histogram
# in feature space. Then, all features that are in the same histogram bin get the
# same label. The `binning` is an instance of [`FixedRectangularBinning`](@ref)
# from Entropies.jl. (the reason to not allow `RectangularBinning` is because
# during continuation we need to ensure that bins remain identical).

# The second argument `min_entries_per_bin` quantifies how many entries each bin
# needs to have to get assigned a unique label. All bins that have _less_ than
# `min_entries_per_bin` get assigned the special label `-1`.
# """
# struct GroupViaHistogram{B<:FixedRectangularBinning, L} <: GroupingConfig
#     binning::B
#     min_entries_per_bin::Int
#     labels_to_bins::L
# end

# GroupViaHistogram(binning) = GroupViaHistogram(binning, 1)

# function group_features(features::Vector{<:AbstractVector}, config::GroupViaHistogram)
#     error("Not yet implemented (waiting for resolution in Entropies.jl)")
# end