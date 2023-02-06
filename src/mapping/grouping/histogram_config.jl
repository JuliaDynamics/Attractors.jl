using ComplexityMeasures: FixedRectangularBinning, probabilities_and_outcomes
export GroupViaHistogram
export FixedRectangularBinning

"""
    GroupViaHistogram(binning::FixedRectangularBinning)

Initialize a struct that contains instructions on how to group features in
[`AttractorsViaFeaturizing`](@ref). `GroupViaHistogram` performs a histogram
in feature space. Then, all features that are in the same histogram bin get the
same label. The `binning` is an instance of [`FixedRectangularBinning`](@ref)
from ComplexityMeasures.jl. (the reason to not allow `RectangularBinning` is because
during continuation we need to ensure that bins remain identical).
"""
struct GroupViaHistogram{E<:RectangularBinEncoding} <: GroupingConfig
    encoding::E
end

function GroupViaHistogram(binning::FixedRectangularBinning)
    enc = RectangularBinEncoder(binning)
    return GroupViaHistogram(enc)
end

function feature_to_group(feature, config::GroupViaHistogram)
    enc = config.encoding
    return encode(enc, feature)
end
