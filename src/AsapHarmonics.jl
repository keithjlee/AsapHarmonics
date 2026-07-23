module AsapHarmonics

using Asap, LinearAlgebra, SparseArrays, StaticArrays
using FastSphericalHarmonics, FFTW
using SpecialFunctions: besselix

include("descriptors.jl")
export spherical_feature_vector, circular_feature_vector
export zonal_coefficients, pairwise_legendre_sums

include("functions.jl")
export spherical_gaussian, circular_gaussian, sampled_force_function
export make_xsphere, make_ysphere, make_zsphere

include("types.jl")
export NodeSignature, feature_vector
export HarmonicAnalysis, HarmonicAnalysis2d

end # module AsapHarmonics
