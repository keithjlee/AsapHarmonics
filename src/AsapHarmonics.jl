module AsapHarmonics

using Asap, LinearAlgebra, SparseArrays
using FastSphericalHarmonics, FFTW

l_max = 25
resolution = 90

include("init.jl")
export Y
export xsphere, ysphere, zsphere

include("functions.jl")
export spherical_feature_vector

include("types.jl")
export NodeForces, HarmonicAnalysis
export NodeForces2d, HarmonicAnalysis2d

end # module AsapHarmonics
