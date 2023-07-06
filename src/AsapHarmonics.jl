module AsapHarmonics

using Asap, LinearAlgebra
using FastSphericalHarmonics

l_max = 25
resolution = 90

include("init.jl")
export Y
export xsphere, ysphere, zsphere

include("functions.jl")
export generate_feature_vector

include("types.jl")
export NodeForce, HarmonicAnalysis

end # module AsapHarmonics
