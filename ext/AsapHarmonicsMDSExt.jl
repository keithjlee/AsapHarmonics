module AsapHarmonicsMDSExt

using AsapHarmonics, MultivariateStats

function AsapHarmonics.embed_nodes(ha::HarmonicAnalysis; maxoutdim::Integer = 2)
    D = distance_matrix(ha)
    M = fit(MDS, D; distances = true, maxoutdim = maxoutdim)
    return predict(M)
end

function AsapHarmonics.embed_nodes(fvs::AbstractVector{<:AbstractVector}; maxoutdim::Integer = 2)
    D = distance_matrix(fvs)
    M = fit(MDS, D; distances = true, maxoutdim = maxoutdim)
    return predict(M)
end

end # module
