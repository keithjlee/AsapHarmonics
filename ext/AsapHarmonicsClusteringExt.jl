module AsapHarmonicsClusteringExt

using AsapHarmonics, Clustering

function AsapHarmonics.cluster_nodes(ha::HarmonicAnalysis, k::Integer; kwargs...)
    return kmeans(feature_matrix(ha), k; kwargs...)
end

function AsapHarmonics.cluster_nodes(fvs::AbstractVector{<:AbstractVector}, k::Integer; kwargs...)
    return kmeans(reduce(hcat, fvs), k; kwargs...)
end

end # module
