# Analysis layer over feature vectors: distances, design-complexity scores,
# and per-cluster complexities. K-means clustering and MDS embedding live in
# package extensions (load Clustering.jl / MultivariateStats.jl).

"""
    feature_matrix(analysis::HarmonicAnalysis) -> Matrix

Feature vectors as a `dims × n_nodes` matrix (one column per node, in model
order) — the layout Clustering.jl and MultivariateStats.jl expect.
"""
feature_matrix(ha::HarmonicAnalysis) = reduce(hcat, ha.featurevectors)

"""
    distance_matrix(analysis::HarmonicAnalysis)
    distance_matrix(featurevectors)

Symmetric `n × n` matrix of pairwise Euclidean distances between nodal
feature vectors — the inter-node dissimilarity measure of Lee, Danhaive &
Mueller (2022), eq. 7.
"""
function distance_matrix(fvs::AbstractVector{<:AbstractVector})
    n = length(fvs)
    T = float(eltype(eltype(fvs)))
    D = zeros(T, n, n)

    for i = 1:n, j = i+1:n
        D[i, j] = D[j, i] = norm(fvs[i] - fvs[j])
    end

    return D
end

distance_matrix(ha::HarmonicAnalysis) = distance_matrix(ha.featurevectors)

# minimal enclosing hypersphere: Welzl's algorithm with move-to-front.
# The circumball of an affinely independent boundary set B = {p₀, …, pₖ} has
# its center in aff(B): c = p₀ + Σₐ λₐ qₐ, qₐ = pₐ - p₀, where the boundary
# conditions |c - pₐ| = |c - p₀| give the SPD system 2⟨qₐ,q_b⟩λ = |qₐ|².
struct _Ball{T}
    center::Vector{T}
    r2::T
end

function _circumball(B::Vector{Vector{T}}) where {T}
    isempty(B) && return _Ball(T[], T(-Inf))

    p0 = B[1]
    k = length(B) - 1
    k == 0 && return _Ball(copy(p0), zero(T))

    Q = [B[a+1] - p0 for a = 1:k]
    G = [2 * dot(Q[a], Q[b]) for a = 1:k, b = 1:k]
    rhs = [dot(q, q) for q in Q]
    λ = G \ rhs

    c = p0 + sum(λ[a] * Q[a] for a = 1:k)
    return _Ball(c, sum(abs2, c - p0))
end

function _inball(b::_Ball, p::AbstractVector, tol::Real)
    isempty(b.center) && return false
    return sum(abs2, p - b.center) <= b.r2 + tol * (1 + b.r2)
end

function _welzl!(points::Vector{Vector{T}}, n::Int, boundary::Vector{Vector{T}}, d::Int, tol::Real) where {T}
    b = _circumball(boundary)
    length(boundary) == d + 1 && return b

    for i = 1:n
        p = points[i]
        if !_inball(b, p, tol)
            b = _welzl!(points, i - 1, vcat(boundary, [p]), d, tol)

            # move-to-front: violators are likely on the final boundary
            for j = i:-1:2
                points[j] = points[j-1]
            end
            points[1] = p
        end
    end

    return b
end

"""
    bounding_sphere(points) -> (center, radius)

Exact minimal enclosing hypersphere of a set of equal-length vectors, via
Welzl's move-to-front algorithm (expected O(n) for fixed dimension). Used by
[`complexity`](@ref); exported for direct use on feature-vector sets.
"""
function bounding_sphere(points_in::AbstractVector{<:AbstractVector}; tol::Real = 1e-10)
    isempty(points_in) && throw(ArgumentError("bounding_sphere of an empty point set"))

    points = [convert(Vector{Float64}, collect(p)) for p in points_in]
    d = length(points[1])
    all(p -> length(p) == d, points) ||
        throw(DimensionMismatch("points must all have the same length"))

    b = _welzl!(points, length(points), Vector{Vector{Float64}}(), d, tol)
    return (center = b.center, radius = sqrt(max(b.r2, 0.0)))
end

"""
    complexity(analysis::HarmonicAnalysis)
    complexity(featurevectors)

Design-complexity score of Lee, Danhaive & Mueller (2022): the radius of the
minimal bounding hypersphere of the nodal feature vectors. Comparative only —
larger radius ⇒ greater variation in nodal force demands ⇒ higher penalty for
connection standardization.

This is a minimax quantity and is not smooth in the underlying forces; for
gradient-based optimization use [`soft_complexity`](@ref).
"""
complexity(fvs::AbstractVector{<:AbstractVector}) = bounding_sphere(fvs).radius
complexity(ha::HarmonicAnalysis) = complexity(ha.featurevectors)

"""
    soft_complexity(analysis::HarmonicAnalysis)
    soft_complexity(featurevectors)

Smooth surrogate for [`complexity`](@ref): the root-mean-square distance of
the nodal feature vectors from their centroid,
`√(1/n · Σᵢ ‖FVᵢ - F̄V‖²)`. Differentiable in the feature vectors (and hence
in force magnitudes/directions), scales linearly with force magnitudes, and is
bounded by `2·complexity`. Intended as the objective for minimizing nodal
demand variation with AsapOptim.
"""
function soft_complexity(fvs::AbstractVector{<:AbstractVector})
    n = length(fvs)
    centroid = sum(fvs) / n
    return sqrt(sum(sum(abs2, fv - centroid) for fv in fvs) / n)
end

soft_complexity(ha::HarmonicAnalysis) = soft_complexity(ha.featurevectors)

"""
    cluster_complexities(analysis, assignments)
    cluster_complexities(featurevectors, assignments)

Per-cluster complexity: the bounding-sphere radius of each cluster's feature
vectors, given integer cluster `assignments` (e.g.
`cluster_nodes(analysis, k).assignments`). Returned in order of cluster label
`1:k`. Identical demands within a cluster give radius 0 — the
standardization-penalty measure used for connection rationalization.
"""
function cluster_complexities(fvs::AbstractVector{<:AbstractVector}, assignments::AbstractVector{<:Integer})
    length(fvs) == length(assignments) ||
        throw(DimensionMismatch("$(length(fvs)) feature vectors for $(length(assignments)) assignments"))

    return [complexity(fvs[findall(==(c), assignments)]) for c in sort!(unique(assignments))]
end

cluster_complexities(ha::HarmonicAnalysis, assignments::AbstractVector{<:Integer}) =
    cluster_complexities(ha.featurevectors, assignments)

"""
    cluster_nodes(analysis, k; kwargs...)

K-means clustering of the nodal feature vectors into `k` groups for connection
standardization. **Requires Clustering.jl** (`using Clustering`) — implemented
in a package extension. Returns the `Clustering.KmeansResult` (cluster labels
in `.assignments`); keyword arguments pass through to `Clustering.kmeans`.
"""
function cluster_nodes end

"""
    embed_nodes(analysis; maxoutdim = 2)

Classical multidimensional-scaling embedding of the nodes into
`maxoutdim`-dimensional space, best preserving the feature-vector distance
matrix — the demand-space visualization of Lee, Danhaive & Mueller (2022).
**Requires MultivariateStats.jl** (`using MultivariateStats`) — implemented in
a package extension. Returns a `maxoutdim × n_nodes` coordinate matrix.
"""
function embed_nodes end

"""
    harmonic_params(p; delta = 20, dims = 16, dimension = 3)

Precompile the constant (non-differentiable) data of a shape-descriptor
evaluation over an `AsapOptim.OptParams`: per-node incident elements and
orientation signs, applied-load bumps, support flags, and the applied-load
matrix for equilibrium reaction recovery. Feed the result to
[`feature_vectors`](@ref), [`soft_complexity`](@ref), and
[`complexity`](@ref). **Requires AsapOptim.jl** (`using AsapOptim`) —
implemented in a package extension.
"""
function harmonic_params end

"""
    feature_vectors(res, p, hp)
    feature_vectors(x, p, hp)

Fully differentiable nodal feature vectors of a design evaluation: member
directions from the evaluated positions (`res.X` + incidence), member forces
via `AsapOptim.axial_force`, and support reactions recovered from nodal
equilibrium (no extra solve). `res = solve_structure(x, p)`; `hp` from
[`harmonic_params`](@ref). Gradients flow with ForwardDiff or Zygote
(`soft_complexity(x, p, hp)` is the intended smooth objective).
**Requires AsapOptim.jl** — implemented in a package extension.
"""
function feature_vectors end
