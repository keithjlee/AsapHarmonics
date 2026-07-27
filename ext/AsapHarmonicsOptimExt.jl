# Differentiable shape-descriptor pipeline over AsapOptim design evaluations:
# design vector x -> solve_structure -> per-node signatures -> closed-form
# feature vectors -> (soft) complexity.
#
# The evaluation is fully BATCHED: reverse-mode engines pay one pullback per
# array operation, and Zygote's adjoint for indexing allocates a parent-sized
# zero buffer per access — a per-node gather loop costs O(n_nodes · n_el)
# memory traffic and dominates the gradient (~97% measured). Instead, all
# constant structure is precompiled into sparse selection/aggregation
# matrices ("bumps" = every (node, force) pair; "pairs" = every within-node
# bump pair), and the whole model evaluates as a fixed, small number of
# matmuls and broadcasts:
#
#   member dirs  D_mem = Ê · G_mem          member mags  m_mem = S_mem · f
#   reactions    R = M .* ((Ê·diag(f))·Cinc)   (elastic member forces summed at
#                fixed DOFs — Asap's reaction() convention: applied loads at
#                supports are NOT subtracted; M masks free DOFs per component)
#   pair dots    t = colsum(D[:,I] .* D[:,J]),  pair weights w = m[I] .* m[J]
#   band sums    Sₗ = A_pair · (w .* Pₗ(t))     (Legendre 3D / Chebyshev 2D)
#   descriptors  FV = prefactor .* √Sₗ          (dims × n_nodes matrix)

module AsapHarmonicsOptimExt

using AsapHarmonics, AsapOptim, Asap
using LinearAlgebra, SparseArrays
import AsapHarmonics: harmonic_params, feature_matrix, feature_vectors, soft_complexity,
    complexity, cluster_projector, zonal_coefficients, _band_norm, _constmul, _mulconst

"""
Constant data of a shape-descriptor evaluation over an `OptParams` — see
`AsapHarmonics.harmonic_params`. Bump order is [members; loads; reactions]
(members only under `weights = :unit`).
"""
struct HarmonicOptParams
    dimension::Int                     # 2 (circle) or 3 (sphere)
    delta::Float64                     # kernel parameter (δ in 3D, σ in 2D)
    dims::Int                          # feature-vector length
    weights::Symbol                    # :force (demands) or :unit (geometry only)
    prefactor::Vector{Float64}         # FVₗ = prefactor[l+1]·√Sₗ
    Gmem::SparseMatrixCSC{Float64,Int} # n_el × M_mem: D_mem = Ê · Gmem (cols = -sign·Êₑ)
    Smem::SparseMatrixCSC{Float64,Int} # M_mem × n_el: m_mem = Smem · f
    Dload::Matrix{Float64}             # dimension × M_load constant unit load directions
    mload::Vector{Float64}             # M_load constant load magnitudes
    Grxn::SparseMatrixCSC{Float64,Int} # n_nodes × M_rxn: raw reactions = R · Grxn
    fixmask::Matrix{Float64}           # 3 × n_nodes: 1.0 at fixed translational DOFs
    pairI::Vector{Int}                 # within-node bump pairs (unordered, i ≤ j)
    pairJ::Vector{Int}
    Apair::SparseMatrixCSC{Float64,Int} # n_nodes × n_pairs, weight 1 (diag) / 2 (off-diag)
    nnodes::Int
end

function harmonic_params(p::OptParams; delta::Real = 20, dims::Integer = 16,
    dimension::Integer = 3, weights::Symbol = :force)
    dimension == 2 || dimension == 3 ||
        throw(ArgumentError("dimension must be 2 or 3, got $dimension"))
    weights == :force || weights == :unit ||
        throw(ArgumentError("weights must be :force or :unit for optimization " *
                            "(:sign is not differentiable), got :$weights"))

    model = p.model
    n = length(model.nodes)
    nel = length(model.elements)
    d = Int(dimension)

    bumps = [Int[] for _ = 1:n] # global bump indices per node
    nbumps = 0

    # members: one bump per (node, incident element)
    gI = Int[]; gJ = Int[]; gV = Float64[]   # Gmem triplets (element, bump, -sign)
    sI = Int[]; sJ = Int[]                   # Smem triplets (bump, element)
    rows = rowvals(p.Cinc)
    vals = nonzeros(p.Cinc)
    for i = 1:n, k in nzrange(p.Cinc, i)
        nbumps += 1
        push!(gI, rows[k]); push!(gJ, nbumps); push!(gV, -vals[k])
        push!(sI, nbumps); push!(sJ, rows[k])
        push!(bumps[i], nbumps)
    end
    Mmem = nbumps
    Gmem = sparse(gI, gJ, gV, nel, Mmem)
    Smem = sparse(sI, sJ, ones(Mmem), Mmem, nel)

    # applied point forces: constant bumps (projected per dimension); none
    # under :unit — the geometry-only signature is members only
    loaddirs = Vector{Float64}[]
    mload = Float64[]
    if weights == :force
        for load in model.loads
            load isa NodeForce || continue
            i = load.node.index
            P = collect(Float64, load.value)[1:d]
            m = norm(P)
            iszero(m) && continue
            nbumps += 1
            push!(loaddirs, P ./ m)
            push!(mload, m)
            push!(bumps[i], nbumps)
        end
    end
    Dload = isempty(loaddirs) ? zeros(d, 0) : reduce(hcat, loaddirs)

    # support reactions: one bump per node with a fixed translational DOF
    # (Asap fixity: true = free); fixmask selects the reacting components
    rI = Int[]; rJ = Int[]
    fixmask = zeros(3, n)
    nrxn = 0
    if weights == :force
        for (i, node) in enumerate(model.nodes)
            fixmask[:, i] .= .!node.fixity[1:3]
            any(.!node.fixity[1:3]) || continue
            nbumps += 1
            nrxn += 1
            push!(rI, i); push!(rJ, nrxn)
            push!(bumps[i], nbumps)
        end
    end
    Grxn = sparse(rI, rJ, ones(nrxn), n, nrxn)

    # within-node bump pairs; off-diagonal pairs carry weight 2 in the full
    # double sum Sₗ = Σᵢⱼ fᵢfⱼPₗ(n̂ᵢ·n̂ⱼ)
    pairI = Int[]; pairJ = Int[]
    aI = Int[]; aV = Float64[]
    for i = 1:n
        b = bumps[i]
        for a = 1:length(b), c = a:length(b)
            push!(pairI, b[a]); push!(pairJ, b[c])
            push!(aI, i); push!(aV, a == c ? 1.0 : 2.0)
        end
    end
    Apair = sparse(aI, 1:length(pairI), aV, n, length(pairI))

    # FVₗ = prefactor[l+1]·√Sₗ
    prefactor = if d == 3
        λ = zonal_coefficients(dims, 2 * Float64(delta))
        [λ[l+1] * sqrt((2l + 1) / (4π)) for l = 0:dims-1]
    else
        σ = Float64(delta)
        [σ / sqrt(2π) * exp(-(σ * k)^2 / 2) for k = 0:dims-1]
    end

    return HarmonicOptParams(
        d, Float64(delta), Int(dims), weights, prefactor,
        Gmem, Smem, Dload, mload, Grxn, fixmask,
        pairI, pairJ, Apair, n,
    )
end

"design positions without a solve: X(x) = X0 + reshape(Sx·x, 3, :)"
_positions(x::AbstractVector, p::OptParams) =
    p.X0 + reshape(_constmul(p.Sx, x), 3, :)

# smooth column normalization; the ε floor keeps zero columns (unloaded free
# DOFs never reach here, but exact zeros must not produce NaN directions)
_colnorms(A) = vec(sqrt.(sum(abs2, A; dims = 1) .+ 1e-18))

"unit bump directions (dimension × M) and signed magnitudes (M) of a design evaluation"
function _bump_data(res::OptResults, p::OptParams, hp::HarmonicOptParams)
    f = axial_force(res, p)                     # n_el, differentiable
    ΔX = _mulconst(res.X, transpose(p.Cinc))    # 3 × n_el, end - start
    Ehat3 = ΔX ./ transpose(res.L)              # unit element vectors

    # support reactions, matching Asap's reaction(): elastic member forces
    # accumulated at fixed translational DOFs, Rᵢ = Σₑ Cinc[e,i]·fₑ·Êₑ (per-DOF
    # masked; applied loads at supports are NOT subtracted)
    Rmat3 = hp.fixmask .* _mulconst(Ehat3 .* transpose(f), p.Cinc)

    if hp.dimension == 3
        Ehat, Rmat = Ehat3, Rmat3
    else # planar: project and renormalize (z ≈ 0 for planar models)
        Exy = Ehat3[1:2, :]
        Ehat = Exy ./ transpose(_colnorms(Exy))
        Rmat = Rmat3[1:2, :]
    end

    Rraw = _mulconst(Rmat, hp.Grxn)             # dimension × M_rxn
    mrxn = _colnorms(Rraw)
    Drxn = Rraw ./ transpose(mrxn)

    D = hcat(_mulconst(Ehat, hp.Gmem), hp.Dload, Drxn)
    m = vcat(_constmul(hp.Smem, f), hp.mload, mrxn)
    return D, m
end

"geometry-only bump data from positions alone (weights = :unit): unit member
directions, unit magnitudes — no forces, no solve"
function _unit_bump_data(X::AbstractMatrix, p::OptParams, hp::HarmonicOptParams)
    ΔX = _mulconst(X, transpose(p.Cinc))        # 3 × n_el, end - start
    E = hp.dimension == 3 ? ΔX : ΔX[1:2, :]
    Ehat = E ./ transpose(_colnorms(E))
    return _mulconst(Ehat, hp.Gmem), ones(size(hp.Gmem, 2))
end

"band sums Sₗ (dims × n_nodes): Legendre (3D) or Chebyshev cos(kΔθ) (2D) recurrence over all pairs"
function _band_sums(t, w, hp::HarmonicOptParams)
    A = hp.Apair
    S = reshape(_constmul(A, w), 1, :)          # l = 0: P₀ = T₀ = 1
    hp.dims == 1 && return S

    Pprev = one.(t)
    Pcur = t                                    # P₁ = T₁ = t
    S = vcat(S, reshape(_constmul(A, w .* t), 1, :))

    for l = 2:hp.dims-1
        Pnext = hp.dimension == 3 ?
                ((2l - 1) .* t .* Pcur .- (l - 1) .* Pprev) ./ l :
                2 .* t .* Pcur .- Pprev         # cos(kΔθ) = T_k(cos Δθ)
        S = vcat(S, reshape(_constmul(A, w .* Pnext), 1, :))
        Pprev, Pcur = Pcur, Pnext
    end

    return S
end

"descriptors from bump data: pair dots → band sums → prefactor·√Sₗ"
function _descriptor_matrix(D::AbstractMatrix, m::AbstractVector, hp::HarmonicOptParams)
    t = vec(sum(D[:, hp.pairI] .* D[:, hp.pairJ]; dims = 1))
    w = m[hp.pairI] .* m[hp.pairJ]

    S = _band_sums(t, w, hp)
    return hp.prefactor .* _band_norm.(S)
end

"""
    feature_matrix(res, p, hp)
    feature_matrix(x, p, hp)

Feature vectors of a design evaluation as a `dims × n_nodes` matrix — the
fully batched, reverse-AD-friendly layout every differentiable objective over
the descriptors should build on (see `AsapHarmonics.feature_matrix`).
`res = solve_structure(x, p)`; `hp` from `AsapHarmonics.harmonic_params`.
Under `weights = :unit` the design-vector method evaluates from positions
alone — no structural solve.
"""
feature_matrix(res::OptResults, p::OptParams, hp::HarmonicOptParams) =
    hp.weights == :unit ? _descriptor_matrix(_unit_bump_data(res.X, p, hp)..., hp) :
    _descriptor_matrix(_bump_data(res, p, hp)..., hp)

function feature_matrix(x::AbstractVector, p::OptParams, hp::HarmonicOptParams)
    hp.weights == :unit &&
        return _descriptor_matrix(_unit_bump_data(_positions(x, p), p, hp)..., hp)
    return feature_matrix(solve_structure(x, p), p, hp)
end

feature_vectors(res::OptResults, p::OptParams, hp::HarmonicOptParams) =
    collect.(eachcol(feature_matrix(res, p, hp)))

feature_vectors(x::AbstractVector, p::OptParams, hp::HarmonicOptParams) =
    collect.(eachcol(feature_matrix(x, p, hp)))

soft_complexity(res::OptResults, p::OptParams, hp::HarmonicOptParams) =
    soft_complexity(feature_matrix(res, p, hp))

"""
    soft_complexity(x, p, hp)
    soft_complexity(x, p, hp, P_or_assignments)

The differentiable design objective: smooth complexity of the design `x` —
global (3-argument), or within-cluster (4-argument, given a
`cluster_projector` or the integer assignments themselves; pass the projector
in optimization loops so it is not rebuilt every evaluation). Minimize to
reduce variation in nodal demands, overall or within each connection family.
"""
soft_complexity(x::AbstractVector, p::OptParams, hp::HarmonicOptParams) =
    soft_complexity(feature_matrix(x, p, hp))

const _ClusterSpec = Union{SparseMatrixCSC,AbstractVector{<:Integer}}

soft_complexity(res::OptResults, p::OptParams, hp::HarmonicOptParams, clusters::_ClusterSpec) =
    soft_complexity(feature_matrix(res, p, hp), clusters)

soft_complexity(x::AbstractVector, p::OptParams, hp::HarmonicOptParams, clusters::_ClusterSpec) =
    soft_complexity(feature_matrix(x, p, hp), clusters)

# exact bounding-sphere complexity: reporting only (not differentiable)
complexity(res::OptResults, p::OptParams, hp::HarmonicOptParams) =
    complexity([Float64.(fv) for fv in feature_vectors(res, p, hp)])

complexity(x::AbstractVector, p::OptParams, hp::HarmonicOptParams) =
    complexity(solve_structure(x, p), p, hp)

end # module
