# Differentiable shape-descriptor pipeline over AsapOptim design evaluations:
# design vector x -> solve_structure -> per-node signatures -> closed-form
# feature vectors -> (soft) complexity. All static gather information is
# precompiled into HarmonicOptParams; the per-evaluation path uses only
# AD-transparent operations (batched incidence matmuls, broadcasts, and the
# closed-form descriptors, whose reverse rule lives in the ChainRules
# extension).

module AsapHarmonicsOptimExt

using AsapHarmonics, AsapOptim, Asap
using LinearAlgebra, SparseArrays
import AsapHarmonics: harmonic_params, feature_vectors, soft_complexity, complexity

"""
Constant data of a shape-descriptor evaluation over an `OptParams` — see
`AsapHarmonics.harmonic_params`.
"""
struct HarmonicOptParams
    dimension::Int                        # 2 (circle) or 3 (sphere)
    delta::Float64                        # kernel parameter (δ in 3D, σ in 2D)
    dims::Int                             # feature-vector length
    elids::Vector{Vector{Int}}            # incident elements per node
    signs::Vector{Vector{Float64}}        # Cinc[e, i] per incident element
    loaddirs::Vector{Vector{Vector{Float64}}} # unit direction per applied load
    loadmags::Vector{Vector{Float64}}     # magnitude per applied load
    issupport::Vector{Bool}               # node has a fixed translational DOF
    Pmat::Matrix{Float64}                 # 3 × n_nodes total applied loads
end

function harmonic_params(p::OptParams; delta::Real = 20, dims::Integer = 16, dimension::Integer = 3)
    dimension == 2 || dimension == 3 ||
        throw(ArgumentError("dimension must be 2 or 3, got $dimension"))

    model = p.model
    n = length(model.nodes)

    elids = [Int[] for _ = 1:n]
    signs = [Float64[] for _ = 1:n]
    rows = rowvals(p.Cinc)
    vals = nonzeros(p.Cinc)
    for i = 1:n, k in nzrange(p.Cinc, i)
        push!(elids[i], rows[k])
        push!(signs[i], vals[k])
    end

    loaddirs = [Vector{Float64}[] for _ = 1:n]
    loadmags = [Float64[] for _ = 1:n]
    Pmat = zeros(3, n)
    for load in model.loads
        load isa NodeForce || continue
        i = load.node.index
        P = collect(Float64, load.value)
        Pmat[:, i] .+= P
        m = norm(P)
        iszero(m) && continue
        push!(loaddirs[i], P ./ m)
        push!(loadmags[i], m)
    end

    # Asap fixity: true = free; a node with any fixed translational DOF reacts
    issupport = [any(.!node.fixity[1:3]) for node in model.nodes]

    return HarmonicOptParams(
        Int(dimension), Float64(delta), Int(dims),
        elids, signs, loaddirs, loadmags, issupport, Pmat,
    )
end

_signature_fv(dirs, mags, hp::HarmonicOptParams) =
    hp.dimension == 3 ?
    spherical_feature_vector(dirs, mags; delta = hp.delta, dims = hp.dims) :
    circular_feature_vector([normalize(d[1:2]) for d in dirs], mags; sigma = hp.delta, dims = hp.dims)

function feature_vectors(res::OptResults, p::OptParams, hp::HarmonicOptParams)
    f = axial_force(res, p)                       # n_el, differentiable
    ΔX = res.X * transpose(p.Cinc)                # 3 × n_el, end - start
    Ehat = ΔX ./ transpose(res.L)                 # unit element vectors

    # support reactions from nodal equilibrium: Rᵢ = Σₑ Cinc[e,i]·fₑ·Êₑ - Pᵢ
    # (at free nodes this is the solver residual ≈ 0 and is masked off)
    Rmat = (Ehat .* transpose(f)) * p.Cinc .- hp.Pmat

    T = promote_type(eltype(f), Float64)

    return map(1:length(hp.elids)) do i
        # members: signed axial force at the outward direction -Cinc[e,i]·Êₑ
        dirs = [(-hp.signs[i][k]) .* Ehat[:, e] for (k, e) in enumerate(hp.elids[i])]
        mags = f[hp.elids[i]]

        if !isempty(hp.loadmags[i])
            dirs = vcat(dirs, [T.(d) for d in hp.loaddirs[i]])
            mags = vcat(mags, hp.loadmags[i])
        end

        if hp.issupport[i]
            R = Rmat[:, i]
            mR = norm(R)
            if mR > 0
                dirs = vcat(dirs, [R ./ mR])
                mags = vcat(mags, mR)
            end
        end

        _signature_fv(dirs, mags, hp)
    end
end

feature_vectors(x::AbstractVector, p::OptParams, hp::HarmonicOptParams) =
    feature_vectors(solve_structure(x, p), p, hp)

soft_complexity(res::OptResults, p::OptParams, hp::HarmonicOptParams) =
    soft_complexity(feature_vectors(res, p, hp))

"""
    soft_complexity(x, p, hp)

The differentiable design objective: smooth complexity of the design `x`
evaluated through `solve_structure`. Minimize to reduce variation in nodal
force demands.
"""
soft_complexity(x::AbstractVector, p::OptParams, hp::HarmonicOptParams) =
    soft_complexity(feature_vectors(x, p, hp))

# exact bounding-sphere complexity: reporting only (not differentiable)
complexity(res::OptResults, p::OptParams, hp::HarmonicOptParams) =
    complexity([Float64.(fv) for fv in feature_vectors(res, p, hp)])

complexity(x::AbstractVector, p::OptParams, hp::HarmonicOptParams) =
    complexity(solve_structure(x, p), p, hp)

end # module
