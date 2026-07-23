# Model-facing layer: extract per-node force signatures from a solved Asap
# model and compute their closed-form feature vectors. The 2D (circle) and 3D
# (sphere) analyses share one implementation, parameterized by D.
#
# Sign convention (uniform across D):
#   - incident members: signed axial force (tension +, compression -) at the
#     outward unit member direction
#   - applied point forces and support reactions: magnitude |P| at direction
#     P/|P| (equivalent to a tension member pointing along the force)
# This convention is equivariant: rotating the whole problem (geometry and
# loads together) rotates every signature rigidly, so feature vectors are
# unchanged.

"""
    NodeSignature{D,T}

Force signature of a single node on the circle (`D = 2`) or sphere (`D = 3`):
unit `directions` and signed `magnitudes` of every force acting at the node —
incident member axial forces, applied `NodeForce`s, and the support reaction.
`index`/`id` identify the node in its model.

Built from a solved model via
`NodeSignature{D}(node, model, C, forces; normalize_forces = false)` where `C`
is `connectivity(model)` and `forces` the member axial forces;
`normalize_forces = true` replaces every magnitude by its sign. Usually
constructed in bulk by [`HarmonicAnalysis`](@ref) / [`HarmonicAnalysis2d`](@ref).
"""
struct NodeSignature{D,T}
    index::Int
    id::Symbol
    directions::Vector{SVector{D,T}}
    magnitudes::Vector{T}
end

_project(v::SVector{3,T}, ::Val{3}) where {T} = v
_project(v::SVector{3,T}, ::Val{2}) where {T} = SVector{2,T}(v[1], v[2])

function NodeSignature{D}(
    node::Node,
    model::Model,
    C::SparseMatrixCSC{<:Integer},
    forces::AbstractVector;
    normalize_forces::Bool = false,
) where {D}
    D == 2 || D == 3 || throw(ArgumentError("NodeSignature dimension must be 2 or 3, got $D"))

    i = node.index
    T = promote_type(eltype(node.position), eltype(forces))

    i_connected = findall(!iszero, C[:, i])
    factors = Vector(C[:, i][i_connected])

    directions = Vector{SVector{D,T}}()
    magnitudes = Vector{T}()

    # incident members: signed axial force at the outward direction
    for (e, factor) in zip(model.elements[i_connected], factors)
        d = -normalize(e.nodeEnd.position - e.nodeStart.position) * factor
        push!(directions, normalize(_project(SVector{3,T}(d), Val(D))))
    end
    append!(magnitudes, forces[i_connected])

    # applied point forces: |P| at P̂
    for load in model.loads
        (load isa NodeForce && load.node.index == i) || continue

        P = _project(SVector{3,T}(load.value), Val(D))
        m = norm(P)
        iszero(m) && continue

        push!(directions, P / m)
        push!(magnitudes, m)
    end

    # support reaction (translational components): |R| at R̂
    r = reaction(model.results, node)
    R = _project(SVector{3,T}(r[1], r[2], r[3]), Val(D))
    mR = norm(R)
    if !iszero(mR)
        push!(directions, R / mR)
        push!(magnitudes, mR)
    end

    normalize_forces && (magnitudes = sign.(magnitudes))

    return NodeSignature{D,T}(i, node.id, directions, magnitudes)
end

"""
    feature_vector(signature::NodeSignature; delta, dims = 16)

Closed-form rotation-invariant feature vector of a node signature — dispatches
to [`spherical_feature_vector`](@ref) (`D = 3`, `delta` defaults to 20) or
[`circular_feature_vector`](@ref) (`D = 2`, where `delta` is the angular
Gaussian width σ, defaulting to 0.1).
"""
feature_vector(sig::NodeSignature{3}; delta::Real = 20, dims::Integer = 16) =
    spherical_feature_vector(sig.directions, sig.magnitudes; delta = delta, dims = dims)

feature_vector(sig::NodeSignature{2}; delta::Real = 0.1, dims::Integer = 16) =
    circular_feature_vector(sig.directions, sig.magnitudes; sigma = delta, dims = dims)

# sampled signatures (visualization / numeric cross-validation), on demand
sampled_force_function(sig::NodeSignature{3}; delta::Real = 20, nlat::Integer = 91) =
    sampled_force_function(sig.directions, sig.magnitudes; delta = delta, nlat = nlat)

circular_gaussian(sig::NodeSignature{2}, σ::Real = 0.1, n::Integer = 90) =
    circular_gaussian(sig.directions, sig.magnitudes, σ, n)

"""
    HarmonicAnalysis{D,T}

Shape-descriptor analysis of every node of a solved model: `signatures` (one
[`NodeSignature{D,T}`](@ref) per node, in model order) and their
rotation-invariant `featurevectors` (each of length `dims`, degrees/frequencies
`0:dims-1`), plus the kernel parameter `delta` used.

Construct with [`HarmonicAnalysis(model)`](@ref) (3D, spherical harmonics) or
[`HarmonicAnalysis2d(model)`](@ref) (planar, Fourier).
"""
struct HarmonicAnalysis{D,T}
    signatures::Vector{NodeSignature{D,T}}
    featurevectors::Vector{Vector{T}}
    delta::Float64
    dims::Int
end

function _harmonic_analysis(
    ::Val{D},
    model::Model,
    delta::Real,
    dims::Integer,
    normalize_forces::Bool,
) where {D}
    C = connectivity(model)
    forces = [axial_force(model.results, el) for el in model.elements]

    signatures = [
        NodeSignature{D}(node, model, C, forces; normalize_forces = normalize_forces) for
        node in model.nodes
    ]
    featurevectors = [feature_vector(sig; delta = delta, dims = dims) for sig in signatures]

    return HarmonicAnalysis(signatures, featurevectors, Float64(delta), Int(dims))
end

"""
    HarmonicAnalysis(model; delta = 20, dims = 16, normalize_forces = false)

Spherical-harmonic shape-descriptor analysis of every node of a solved 3D
`Model`. `delta` is the sharpness of the Gaussian force bumps (each bump
approaches a Dirac spike as `delta → ∞`); `dims` is the feature-vector length
(spherical-harmonic degrees `l = 0:dims-1`); `normalize_forces = true`
discards force magnitudes and keeps only orientations and tension/compression
signs.
"""
HarmonicAnalysis(model::Model; delta::Real = 20, dims::Integer = 16, normalize_forces::Bool = false) =
    _harmonic_analysis(Val(3), model, delta, dims, normalize_forces)

"""
    HarmonicAnalysis2d(model; delta = 0.1, dims = 16, normalize_forces = false)

Fourier shape-descriptor analysis of every node of a solved **planar (XY)**
`Model`, returning a `HarmonicAnalysis{2}`. `delta` is the angular width σ of
the Gaussian force bumps; `dims` is the feature-vector length (frequencies
`k = 0:dims-1`).
"""
HarmonicAnalysis2d(model::Model; delta::Real = 0.1, dims::Integer = 16, normalize_forces::Bool = false) =
    _harmonic_analysis(Val(2), model, delta, dims, normalize_forces)

function Base.show(io::IO, sig::NodeSignature{D}) where {D}
    print(io, "NodeSignature{$D}(node $(sig.index) [:$(sig.id)], $(length(sig.magnitudes)) forces)")
end

function Base.show(io::IO, ha::HarmonicAnalysis{D}) where {D}
    print(io, "HarmonicAnalysis{$D}: $(length(ha.signatures)) nodes, dims = $(ha.dims), delta = $(ha.delta)")
end
