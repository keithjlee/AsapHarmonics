# Analytical (closed-form) shape descriptors of Gaussian-bump force signatures.
#
# The nodal force signature f(x̂) = Σᵢ fᵢ·exp(-δ‖x̂ - n̂ᵢ‖²) uses squared chord
# distance on the unit sphere, so each bump is exactly a von Mises–Fisher zonal
# kernel e^(-κ)·e^(κ x̂·n̂ᵢ) with κ = 2δ. Zonal kernels expand in closed form
# (Funk–Hecke theorem), which collapses the rotation-invariant feature vector to
# pairwise Legendre sums over force directions — no grids, no transforms.

"""
    zonal_coefficients(dims::Integer, κ::Real) -> Vector{Float64}

Legendre coefficients `λₗ`, `l = 0:dims-1`, of the unit-height spherical
Gaussian (von Mises–Fisher) kernel `g(t) = e^(κ(t-1))`, `t = x̂·n̂`:

    g(x̂·n̂) = Σₗ λₗ·(2l+1)/(4π)·Pₗ(x̂·n̂),    λₗ = 4π·e^(-κ)·iₗ(κ)

where `iₗ` is the modified spherical Bessel function of the first kind, via the
Gegenbauer expansion `e^(κt) = Σₗ (2l+1)·iₗ(κ)·Pₗ(t)` (DLMF §10.60.7). The
scaled product `e^(-κ)·iₗ(κ)` is evaluated overflow-free as
`√(π/2κ)·besselix(l+1/2, κ)`.

`κ` is a constant hyperparameter of the descriptor: it is evaluated in
`Float64` and gradients are never taken through it.
"""
function zonal_coefficients(dims::Integer, κ::Real)
    κ64 = Float64(κ)
    iszero(κ64) && return [l == 0 ? 4π : 0.0 for l = 0:dims-1]

    scale = 4π * sqrt(π / (2κ64))
    return [scale * besselix(l + 0.5, κ64) for l = 0:dims-1]
end

"""
    pairwise_legendre_sums(directions, magnitudes, dims::Integer) -> Vector

`Sₗ = Σᵢⱼ fᵢ·fⱼ·Pₗ(n̂ᵢ·n̂ⱼ)` for `l = 0:dims-1`, accumulated for all degrees in
a single pass over force pairs using the Legendre three-term recurrence.
Generic over `Real` element types (ForwardDiff duals flow through).
"""
function pairwise_legendre_sums(
    directions::AbstractVector{<:AbstractVector},
    magnitudes::AbstractVector,
    dims::Integer,
)
    T = promote_type(eltype(eltype(directions)), eltype(magnitudes))
    S = zeros(T, dims)
    n = length(magnitudes)

    @inbounds for i = 1:n
        dᵢ = directions[i]
        fᵢ = magnitudes[i]

        for j = i:n
            # off-diagonal pairs appear twice in the full double sum
            w = (i == j ? one(T) : T(2)) * fᵢ * magnitudes[j]
            t = dot(dᵢ, directions[j])

            S[1] += w
            dims == 1 && continue

            Pprev = one(T) # P₀
            Pcurr = t      # P₁
            S[2] += w * t

            for l = 2:dims-1
                Pnext = ((2l - 1) * t * Pcurr - (l - 1) * Pprev) / l
                S[l+1] += w * Pnext
                Pprev, Pcurr = Pcurr, Pnext
            end
        end
    end

    return S
end

"""
    spherical_feature_vector(directions, magnitudes; delta = 20, dims = 16)

Rotation-invariant spherical-harmonic shape descriptor of the nodal force
signature `f(x̂) = Σᵢ fᵢ·exp(-δ‖x̂ - n̂ᵢ‖²)`, computed in closed form.

`directions` are unit vectors on S² (a vector of 3-vectors, or a 3×n matrix
whose columns are directions); `magnitudes` are the associated signed force
values. The component `FVₗ` is the L² norm of the degree-`l` band of the
signature's spherical-harmonic expansion:

    FVₗ = λₗ·√( (2l+1)/(4π) · Σᵢⱼ fᵢ·fⱼ·Pₗ(n̂ᵢ·n̂ⱼ) ),    λₗ = 4π·e^(-2δ)·iₗ(2δ)

obtained from the Funk–Hecke theorem and the spherical-harmonic addition
theorem (the chord-distance Gaussian is the von Mises–Fisher kernel with
κ = 2δ). Cost is O(n²·dims); the result is exactly invariant to rotations of
the direction set, satisfies `FV(c·f) = |c|·FV(f)` and `FV(-f) = FV(f)`, and is
smooth in both magnitudes and directions.

References: Kazhdan, Funkhouser & Rusinkiewicz (SGP 2003); Lee, Danhaive &
Mueller (2022); Funk–Hecke: Müller, *Spherical Harmonics* (1966), Atkinson &
Han (2012); von Mises–Fisher kernel: Fisher (1953), Mardia & Jupp (2000);
Bessel expansion: DLMF §10.60.7.
"""
# NaN-safe root of a band energy. Sₗ ≥ 0 in exact arithmetic, but roundoff
# can leave it slightly negative — in particular the l = 1 (k = 1) band of a
# COMPLETE nodal signature is the squared resultant force, which is ~0 at
# every equilibrated node. A branch (not a clamp) keeps AD clean: clamping to
# exactly 0.0 feeds sqrt a zero with nonzero perturbation → 0/0 = NaN partials
# under ForwardDiff, while the zero branch has the correct zero derivative.
_band_norm(s) = s > 0 ? sqrt(s) : zero(s)

function spherical_feature_vector(
    directions::AbstractVector{<:AbstractVector},
    magnitudes::AbstractVector;
    delta::Real = 20,
    dims::Integer = 16,
)
    length(directions) == length(magnitudes) ||
        throw(DimensionMismatch("$(length(directions)) directions for $(length(magnitudes)) magnitudes"))
    dims ≥ 1 || throw(ArgumentError("dims must be ≥ 1"))

    λ = zonal_coefficients(dims, 2delta)
    S = pairwise_legendre_sums(directions, magnitudes, dims)

    return [λ[l+1] * _band_norm((2l + 1) / (4π) * S[l+1]) for l = 0:dims-1]
end

function spherical_feature_vector(directions::AbstractMatrix, magnitudes::AbstractVector; kwargs...)
    size(directions, 1) == 3 ||
        throw(DimensionMismatch("expected a 3×n matrix of direction columns, got $(size(directions))"))
    return spherical_feature_vector(collect(eachcol(directions)), magnitudes; kwargs...)
end

"""
    circular_feature_vector(angles, magnitudes; sigma = 0.1, dims = 16)
    circular_feature_vector(directions, magnitudes; sigma = 0.1, dims = 16)

Rotation-invariant Fourier shape descriptor of the planar force signature
`s(θ) = Σᵢ fᵢ·exp(-(θ - θᵢ)²/(2σ²))` (geodesic angle on the circle), computed
in closed form. Force directions may be given as angles `θᵢ` or unit 2-vectors.

The component `FVₖ` is the magnitude of the signature's k-th complex Fourier
coefficient; a Gaussian's Fourier coefficients are `σ/√(2π)·e^(-σ²k²/2)`
(wrapping/truncation corrections are O(e^(-π²/2σ²)), negligible for σ ≪ π):

    FVₖ = σ/√(2π) · e^(-σ²k²/2) · √( Σᵢⱼ fᵢ·fⱼ·cos(k(θᵢ - θⱼ)) )

Exactly invariant to rotation of the direction set (a rotation only shifts
Fourier phases); `FV(c·f) = |c|·FV(f)`; `FV(-f) = FV(f)`; O(n·dims).
A raw `rfft` of the n-point sampled signature returns ≈ `n·FVₖ`.

References: Zahn & Roskies (1972); wrapped-Gaussian coefficients: Mardia &
Jupp, *Directional Statistics* (2000).
"""
function circular_feature_vector(
    angles::AbstractVector{<:Real},
    magnitudes::AbstractVector;
    sigma::Real = 0.1,
    dims::Integer = 16,
)
    length(angles) == length(magnitudes) ||
        throw(DimensionMismatch("$(length(angles)) angles for $(length(magnitudes)) magnitudes"))
    dims ≥ 1 || throw(ArgumentError("dims must be ≥ 1"))

    T = promote_type(eltype(angles), eltype(magnitudes))

    return map(0:dims-1) do k
        a = zero(T)
        b = zero(T)
        for (θ, f) in zip(angles, magnitudes)
            s, c = sincos(k * θ)
            a += f * c
            b += f * s
        end
        # _band_norm: the k = 1 band of a complete signature is the resultant
        # force, ~0 at equilibrated nodes (see the 3D note above)
        sigma / sqrt(2π) * exp(-(sigma * k)^2 / 2) * _band_norm(a^2 + b^2)
    end
end

function circular_feature_vector(
    directions::AbstractVector{<:AbstractVector},
    magnitudes::AbstractVector;
    kwargs...,
)
    return circular_feature_vector([atan(d[2], d[1]) for d in directions], magnitudes; kwargs...)
end
