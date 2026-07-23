# Sampled-signature path: numerical cross-validation of the closed-form
# descriptors (see descriptors.jl) and surface visualization. The analytical
# pipeline never touches this file.

"""
    spherical_gaussian(directions, magnitudes, θ, ϕ; delta = 20)

Evaluate the spherical force signature `Σᵢ fᵢ·exp(-δ‖x̂ - n̂ᵢ‖²)` at the point
with polar angle `θ` (from +z) and azimuth `ϕ`.
"""
function spherical_gaussian(
    directions::AbstractVector{<:AbstractVector},
    magnitudes::AbstractVector,
    θ::Real,
    ϕ::Real;
    delta::Real = 20,
)
    px = cos(ϕ) * sin(θ)
    py = sin(ϕ) * sin(θ)
    pz = cos(θ)

    return sum(
        f * exp(-delta * ((px - d[1])^2 + (py - d[2])^2 + (pz - d[3])^2)) for
        (d, f) in zip(directions, magnitudes)
    )
end

"""
    sampled_force_function(directions, magnitudes; delta = 20, nlat = 91)

Sample the spherical force signature on the `nlat × (2nlat - 1)` grid returned
by `FastSphericalHarmonics.sph_points(nlat)` — the grid `sph_transform`
expects (Gauss-type latitudes without the poles; azimuths excluding 2π). Do
NOT feed `sph_transform` a signature sampled on an endpoint-inclusive θ/ϕ
grid.
"""
function sampled_force_function(
    directions::AbstractVector{<:AbstractVector},
    magnitudes::AbstractVector;
    delta::Real = 20,
    nlat::Integer = 91,
)
    Θ, Φ = sph_points(nlat)
    return [spherical_gaussian(directions, magnitudes, θ, ϕ; delta = delta) for θ in Θ, ϕ in Φ]
end

"""
    spherical_feature_vector(forcefunction::AbstractMatrix; dims = 16)

Numeric-path descriptor for cross-validation: real spherical-harmonic
transform of a sampled signature (`sampled_force_function`), then the ℓ² norm
of the orthonormal coefficients of each degree band — the discrete estimate of
the L² band norm that `spherical_feature_vector(directions, magnitudes)`
computes in closed form.
"""
function spherical_feature_vector(forcefunction::AbstractMatrix{Float64}; dims::Integer = 16)
    coefficients = sph_transform(forcefunction)
    return [norm(coefficients[sph_mode(l, m)] for m = -l:l) for l = 0:dims-1]
end

"""
    circular_gaussian(positions, magnitudes, σ, n = 90)

Sample the planar force signature `Σᵢ fᵢ·exp(-θᵢ(x̂)²/(2σ²))` (θᵢ = angle
between the sample direction and force direction `positions[i]`) at `n`
equally spaced angles on `[0, 2π)`.
"""
function circular_gaussian(
    positions::AbstractVector{<:AbstractVector},
    magnitudes::AbstractVector,
    σ::Real,
    n::Integer = 90,
)
    angles = range(0, 2π, n + 1)[1:end-1]

    T = promote_type(eltype(eltype(positions)), eltype(magnitudes), typeof(float(σ)))
    force_function = zeros(T, n)

    for i = 1:n
        s, c = sincos(angles[i])

        for (position, magnitude) in zip(positions, magnitudes)
            θ = acos(clamp(position[1] * c + position[2] * s, -1, 1))
            force_function[i] += exp(-(θ / σ)^2 / 2) * magnitude
        end
    end

    return force_function
end

"""
    circular_feature_vector(forcefunction::AbstractVector; dims = 16)

Numeric-path descriptor for cross-validation: magnitudes of the first `dims`
real-FFT coefficients of a sampled signature. Relative to the closed-form
`circular_feature_vector(angles, magnitudes)` this is scaled by the sample
count: `rfft`-based values ≈ `n·FVₖ`.
"""
function circular_feature_vector(forcefunction::AbstractVector{Float64}; dims::Integer = 16)
    return abs.(rfft(forcefunction)[1:dims])
end

# visualization helpers: Cartesian coordinates of a (θ, ϕ) grid, for plotting
# force signatures as radially-scaled sphere surfaces.
make_xsphere(thetas::AbstractVector, phis::AbstractVector) = [cos(ϕ) * sin(θ) for θ in thetas, ϕ in phis]
make_ysphere(thetas::AbstractVector, phis::AbstractVector) = [sin(ϕ) * sin(θ) for θ in thetas, ϕ in phis]
make_zsphere(thetas::AbstractVector, phis::AbstractVector) = [cos(θ) for θ in thetas, ϕ in phis]
