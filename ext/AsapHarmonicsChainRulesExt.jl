# Reverse-mode AD rules. The Legendre kernel mutates a local accumulator (fast
# and ForwardDiff/Enzyme-safe), so ChainRules-based engines (Zygote) get an
# analytic pullback instead:
#
#   Sₗ = Σᵢⱼ fᵢ fⱼ Pₗ(tᵢⱼ),  tᵢⱼ = dᵢ·dⱼ
#   ∂Sₗ/∂fₖ = 2 Σⱼ fⱼ Pₗ(tₖⱼ)
#   ∂Sₗ/∂dₖ = 2 fₖ Σⱼ fⱼ Pₗ′(tₖⱼ) dⱼ      (j = k term included: tₖₖ = dₖ·dₖ)
#
# with Pₗ′ by the derivative recurrence Pₗ′ = Pₗ₋₂′ + (2l-1)·Pₗ₋₁.

module AsapHarmonicsChainRulesExt

using AsapHarmonics, ChainRulesCore, LinearAlgebra
import AsapHarmonics: pairwise_legendre_sums, zonal_coefficients, _constmul, _mulconst

# constant hyperparameter data by design (documented): κ is never a
# differentiation target, and besselix must not be traced
ChainRulesCore.@non_differentiable zonal_coefficients(::Any, ::Any)

# products with constant selection/aggregation/incidence operators: no
# cotangent for the constant matrix (the generic mul rrule would materialize
# a dense outer product for it)
function ChainRulesCore.rrule(::typeof(_constmul), A::AbstractMatrix, x)
    _constmul_pullback(Δ) = (NoTangent(), NoTangent(), A' * unthunk(Δ))
    return A * x, _constmul_pullback
end

function ChainRulesCore.rrule(::typeof(_mulconst), X, A::AbstractVecOrMat)
    _mulconst_pullback(Δ) = (NoTangent(), unthunk(Δ) * A', NoTangent())
    return X * A, _mulconst_pullback
end

function ChainRulesCore.rrule(
    ::typeof(pairwise_legendre_sums),
    directions::AbstractVector{<:AbstractVector},
    magnitudes::AbstractVector,
    dims::Integer,
)
    S = pairwise_legendre_sums(directions, magnitudes, dims)
    project_dirs = ProjectTo(directions)
    project_mags = ProjectTo(magnitudes)

    function pairwise_legendre_sums_pullback(ΔS_raw)
        ΔS = unthunk(ΔS_raw)
        T = eltype(S)
        n = length(magnitudes)

        f̄ = zeros(T, n)
        d̄ = [zeros(T, length(d)) for d in directions]

        for k = 1:n, j = 1:n
            t = dot(directions[k], directions[j])

            # A = Σₗ ΔSₗ Pₗ(t), B = Σₗ ΔSₗ Pₗ′(t) in one recurrence pass
            A = ΔS[1] # P₀ = 1, P₀′ = 0
            B = zero(T)
            if dims > 1
                Pp, Pc = one(t), t   # P₀, P₁
                Dp, Dc = zero(t), one(t) # P₀′, P₁′
                A += ΔS[2] * Pc
                B += ΔS[2] * Dc
                for l = 2:dims-1
                    Pn = ((2l - 1) * t * Pc - (l - 1) * Pp) / l
                    Dn = Dp + (2l - 1) * Pc
                    A += ΔS[l+1] * Pn
                    B += ΔS[l+1] * Dn
                    Pp, Pc = Pc, Pn
                    Dp, Dc = Dc, Dn
                end
            end

            f̄[k] += 2 * magnitudes[j] * A
            d̄[k] .+= (2 * magnitudes[k] * magnitudes[j] * B) .* directions[j]
        end

        return (NoTangent(), project_dirs(d̄), project_mags(f̄), NoTangent())
    end

    return S, pairwise_legendre_sums_pullback
end

end # module
