# AsapHarmonics v2 test suite:
# 1. closed-form descriptors cross-validated against the sampled/transform path
# 2. invariance properties (rotation, scaling, negation, permutation)
# 3. differentiability (ForwardDiff vs central differences)
# 4. end-to-end model-facing smoke tests on solved Asap models

using AsapHarmonics
using Asap
using Test
using LinearAlgebra
using ForwardDiff

# Rodrigues rotation matrix
function rotmat(axis::Vector{Float64}, angle::Float64)
    k = normalize(axis)
    K = [0.0 -k[3] k[2]; k[3] 0.0 -k[1]; -k[2] k[1] 0.0]
    return I(3) + sin(angle) * K + (1 - cos(angle)) * K^2
end

function small_truss()
    mat = Material(200e6, 1.0, 80.0, 0.3)
    sec = Section(mat, 1e-2)
    rot = [true, true, true]
    n1 = Node([0.0, 0.0, 0.0], vcat([false, false, false], rot))
    n2 = Node([4.0, 0.0, 0.0], vcat([false, false, false], rot))
    n3 = Node([2.0, 3.0, 0.0], vcat([true, true, false], rot))
    n4 = Node([2.0, 3.0, 4.0], vcat([true, true, true], rot))
    els = AbstractElement{Float64}[
        TrussElement(n1, n3, sec), TrussElement(n2, n3, sec),
        TrussElement(n1, n4, sec), TrussElement(n2, n4, sec),
        TrussElement(n3, n4, sec)]
    loads = AbstractLoad{Float64}[
        NodeForce(n3, [0.0, -50.0, 0.0]), NodeForce(n4, [10.0, -20.0, 5.0])]
    model = Model([n1, n2, n3, n4], els, loads)
    solve!(model)
    return model
end

# a generic, non-symmetric force set on the sphere
const DIRS = normalize.([
    [1.0, 0.2, -0.1],
    [-0.3, 1.0, 0.4],
    [0.1, -0.6, 1.0],
    [-1.0, -0.8, -0.5],
    [0.0, 0.0, -1.0],
])
const MAGS = [2.0, -1.5, 0.8, 3.0, -2.2]

@testset "AsapHarmonics" begin

    @testset "spherical descriptor: closed form vs sampled transform" begin
        for delta in (5.0, 20.0)
            fv = spherical_feature_vector(DIRS, MAGS; delta = delta, dims = 12)
            F = sampled_force_function(DIRS, MAGS; delta = delta, nlat = 128)
            fv_num = spherical_feature_vector(F; dims = 12)

            @test length(fv) == 12
            @test all(isfinite, fv)
            @test isapprox(fv, fv_num; rtol = 1e-6)
        end

        # single force along +z: band energies are λₗ²(2l+1)/4π exactly
        fv1 = spherical_feature_vector([[0.0, 0.0, 1.0]], [1.0]; delta = 20, dims = 8)
        λ = zonal_coefficients(8, 40.0)
        @test isapprox(fv1, [λ[l+1] * sqrt((2l + 1) / (4π)) for l = 0:7]; rtol = 1e-12)
    end

    @testset "spherical descriptor: invariances" begin
        fv = spherical_feature_vector(DIRS, MAGS; delta = 20, dims = 10)

        R = rotmat([1.0, 2.0, 0.5], 1.234)
        fv_rot = spherical_feature_vector([R * d for d in DIRS], MAGS; delta = 20, dims = 10)
        @test isapprox(fv, fv_rot; rtol = 1e-12)

        @test isapprox(spherical_feature_vector(DIRS, 3.5 .* MAGS; delta = 20, dims = 10), 3.5 .* fv; rtol = 1e-12)
        @test isapprox(spherical_feature_vector(DIRS, -MAGS; delta = 20, dims = 10), fv; rtol = 1e-12)

        perm = [3, 1, 5, 2, 4]
        @test isapprox(spherical_feature_vector(DIRS[perm], MAGS[perm]; delta = 20, dims = 10), fv; rtol = 1e-12)

        # 3×n matrix method matches vector-of-vectors method
        @test spherical_feature_vector(reduce(hcat, DIRS), MAGS; delta = 20, dims = 10) == fv
    end

    @testset "circular descriptor: closed form vs rfft" begin
        θs = [0.2, 1.7, 3.5, 5.0]
        mags = [1.0, -2.0, 0.5, 1.5]
        dirs = [[cos(t), sin(t)] for t in θs]
        σ = 0.1
        n = 720

        fv = circular_feature_vector(θs, mags; sigma = σ, dims = 10)
        fv_dir = circular_feature_vector(dirs, mags; sigma = σ, dims = 10)
        fv_num = circular_feature_vector(circular_gaussian(dirs, mags, σ, n); dims = 10) ./ n

        @test isapprox(fv, fv_dir; rtol = 1e-12)
        @test isapprox(fv, fv_num; rtol = 1e-6)
    end

    @testset "circular descriptor: invariances" begin
        θs = [0.2, 1.7, 3.5, 5.0]
        mags = [1.0, -2.0, 0.5, 1.5]
        fv = circular_feature_vector(θs, mags; sigma = 0.1, dims = 10)

        @test isapprox(circular_feature_vector(θs .+ 0.987, mags; sigma = 0.1, dims = 10), fv; rtol = 1e-12)
        @test isapprox(circular_feature_vector(θs, 2.5 .* mags; sigma = 0.1, dims = 10), 2.5 .* fv; rtol = 1e-12)
        @test isapprox(circular_feature_vector(θs, -mags; sigma = 0.1, dims = 10), fv; rtol = 1e-12)
    end

    @testset "differentiability and type stability" begin
        f_sph(m) = sum(spherical_feature_vector(DIRS, m; delta = 20, dims = 8))
        g = ForwardDiff.gradient(f_sph, MAGS)

        h = 1e-6
        g_fd = [
            (f_sph([MAGS[1:i-1]; MAGS[i] + h; MAGS[i+1:end]]) -
             f_sph([MAGS[1:i-1]; MAGS[i] - h; MAGS[i+1:end]])) / 2h for i in eachindex(MAGS)
        ]
        @test isapprox(g, g_fd; rtol = 1e-6)

        θs = [0.2, 1.7, 3.5, 5.0]
        mags = [1.0, -2.0, 0.5, 1.5]
        f_circ(m) = sum(circular_feature_vector(θs, m; sigma = 0.1, dims = 8))
        g2 = ForwardDiff.gradient(f_circ, mags)
        g2_fd = [
            (f_circ([mags[1:i-1]; mags[i] + h; mags[i+1:end]]) -
             f_circ([mags[1:i-1]; mags[i] - h; mags[i+1:end]])) / 2h for i in eachindex(mags)
        ]
        @test isapprox(g2, g2_fd; rtol = 1e-6)

        @test (@inferred spherical_feature_vector(DIRS, MAGS; delta = 20, dims = 8)) isa Vector{Float64}
        @test (@inferred circular_feature_vector(θs, mags; sigma = 0.1, dims = 8)) isa Vector{Float64}
    end

    @testset "HarmonicAnalysis (3D end-to-end)" begin
        model = small_truss()

        ha = HarmonicAnalysis(model; dims = 8)
        @test length(ha.nodeforces) == length(model.nodes)
        @test all(fv -> length(fv) == 8 && all(isfinite, fv), ha.featurevectors)
        @test all(nf -> isempty(nf.forcefunction), ha.nodeforces)

        # loaded node's signature includes its external load
        loaded = ha.nodeforces[3]
        @test length(loaded.forcemagnitudes) >= 3 # connected members + load

        # stored sampled signatures reproduce the closed-form feature vectors
        ha_f = HarmonicAnalysis(model; dims = 8, store_functions = true, resolution = 90)
        @test all(nf -> !isempty(nf.forcefunction) && all(isfinite, nf.forcefunction), ha_f.nodeforces)
        for (ff, fv) in zip(ha_f.forcefunctions, ha_f.featurevectors)
            @test isapprox(spherical_feature_vector(ff; dims = 8), fv; rtol = 1e-6)
        end
    end

    @testset "HarmonicAnalysis2d (planar end-to-end)" begin
        mat = Material(200e6, 1.0, 80.0, 0.3)
        sec = Section(mat, 1e-2)
        rot = [true, true, true]
        p1 = Node([0.0, 0.0, 0.0], vcat([false, false, false], rot))
        p2 = Node([4.0, 0.0, 0.0], vcat([false, false, false], rot))
        p3 = Node([2.0, 3.0, 0.0], vcat([true, true, false], rot))
        pels = AbstractElement{Float64}[TrussElement(p1, p3, sec), TrussElement(p2, p3, sec)]
        pmodel = Model([p1, p2, p3], pels,
            AbstractLoad{Float64}[NodeForce(p3, [0.0, -50.0, 0.0])])
        solve!(pmodel)

        ha2 = HarmonicAnalysis2d(pmodel; dims = 8, n = 360, store_functions = true)
        @test length(ha2.nodeforces) == length(pmodel.nodes)
        @test all(fv -> all(isfinite, fv), ha2.featurevectors)

        # sampled rfft path ≈ n × closed form
        for (ff, fv) in zip(ha2.forcefunctions, ha2.featurevectors)
            @test isapprox(circular_feature_vector(ff; dims = 8) ./ 360, fv; rtol = 1e-4)
        end
    end
end
