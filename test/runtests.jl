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
        @test ha isa HarmonicAnalysis{3,Float64}
        @test length(ha.signatures) == length(model.nodes)
        @test all(fv -> length(fv) == 8 && all(isfinite, fv), ha.featurevectors)

        # loaded node's signature includes its external load
        loaded = ha.signatures[3]
        @test length(loaded.magnitudes) >= 3 # connected members + load
        @test all(d -> norm(d) ≈ 1, loaded.directions)

        # on-demand sampled signatures reproduce the closed-form feature vectors
        for (sig, fv) in zip(ha.signatures, ha.featurevectors)
            F = sampled_force_function(sig; nlat = 91)
            @test isapprox(spherical_feature_vector(F; dims = 8), fv; rtol = 1e-6)
        end
    end

    function planar_truss(α::Real) # in-plane rotation of geometry and loads by α
        R = [cos(α) -sin(α) 0.0; sin(α) cos(α) 0.0; 0.0 0.0 1.0]
        mat = Material(200e6, 1.0, 80.0, 0.3)
        sec = Section(mat, 1e-2)
        rot = [true, true, true]
        p1 = Node(R * [0.0, 0.0, 0.0], vcat([false, false, false], rot))
        p2 = Node(R * [4.0, 0.0, 0.0], vcat([false, false, false], rot))
        p3 = Node(R * [2.0, 3.0, 0.0], vcat([true, true, false], rot))
        # loads at the supports keep every signature non-degenerate: with a
        # bare pinned 2-bar truss each reaction is exactly collinear with its
        # member and the support signatures legitimately cancel to zero
        pels = AbstractElement{Float64}[
            TrussElement(p1, p3, sec), TrussElement(p2, p3, sec), TrussElement(p1, p2, sec)]
        pmodel = Model([p1, p2, p3], pels,
            AbstractLoad{Float64}[
                NodeForce(p3, R * [0.0, -50.0, 0.0]),
                NodeForce(p1, R * [7.0, 3.0, 0.0]),
                NodeForce(p2, R * [-4.0, 6.0, 0.0])])
        solve!(pmodel)
        return pmodel
    end

    @testset "HarmonicAnalysis2d (planar end-to-end)" begin
        pmodel = planar_truss(0.0)

        ha2 = HarmonicAnalysis2d(pmodel; dims = 8)
        @test ha2 isa HarmonicAnalysis{2,Float64}
        @test length(ha2.signatures) == length(pmodel.nodes)
        @test all(fv -> all(isfinite, fv), ha2.featurevectors)

        # on-demand sampled rfft path ≈ n × closed form
        for (sig, fv) in zip(ha2.signatures, ha2.featurevectors)
            ff = circular_gaussian(sig, 0.1, 360)
            @test isapprox(circular_feature_vector(ff; dims = 8) ./ 360, fv; rtol = 1e-4, atol = 1e-9)
        end
        @test all(fv -> any(>(1e-3), fv), ha2.featurevectors) # non-degenerate signatures
    end

    @testset "whole-model rotation equivariance" begin
        # rotating geometry AND loads together must leave all descriptors
        # unchanged (the v1 2D load/reaction sign convention depended on the
        # global vertical and violated this)
        ha_ref = HarmonicAnalysis2d(planar_truss(0.0); dims = 8)
        for α in (0.6, π / 2, π, 4.0)
            ha_rot = HarmonicAnalysis2d(planar_truss(α); dims = 8)
            for (fv, fv_rot) in zip(ha_ref.featurevectors, ha_rot.featurevectors)
                @test isapprox(fv, fv_rot; rtol = 1e-8, atol = 1e-9)
            end
        end

        # 3D: rotate the full model about z (supports stay axis-consistent)
        function rotated_truss(α)
            R = [cos(α) -sin(α) 0.0; sin(α) cos(α) 0.0; 0.0 0.0 1.0]
            mat = Material(200e6, 1.0, 80.0, 0.3)
            sec = Section(mat, 1e-2)
            rot = [true, true, true]
            n1 = Node(R * [0.0, 0.0, 0.0], vcat([false, false, false], rot))
            n2 = Node(R * [4.0, 0.0, 0.0], vcat([false, false, false], rot))
            n3 = Node(R * [2.0, 3.0, 0.0], vcat([true, true, false], rot))
            n4 = Node(R * [2.0, 3.0, 4.0], vcat([true, true, true], rot))
            els = AbstractElement{Float64}[
                TrussElement(n1, n3, sec), TrussElement(n2, n3, sec),
                TrussElement(n1, n4, sec), TrussElement(n2, n4, sec),
                TrussElement(n3, n4, sec)]
            loads = AbstractLoad{Float64}[
                NodeForce(n3, R * [0.0, -50.0, 0.0]), NodeForce(n4, R * [10.0, -20.0, 5.0])]
            model = Model([n1, n2, n3, n4], els, loads)
            solve!(model)
            return model
        end

        ha3_ref = HarmonicAnalysis(rotated_truss(0.0); dims = 8)
        ha3_rot = HarmonicAnalysis(rotated_truss(1.1); dims = 8)
        for (fv, fv_rot) in zip(ha3_ref.featurevectors, ha3_rot.featurevectors)
            @test isapprox(fv, fv_rot; rtol = 1e-8, atol = 1e-9)
        end
    end

    @testset "bounding sphere (Welzl)" begin
        # two points: radius = half the distance, center = midpoint
        b = bounding_sphere([[0.0, 0.0], [4.0, 0.0]])
        @test b.radius ≈ 2.0 && b.center ≈ [2.0, 0.0]

        # equilateral triangle, side s: circumradius s/√3
        s = 2.0
        tri = [[0.0, 0.0], [s, 0.0], [s / 2, s * √3 / 2]]
        @test bounding_sphere(tri).radius ≈ s / √3

        # obtuse triangle: minimal ball is on the longest side, NOT the
        # circumcircle (catches naive implementations)
        obtuse = [[0.0, 0.0], [10.0, 0.0], [5.0, 0.5]]
        bo = bounding_sphere(obtuse)
        @test bo.radius ≈ 5.0 && bo.center ≈ [5.0, 0.0]

        # unit square: radius = half diagonal
        sq = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
        @test bounding_sphere(sq).radius ≈ √2 / 2

        # high-dimensional sanity (deterministic pseudo-random points)
        pts = [[sin(17.3i + 3.1j) + 0.1j for j = 1:17] for i = 1:60]
        bh = bounding_sphere(pts)
        @test all(p -> norm(p - bh.center) <= bh.radius * (1 + 1e-8), pts)
        maxpair = maximum(norm(pts[i] - pts[j]) for i = 1:60 for j = i+1:60)
        centroid = sum(pts) / 60
        @test bh.radius >= maxpair / 2 - 1e-8
        @test bh.radius <= maximum(norm(p - centroid) for p in pts) + 1e-8
        # at least two points on the boundary of the minimal ball
        @test count(p -> norm(p - bh.center) >= bh.radius * (1 - 1e-6), pts) >= 2
    end

    @testset "analysis layer" begin
        ha = HarmonicAnalysis(small_truss(); dims = 8)
        n = length(ha.signatures)

        F = feature_matrix(ha)
        @test size(F) == (8, n)

        D = distance_matrix(ha)
        @test size(D) == (n, n) && issymmetric(D) && all(iszero, D[i, i] for i = 1:n)
        @test D[1, 2] ≈ norm(ha.featurevectors[1] - ha.featurevectors[2])

        c = complexity(ha)
        sc = soft_complexity(ha)
        @test c > 0 && sc > 0
        @test sc <= 2c + 1e-12

        # complexity scales linearly with force magnitudes: FV(cf) = c FV(f)
        fvs2 = [2.0 .* fv for fv in ha.featurevectors]
        @test complexity(fvs2) ≈ 2c
        @test soft_complexity(fvs2) ≈ 2sc

        # differentiable chain: forces -> feature vectors -> soft complexity
        dirs_b = normalize.([[1.0, 0.0, 0.1], [0.0, -1.0, 0.4], [-0.7, 0.7, 0.0]])
        function chain(m)
            fv1 = spherical_feature_vector(DIRS, m; delta = 20, dims = 8)
            fv2 = spherical_feature_vector(dirs_b, m[1:3]; delta = 20, dims = 8)
            return soft_complexity([fv1, fv2])
        end
        g = ForwardDiff.gradient(chain, MAGS)
        h = 1e-6
        g_fd = [
            (chain([MAGS[1:i-1]; MAGS[i] + h; MAGS[i+1:end]]) -
             chain([MAGS[1:i-1]; MAGS[i] - h; MAGS[i+1:end]])) / 2h for i in eachindex(MAGS)
        ]
        @test isapprox(g, g_fd; rtol = 1e-5)

        # per-cluster complexity: cluster of identical demands has radius 0
        fvs = [ha.featurevectors; [copy(ha.featurevectors[1])]]
        assignments = [1, 2, 2, 2, 1]
        cc = cluster_complexities(fvs, assignments)
        @test length(cc) == 2
        @test cc[1] ≈ 0 atol = 1e-12 # nodes 1 and 5 are identical
        @test cc[2] ≈ complexity(fvs[2:4])
    end

    @testset "extensions: clustering and MDS" begin
        using Clustering, MultivariateStats

        ha = HarmonicAnalysis(small_truss(); dims = 8)
        n = length(ha.signatures)

        km = cluster_nodes(ha, 2; maxiter = 500)
        @test km isa Clustering.KmeansResult
        @test length(km.assignments) == n && sort(unique(km.assignments)) == [1, 2]

        cc = cluster_complexities(ha, km.assignments)
        @test length(cc) == 2 && all(>=(0), cc)
        @test maximum(cc) <= complexity(ha) + 1e-10 # clustering can only shrink spheres

        E = embed_nodes(ha; maxoutdim = 2)
        @test size(E) == (2, n)
        # MDS preserves distances up to the embedding truncation
        D = distance_matrix(ha)
        Demb = [norm(E[:, i] - E[:, j]) for i = 1:n, j = 1:n]
        @test maximum(abs.(Demb .- D)) / maximum(D) < 0.5
    end
end
