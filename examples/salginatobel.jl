# # Salginatobel — 2D Fourier shape descriptors of nodal force demands
#
# The complete 2D pipeline on a trussed interpretation of Maillart's
# Salginatobel bridge, following chapter 4 of the author's dissertation
# (`../resources/kjl_dissertation_chapter4.pdf`) and reproducing its figure
# suite with the AsapHarmonics v2 API:
#
#   1. build and solve the truss
#   2. per-node circular force signatures (sums of angular Gaussian bumps)
#   3. rotation-invariant Fourier feature vectors (closed form)
#   4. distance matrix and demand-space (MDS) embedding
#   5. k-means clustering for connection standardization
#   6. complexity scores, per cluster and vs. number of clusters
#
# Run from this folder with `julia --project=. salginatobel.jl`.

using Asap, AsapHarmonics
using LinearAlgebra, Statistics, Random, JSON
using CairoMakie, Clustering, MultivariateStats

include("theme.jl")

# ## 1. Model from data
#
# The JSON stores node coordinates, element connectivity, supports, and
# loaded nodes (0-indexed). The truss is planar in XY: every node has its
# out-of-plane translation fixed; supports are pinned.

data = JSON.parsefile(joinpath(@__DIR__, "data", "salginatobel.json"))

X = Float64.(data["X"])
Y = Float64.(data["Y"])
ISTART = Int.(data["ISTART"]) .+ 1
IEND = Int.(data["IEND"]) .+ 1
ISUPPORT = Int.(data["ISUPPORT"]) .+ 1
ILOAD = Int.(data["ILOAD"]) .+ 1

free_rotations = [true, true, true] # inactive for truss-only models
planar_free = vcat([true, true, false], free_rotations)
pinned = vcat([false, false, false], free_rotations)

nodes = [
    Node([x, y, 0.0], i in ISUPPORT ? pinned : planar_free, i in ISUPPORT ? :support : :node)
    for (i, (x, y)) in enumerate(zip(X, Y))
]

d = 8 * 0.0254 # 8 in square section
section = Section(Material(15e6, 1.0, 80.0, 0.3), d^2)
elements = AbstractElement{Float64}[
    TrussElement(nodes[i], nodes[j], section) for (i, j) in zip(ISTART, IEND)
]

P = 100.0 # deck load per loaded node [kN]
loads = AbstractLoad{Float64}[NodeForce(nodes[i], [0.0, -P, 0.0]) for i in ILOAD]

model = Model(nodes, elements, loads)
solve!(model)

# ## 2. Harmonic analysis
#
# σ (the `delta` keyword, matching the dissertation) is the angular width of
# each force bump; `dims` Fourier frequencies k = 0:dims-1 form the feature
# vector. Feature vectors are computed in closed form — no sampling, no FFT.

σ = 0.05
ha = HarmonicAnalysis2d(model; delta = σ, dims = 16)

println("Salginatobel: $(length(ha.signatures)) nodes")
println("  design complexity (bounding-sphere radius): ", round(complexity(ha), digits = 2))
println("  soft complexity (RMS from centroid):        ", round(soft_complexity(ha), digits = 2))

# ## 3. Structure and internal forces

pts = node_points(model)
segments = element_segments(model)
forces, crange = force_colors(model; frac = 0.25)

begin
    fig = Figure(size = (900, 500))

    ax = Axis(fig[1, 1], title = "SALGINATOBEL")
    structurestyle!(ax)
    linesegments!(ax, segments, color = :black, linewidth = 1)
    scatter!(ax, pts, color = :white, strokecolor = :black, strokewidth = 1, markersize = 7)

    ax2 = Axis(fig[2, 1], title = "AXIAL FORCES")
    structurestyle!(ax2)
    linesegments!(
        ax2, segments;
        color = repeat(forces, inner = 2), colorrange = crange, colormap = FORCE_CMAP,
        linewidth = 2,
    )

    savefig("salginatobel_forces", fig)
end

# ## 4. Force-signature array
#
# Each node's signature drawn as a closed polar curve: the unit circle,
# radially offset by the normalized signature (protrusions = tension demand,
# indentations = compression demand) — the dissertation's signature array.
# Signatures are sampled on demand purely for drawing.

"closed polar curve of a signature, offset from the unit circle"
function signature_curve(sig; sigma = σ, n = 240)
    ff = circular_gaussian(sig, sigma, n)
    fmax = maximum(abs, ff)
    r = fmax > 0 ? 1.0 .+ ff ./ fmax : fill(1.0, n)
    θ = range(0, 2π, n + 1)[1:end-1]
    curve = Point2f.(r .* cos.(θ), r .* sin.(θ))
    return push!(curve, curve[1]) # close the loop
end

begin
    ncols = 11
    nrows = ceil(Int, length(ha.signatures) / ncols)
    fig = Figure(size = (110 * ncols, 110 * nrows))

    for (i, sig) in enumerate(ha.signatures)
        r, c = fldmod1(i, ncols)
        local a = Axis(fig[r, c])
        structurestyle!(a)
        limits!(a, -2.2, 2.2, -2.2, 2.2)
        lines!(a, signature_curve(sig), color = :black, linewidth = 1)
    end

    savefig("salginatobel_signatures", fig)
end

# ## 5. Feature-vector array
#
# The fixed-length, rotation-invariant descriptor of every node. Note k = 1
# is always ≈ 0: that band measures the resultant of a complete signature,
# which vanishes at equilibrated nodes.

begin
    fvmax = maximum(maximum, ha.featurevectors)
    ncols = 11
    nrows = ceil(Int, length(ha.featurevectors) / ncols)
    fig = Figure(size = (110 * ncols, 90 * nrows))

    for (i, fv) in enumerate(ha.featurevectors)
        r, c = fldmod1(i, ncols)
        local a = Axis(fig[r, c])
        hidedecorations!(a)
        hidespines!(a)
        ylims!(a, 0, 1.05 * fvmax)
        barplot!(a, 0:length(fv)-1, fv, color = :black, gap = 0.15)
    end

    savefig("salginatobel_featurevectors", fig)
end

# ## 6. Distance matrix

D = distance_matrix(ha)

begin
    fig = Figure(size = (560, 520))
    ax = Axis(
        fig[1, 1];
        title = "NODAL DISTANCE MATRIX", xlabel = "node", ylabel = "node",
        yreversed = true, aspect = DataAspect(),
    )
    hm = heatmap!(ax, D, colormap = cgrad([:white, TENSION_BLUE]))
    Colorbar(fig[1, 2], hm, label = "‖FVᵢ − FVⱼ‖")
    savefig("salginatobel_distance_matrix", fig)
end

# ## 7. Clustering and the demand space
#
# k-means over the feature vectors groups nodes with similar force demands —
# candidate groups for standardized connections. Classical MDS embeds the
# distance matrix in 2D for inspection (the dissertation's "demand space").

Random.seed!(1)
k = 12
km = cluster_nodes(ha, k; maxiter = 5000)
E = embed_nodes(ha; maxoutdim = 2)

begin
    fig = Figure(size = (1000, 480))

    ax = Axis(fig[1, 1], title = "NODE CLUSTERS")
    structurestyle!(ax)
    linesegments!(ax, segments, color = :black, linewidth = 1)
    scatter!(
        ax, pts;
        color = km.assignments, colormap = CLUSTER_CMAP, colorrange = (1, k),
        strokecolor = :black, strokewidth = 1, markersize = 11,
    )

    ax2 = Axis(
        fig[1, 2];
        title = "2D PROJECTION", xlabel = "dim 1", ylabel = "dim 2",
        backgroundcolor = RGBf(0.96, 0.96, 0.96), aspect = DataAspect(),
    )
    scatter!(
        ax2, Point2f.(eachcol(E));
        color = km.assignments, colormap = CLUSTER_CMAP, colorrange = (1, k),
        strokecolor = :black, strokewidth = 1, markersize = 11,
    )

    savefig("salginatobel_clustering", fig)
end

# ## 8. Complexity: per cluster, and vs. number of clusters
#
# Per-cluster bounding-sphere radii are the residual standardization penalty
# after grouping; sweeping k shows the diminishing returns of adding more
# unique connection designs.

cc = cluster_complexities(ha, km.assignments)

krange = 1:20
mean_cc = map(krange) do kk
    Random.seed!(1)
    a = kk == 1 ? ones(Int, length(ha.signatures)) : cluster_nodes(ha, kk; maxiter = 5000).assignments
    mean(cluster_complexities(ha, a))
end

begin
    fig = Figure(size = (1000, 400))

    ax = Axis(
        fig[1, 1];
        title = "CLUSTER COMPLEXITY (k = $k)", xlabel = "cluster", ylabel = "radius",
        xticks = 1:k,
    )
    barplot!(ax, 1:k, cc, color = 1:k, colormap = CLUSTER_CMAP, colorrange = (1, k),
        strokecolor = :black, strokewidth = 1)

    ax2 = Axis(
        fig[1, 2];
        title = "MEAN CLUSTER COMPLEXITY vs k", xlabel = "number of clusters",
        ylabel = "mean radius",
    )
    scatterlines!(ax2, krange, mean_cc, color = :black, markercolor = TENSION_BLUE)

    savefig("salginatobel_complexity", fig)
end

println("  per-cluster complexity (k = $k): ", round.(cc, digits = 2))
println("figures written to examples/figures/")
