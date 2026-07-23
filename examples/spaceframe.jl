# # Spaceframe — 3D spherical-harmonic shape descriptors
#
# The 3D pipeline of Lee, Danhaive & Mueller (2022) on a doubly-curved
# spaceframe roof (cf. the paper's example structure): spherical Gaussian
# force signatures per node, closed-form spherical-harmonic feature vectors,
# demand-space clustering, and complexity scoring.
#
# Run from this folder with `julia --project=. spaceframe.jl`.

using Asap, AsapHarmonics
using LinearAlgebra, Statistics, Random
using CairoMakie, Clustering, MultivariateStats

include("theme.jl")

# ## 1. A doubly-curved spaceframe
#
# Asap's `SpaceFrame` generator: 10 × 10 bays of 1 m, ~1 m structural depth,
# top chord lofted by a sine-bump height field, corner-supported, uniform
# downward load on the top nodes.

Random.seed!(1)

section = Section(Material(200e6, 1.0, 80.0, 0.3), 1e-3)
sf = SpaceFrame(10, 1.0, 10, 1.0, 1.0, (u, v) -> 1.5 * sinpi(u) * sinpi(v), section;
    support = :corner, load = [0.0, 0.0, -20.0])
model = sf.model
solve!(model)

# ## 2. Harmonic analysis (closed-form descriptors)

δ = 20 # bump sharpness, as in the paper
ha = HarmonicAnalysis(model; delta = δ, dims = 16)

println("Spaceframe: $(length(ha.signatures)) nodes, $(length(model.elements)) elements")
println("  design complexity (bounding-sphere radius): ", round(complexity(ha), digits = 2))
println("  soft complexity (RMS from centroid):        ", round(soft_complexity(ha), digits = 2))

# ## 3. Structure and internal forces

pts = node_points3(model)
segments = element_segments3(model)
forces, crange = force_colors(model; frac = 0.5)

function structure_axis(pos; title = "")
    ax = Axis3(pos; title = title, aspect = :data, protrusions = 0, elevation = 0.25, azimuth = -1.1π)
    hidedecorations!(ax)
    hidespines!(ax)
    return ax
end

begin
    fig = Figure(size = (900, 450))
    ax = structure_axis(fig[1, 1]; title = "AXIAL FORCES")
    linesegments!(
        ax, segments;
        color = repeat(forces, inner = 2), colorrange = crange, colormap = FORCE_CMAP,
        linewidth = 1.5,
    )
    savefig("spaceframe_forces", fig)
end

# ## 4. Spherical force signatures
#
# Signatures drawn as radially-offset unit spheres (sampled on demand purely
# for drawing, on the same grid convention as the numeric validation path):
# protrusions = tension demand, indentations = compression. A sample of nodes
# is shown — a support corner, top-chord apex nodes, and interior bottom
# nodes — in the manner of the paper's Fig. 5.

"radially-scaled signature surface of a node signature"
function signature_surface(sig; delta = δ, nlat = 60)
    F = sampled_force_function(sig; delta = delta, nlat = nlat)
    Θ, Φ = sphere_points(nlat)
    Φw = vcat(Φ, 2π)          # close the seam for plotting
    Fw = hcat(F, F[:, 1])
    fmax = maximum(abs, Fw)
    R = fmax > 0 ? 1.0 .+ Fw ./ fmax : ones(size(Fw))
    x = make_xsphere(Θ, Φw) .* R
    y = make_ysphere(Θ, Φw) .* R
    z = make_zsphere(Θ, Φw) .* R
    return x, y, z, Fw
end

begin
    isample = unique([sf.isupport[1]; vec(sf.itop)[[1, 45, 55, 100]]; vec(sf.ibottom)[[13, 61]]])
    n = length(isample)
    fig = Figure(size = (220 * n, 240))

    for (j, i) in enumerate(isample)
        local sig = ha.signatures[i]
        local x, y, z, F = signature_surface(sig)
        local a = Axis3(fig[1, j]; aspect = :data, protrusions = 0,
            title = "node $(i) (:$(sig.id))", titlesize = 11)
        hidedecorations!(a)
        hidespines!(a)
        surface!(
            a, x, y, z;
            color = F, colormap = FORCE_CMAP,
            colorrange = (-maximum(abs, F), maximum(abs, F)),
        )
    end

    savefig("spaceframe_signatures", fig)
end

# ## 5. Distance matrix, clustering, demand space

D = distance_matrix(ha)
k = 10
km = cluster_nodes(ha, k; maxiter = 5000)
E = embed_nodes(ha; maxoutdim = 2)

begin
    fig = Figure(size = (1200, 420))

    ax0 = Axis(
        fig[1, 1];
        title = "NODAL DISTANCE MATRIX", xlabel = "node", ylabel = "node",
        yreversed = true, aspect = DataAspect(),
    )
    heatmap!(ax0, D, colormap = cgrad([:white, TENSION_BLUE]))

    ax = structure_axis(fig[1, 2]; title = "NODE CLUSTERS")
    linesegments!(ax, segments, color = (:black, 0.25), linewidth = 1)
    scatter!(
        ax, pts;
        color = km.assignments, colormap = CLUSTER_CMAP, colorrange = (1, k),
        strokecolor = :black, strokewidth = 1, markersize = 9,
    )

    ax2 = Axis(
        fig[1, 3];
        title = "2D PROJECTION", xlabel = "dim 1", ylabel = "dim 2",
        backgroundcolor = RGBf(0.96, 0.96, 0.96), aspect = DataAspect(),
    )
    scatter!(
        ax2, Point2f.(eachcol(E));
        color = km.assignments, colormap = CLUSTER_CMAP, colorrange = (1, k),
        strokecolor = :black, strokewidth = 1, markersize = 9,
    )

    savefig("spaceframe_clustering", fig)
end

cc = cluster_complexities(ha, km.assignments)
println("  per-cluster complexity (k = $k): ", round.(cc, digits = 2))
println("figures written to examples/figures/")
