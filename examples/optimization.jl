# # Minimizing nodal complexity by gradient descent
#
# The capability the closed-form descriptors unlock: the design-complexity
# surrogate is smooth in the design variables, so a spaceframe's shape can be
# optimized to REDUCE the variation in its nodal force demands — fewer unique
# connection designs for the same topology — with gradients through the full
# chain (design vector → linear solve → axial forces → signatures → feature
# vectors → complexity) via AsapOptim + Zygote.
#
# Run from this folder with `julia --project=. optimization.jl`.

using Asap, AsapOptim, AsapHarmonics
using LinearAlgebra, Statistics, Random
using Zygote
using CairoMakie

include("theme.jl")

# ## 1. Base design: a shallow vault spaceframe

Random.seed!(1)

section = Section(Material(200e6, 1.0, 80.0, 0.3), 1e-3)
sf = SpaceFrame(8, 1.0, 8, 1.0, 0.8, (u, v) -> 1.0 * sinpi(u) * sinpi(v), section;
    support = :corner, load = [0.0, 0.0, -20.0])
model = sf.model
solve!(model)

# NOTE: `updatemodel(p, x)` (step 4) materializes the optimized design INTO
# the reference model held by OptParams — it mutates `model` in place and
# returns it. Keep an untouched copy for the before/after comparison.
model_base = deepcopy(model)

# ## 2. Design variables and the differentiable objective
#
# Every top-chord node may move vertically. `harmonic_params` precompiles the
# constant per-node gather data; `soft_complexity(x, p, hp)` is then smooth in
# the design vector.

vars = [SpatialVariable(model.nodes[i], 0.0, -0.8, 1.5, :Z) for i in vec(sf.itop)]
p = OptParams(model, vars)
hp = harmonic_params(p; delta = 20, dims = 16)

x0 = copy(p.values)
obj(x) = soft_complexity(x, p, hp)

println("initial soft complexity: ", round(obj(x0), digits = 3))
println("initial exact complexity: ", round(complexity(x0, p, hp), digits = 3))

# ## 3. Projected gradient descent with a fixed step

x = copy(x0)
history = [obj(x)]
steps = 300
η = 0.02

for _ = 1:steps
    g = Zygote.gradient(obj, x)[1]
    global x = clamp.(x - η * g / (norm(g) + 1e-12), p.lb, p.ub)
    push!(history, obj(x))
end

println("final soft complexity:   ", round(history[end], digits = 3),
    "  (", round(100 * (1 - history[end] / history[1]), digits = 1), "% reduction)")
println("final exact complexity:  ", round(complexity(x, p, hp), digits = 3))

# ## 4. Before / after
#
# `updatemodel` materializes the optimized design back into the (solved)
# reference model — mutating `model`; `model_base` holds the original.

model_opt = updatemodel(p, x)

function force_view(pos, m; title = "")
    ax = Axis3(pos; title = title, aspect = :data, protrusions = 0,
        elevation = 0.2, azimuth = -1.1π)
    hidedecorations!(ax)
    hidespines!(ax)
    f, crange = force_colors(m; frac = 0.5)
    linesegments!(
        ax, element_segments3(m);
        color = repeat(f, inner = 2), colorrange = crange, colormap = FORCE_CMAP,
        linewidth = 1.5,
    )
    return ax
end

begin
    fig = Figure(size = (1100, 420))

    force_view(fig[1, 1], model_base; title = "BASE DESIGN")
    force_view(fig[1, 2], model_opt; title = "COMPLEXITY-OPTIMIZED")

    ax = Axis(
        fig[1, 3];
        title = "OBJECTIVE", xlabel = "iteration", ylabel = "soft complexity",
    )
    lines!(ax, 0:steps, history, color = TENSION_BLUE, linewidth = 2)

    savefig("optimization", fig)
end

# ## 5. Demand-space contraction
#
# The point of the exercise: after optimization the nodal feature vectors
# huddle closer together — nodes have become more similar, and fewer unique
# connections cover the design.
#
# MDS coordinates are only meaningful up to rotation/reflection/translation,
# so two separate embeddings cannot be compared point to point. Instead, BOTH
# feature-vector sets are pooled into one distance matrix and embedded
# JOINTLY: all base-vs-base, optimized-vs-optimized, AND base-vs-optimized
# distances shape the shared space, so the two clouds sit where they actually
# live relative to each other.

ha0 = HarmonicAnalysis(model_base; delta = 20, dims = 16)
ha1 = HarmonicAnalysis(model_opt; delta = 20, dims = 16)

begin
    using MultivariateStats

    nn = length(ha0.featurevectors)
    E = embed_nodes(vcat(ha0.featurevectors, ha1.featurevectors); maxoutdim = 2)

    fig = Figure(size = (760, 480))
    ax = Axis(
        fig[1, 1];
        title = "DEMAND SPACE (joint embedding)",
        xlabel = "dim 1", ylabel = "dim 2",
        backgroundcolor = RGBf(0.96, 0.96, 0.96), aspect = DataAspect(),
    )

    scatter!(
        ax, Point2f.(eachcol(E[:, 1:nn]));
        color = (:gray60, 0.9), strokecolor = :black, strokewidth = 0.5, markersize = 8,
        label = "base — complexity $(round(complexity(ha0), digits = 2))",
    )
    scatter!(
        ax, Point2f.(eachcol(E[:, nn+1:end]));
        color = TENSION_BLUE, strokecolor = :black, strokewidth = 0.5, markersize = 8,
        label = "optimized — complexity $(round(complexity(ha1), digits = 2))",
    )
    axislegend(ax; position = :rb, framevisible = false)

    savefig("optimization_demandspace", fig)
end

println("figures written to examples/figures/")
