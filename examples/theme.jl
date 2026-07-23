# Shared plotting style for the AsapHarmonics examples.
#
# Approximates the visual language of the dissertation figures (kjlMakie's
# light mono theme) without the dependency: pink = compression, blue =
# tension, thin black structural wireframes, undecorated axes.

using CairoMakie

const COMPRESSION_PINK = RGBf(0.961, 0.145, 0.396) # #f52565
const TENSION_BLUE = RGBf(0.145, 0.349, 0.961)     # #2559f5
const FORCE_CMAP = cgrad([COMPRESSION_PINK, RGBf(0.85, 0.85, 0.85), TENSION_BLUE])
const CLUSTER_CMAP = :lightrainbow

set_theme!(Theme(
    font = "JuliaMono", # dissertation figures use monospace labels
    fontsize = 12,
    figure_padding = 20,
    Axis = (
        titlefont = :bold,
        titlesize = 14,
        xgridvisible = false,
        ygridvisible = false,
    ),
))

"undecorated structural view: no spines, ticks, or grid; true aspect"
function structurestyle!(ax::Axis)
    hidedecorations!(ax)
    hidespines!(ax)
    ax.aspect = DataAspect()
    return ax
end

"line segments (pairs of Point2f) of every element of a planar model"
function element_segments(model)
    return [
        Point2f(n.position[1], n.position[2]) for
        e in model.elements for n in (e.nodeStart, e.nodeEnd)
    ]
end

"node coordinates of a planar model"
node_points(model) = [Point2f(n.position[1], n.position[2]) for n in model.nodes]

"line segments (pairs of Point3f) of every element of a spatial model"
function element_segments3(model)
    return [
        Point3f(n.position...) for
        e in model.elements for n in (e.nodeStart, e.nodeEnd)
    ]
end

"node coordinates of a spatial model"
node_points3(model) = [Point3f(n.position...) for n in model.nodes]

"per-element axial forces, and a symmetric colorrange saturated at `frac`"
function force_colors(model; frac = 0.5)
    f = [axial_force(model.results, e) for e in model.elements]
    fmax = maximum(abs, f)
    return f, (-fmax, fmax) .* frac
end

savefig(name, fig) = save(joinpath(@__DIR__, "figures", name * ".png"), fig; px_per_unit = 2)
