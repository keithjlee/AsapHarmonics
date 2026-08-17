# Build the AsapHarmonics.jl documentation site.
#
# Local build (the [sources] section in docs/Project.toml dev's the sibling
# repos, so no manual Pkg.develop step is needed):
#   julia --project=docs -e 'using Pkg; Pkg.instantiate()'
#   julia --project=docs docs/make.jl
# then open docs/build/index.html (prettyurls are off outside CI, so direct
# file browsing works).
#
# The package is unregistered and has no CI for now, so there is no
# deploydocs/gh-pages step — the built site is browsed locally. When the
# package goes public, add deploydocs(...) and a Documentation workflow
# mirroring Asap.jl's.

using Documenter, AsapHarmonics

# Load the weak dependencies so the package extensions (and their docstrings
# and guide @examples) are active during the build.
using Asap, AsapOptim, Zygote, Clustering, MultivariateStats

const OptimExt = Base.get_extension(AsapHarmonics, :AsapHarmonicsOptimExt)

makedocs(
    sitename = "AsapHarmonics.jl",
    authors = "Keith J. Lee",
    modules = [AsapHarmonics, OptimExt],
    checkdocs = :exports,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "v2-dev",
        size_threshold = 400 * 2^10,
        size_threshold_warn = 250 * 2^10,
    ),
    pages = [
        "Home" => "index.md",
        "Guide" => [
            "The method" => "guide/method.md",
            "Signatures and feature vectors" => "guide/signatures.md",
            "Distances, complexity, clustering" => "guide/analysis.md",
            "Differentiable optimization" => "guide/optimization.md",
        ],
        "API reference" => [
            "Signatures and drivers" => "api/types.md",
            "Closed-form descriptors" => "api/descriptors.md",
            "Analysis layer" => "api/analysis.md",
            "Sampled path and visualization" => "api/sampled.md",
            "AsapOptim extension" => "api/optim.md",
        ],
    ],
)
