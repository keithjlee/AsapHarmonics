# AsapOptim extension

Methods added when AsapOptim is loaded (`using AsapOptim`): differentiable
design-to-complexity evaluation over an `OptParams`. The user-facing entry
points ([`harmonic_params`](@ref), [`feature_vectors`](@ref),
design methods of [`feature_matrix`](@ref), [`complexity`](@ref), and
[`soft_complexity`](@ref)) are documented in the
[Analysis layer](analysis.md) page; this page collects the
extension-specific method documentation. See
[Differentiable optimization](../guide/optimization.md) for the walkthrough.

```@autodocs
Modules = [Base.get_extension(AsapHarmonics, :AsapHarmonicsOptimExt)]
```
