# AsapHarmonics.jl

Connection analysis over Asap.jl (`../Asap`): builds spherical/circular Gaussian force "signatures" per node and FFT / spherical-harmonic feature vectors (`NodeForces`, `HarmonicAnalysis`). Deps: FFTW, FastSphericalHarmonics.

**Asap is undergoing a v1.0 modernization (`../Asap/docs/MODERNIZATION.md`); this package migrates in lockstep at Phase 5c** — a light rename pass. Coupling is narrow: `connectivity`, `axial_force`, `TrussNode`/`TrussModel`, and direct reads of `node.reaction`, `node.nodeID`, `load.value`, `element.LCS[1]`, element end positions. In the new core, results move off structs (`node.reaction` → `reaction(results, node)`), and `TrussNode`/`TrussModel` are deleted in favor of unified `Node`/`Model`.

A newer rewrite of the same ideas lives in `../AsapHarmonics2` — check both when Asap's reaction/LCS access patterns change.

## Commands

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
