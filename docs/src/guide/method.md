# The method

This page walks through the full process — from a solved structure to a
design-complexity score — and derives the closed-form descriptors that make
the pipeline exact, fast, and differentiable. It condenses the derivation in
the repository README (which carries the complete set of primary citations);
the numbered references at the bottom match the README's.

## 1. From a solved model to nodal force signatures

Start from a solved Asap `Model`. At each node, collect every force demand
acting there:

- **Members**: each incident element contributes its signed axial force
  ``f_i`` (tension positive) placed at the *outward* unit direction
  ``\hat{n}_i`` — the direction pointing from the node along the member.
- **Applied loads** and **support reactions**: each contributes its magnitude
  ``|P|`` placed at its direction ``\hat{P}``.

This demand set ``\{(f_i, \hat{n}_i)\}`` is what a connection at that node
must resolve. To compare demand sets *as shapes* — independent of how many
forces there are or how they are ordered — each set is converted to a smooth
function on the unit sphere, the **force signature**:

```math
f(\hat{x}) = \sum_i f_i\, e^{-\delta \lVert \hat{x} - \hat{n}_i \rVert^2},
\qquad \hat{x} \in S^2,
```

a sum of Gaussian bumps measured by squared *chord* distance, with sharpness
``\delta`` (each bump approaches a Dirac spike as ``\delta \to \infty``; the
paper uses ``\delta = 20``). This is the deterministic force-to-shape
conversion of Lee, Danhaive & Mueller [1], adapting the spherical-harmonic
shape descriptors of Kazhdan, Funkhouser & Rusinkiewicz [2] to structural
force demands. For planar trusses the analog lives on the circle, with bumps
in geodesic (arc) angle following the planar Fourier descriptors of Zahn &
Roskies [3]:

```math
s(\theta) = \sum_i f_i\, e^{-(\theta - \theta_i)^2 / 2\sigma^2}.
```

In the package, signature extraction is [`NodeSignature`](@ref) and the
whole-model drivers are [`HarmonicAnalysis`](@ref) (3D) and
[`HarmonicAnalysis2d`](@ref) (2D) — see
[Signatures and feature vectors](signatures.md).

## 2. From signature to descriptor — why there is a closed form

The classical route (and the original implementation of [1]) samples
``f(\hat{x})`` on a latitude–longitude grid, runs a discrete
spherical-harmonic transform, and reads off band energies. AsapHarmonics
skips all of that, because the chord-distance Gaussian is secretly a **zonal
kernel**. For unit vectors,
``\lVert \hat{x} - \hat{n} \rVert^2 = 2 - 2\, \hat{x}\cdot\hat{n}``, so each
bump is exactly

```math
e^{-\delta \lVert \hat{x} - \hat{n} \rVert^2}
  = e^{-\kappa}\, e^{\kappa\, \hat{x} \cdot \hat{n}},
\qquad \kappa = 2\delta
```

— the kernel of the **von Mises–Fisher distribution** on the sphere (Fisher
[4]; Mardia & Jupp [5]). It depends on ``\hat{x}`` only through
``\hat{x}\cdot\hat{n}``, and functions of that form are precisely the class
with closed-form spherical-harmonic expansions.

Two classical results finish the job:

**Funk–Hecke / Gegenbauer** (Funk [8], Hecke [9]; DLMF §10.60.7 [6]): a
single bump of weight ``f_i`` at ``\hat{n}_i`` has real orthonormal
spherical-harmonic coefficients

```math
a_{lm} = f_i\, \lambda_l\, Y_{lm}(\hat{n}_i),
\qquad
\lambda_l = 4\pi\, e^{-\kappa}\, i_l(\kappa),
```

with ``i_l`` the modified spherical Bessel functions — evaluated stably as
``e^{-\kappa} i_l(\kappa) = \sqrt{\pi/2\kappa}\,\mathrm{besselix}(l+\tfrac12,
\kappa)``. In the package: [`zonal_coefficients`](@ref). Because ``\delta``
is a constant hyperparameter, these coefficients are computed once — no
differentiation ever passes through a Bessel function.

**Addition theorem** (Atkinson & Han [11]): summing
``Y_{lm}(\hat{n}_i)\,Y_{lm}(\hat{n}_j)`` over ``m`` collapses to a Legendre
polynomial of the angle between the two directions,

```math
\sum_{m=-l}^{l} Y_{lm}(\hat{n}_i)\, Y_{lm}(\hat{n}_j)
  = \frac{2l+1}{4\pi}\, P_l(\hat{n}_i \cdot \hat{n}_j).
```

The **feature vector** is the ``L^2`` norm of each degree band of the full
signature, and combining the two results turns it into a pairwise sum over
the node's force directions — no grid, no quadrature, no transform:

```math
\mathrm{FV}_l^2 = \sum_{m=-l}^{l} a_{lm}^2
  = \lambda_l^2\, \frac{2l+1}{4\pi}
    \sum_{i,j} f_i f_j\, P_l(\hat{n}_i \cdot \hat{n}_j).
```

In the package: [`spherical_feature_vector`](@ref), with
[`pairwise_legendre_sums`](@ref) evaluating every degree in a single upward
recurrence pass over the force pairs. Cost is
``O(\text{valence}^2 \cdot \text{dims})`` per node. On the circle the same
logic goes through the **wrapped Gaussian**, whose Fourier coefficients are
``\tfrac{\sigma}{\sqrt{2\pi}} e^{-\sigma^2 k^2/2}`` (Mardia & Jupp [5]):

```math
\mathrm{FV}_k
  = \frac{\sigma}{\sqrt{2\pi}}\, e^{-\sigma^2 k^2/2}
    \sqrt{\textstyle\sum_{i,j} f_i f_j \cos k(\theta_i - \theta_j)}
```

— [`circular_feature_vector`](@ref).

!!! note "The sampled path still exists"
    The grid/transform pipeline ([`sampled_force_function`](@ref) +
    spherical transform, [`circular_gaussian`](@ref) + FFT) is kept for
    cross-validation tests and for visualizing signatures as surfaces. It
    reproduces the closed forms to machine precision — see
    [Signatures and feature vectors](signatures.md#Visualizing-signatures).

## 3. What the descriptor guarantees

- **Exact rotation invariance.** A rotation of the whole demand set only
  mixes spherical-harmonic coefficients within each degree (or shifts
  Fourier phases); band energies are unchanged. Two nodes whose demand
  geometries differ by any rigid rotation are *identical* to the descriptor.
  This is exact, not approximate: the angle-coordinate kernel as written in
  [1] is only approximately zonal near the poles, whereas the chord-distance
  kernel actually implemented is exactly zonal.
- **Scaling and polarity.** ``\mathrm{FV}(cf) = |c|\,\mathrm{FV}(f)`` and
  ``\mathrm{FV}(-f) = \mathrm{FV}(f)``; the descriptor is also invariant to
  force ordering. Doubling every force doubles the descriptor; a mirror-image
  demand pattern matches its original.
- **Smoothness.** ``\mathrm{FV}`` is smooth in magnitudes *and* directions —
  the property that lets design gradients flow through the whole pipeline
  (see [Differentiable optimization](optimization.md)).
- **``\mathrm{FV}_1`` is the equilibrium residual.** For a *complete*
  signature (members + loads + reactions), the degree-1 band is proportional
  to ``\lVert \sum_i f_i \hat{n}_i \rVert`` — the nodal force residual, which
  vanishes at every equilibrated node. ``\mathrm{FV}_1`` therefore carries no
  design information; it is a built-in equilibrium check.

!!! warning "A known degeneracy"
    Where a reaction (or load) is *exactly* collinear with a lone balancing
    member — e.g. the support of a bare pinned two-bar truss — the opposing
    bumps cancel and the signature is legitimately zero. This is a faithful
    statement of the sign convention (demands that cancel along a line leave
    no net shape), but it means such support nodes look "demand-free". If the
    distinction matters for a study, note it when interpreting support-node
    clusters. A physically-placed support direction (a proposed Asap `Node`
    extension) would resolve this in a rotation-equivariant way.

## 4. From descriptors to engineering measures

With every node reduced to a feature vector, comparisons become geometry in
descriptor space (see [Distances, complexity, clustering](analysis.md)):

- **Dissimilarity** between two nodes is the Euclidean distance between
  their feature vectors ([`distance_matrix`](@ref)).
- **Design complexity** [1] is the radius of the minimal bounding hypersphere
  of all nodal feature vectors ([`complexity`](@ref), exact Welzl algorithm
  [12]) — how *spread out* the structure's connection demands are.
- The bounding radius is a minimax and not smooth, so
  [`soft_complexity`](@ref) — the RMS distance of feature vectors from their
  centroid, bounded by twice the exact radius — serves as the
  **differentiable objective** for optimization.
- **Connection standardization**: k-means clustering of feature vectors
  ([`cluster_nodes`](@ref)) groups nodes that could share a connection
  design, with per-cluster residual radii ([`cluster_complexities`](@ref))
  quantifying what standardization leaves unresolved.

## References

Numbering matches the repository README, which lists the complete citations:
[1] Lee, Danhaive & Mueller (2022); [2] Kazhdan, Funkhouser & Rusinkiewicz
(2003); [3] Zahn & Roskies (1972); [4] Fisher (1953); [5] Mardia & Jupp
(2000); [6] DLMF §10.60; [8] Funk (1916); [9] Hecke (1918); [10] Müller
(1966); [11] Atkinson & Han (2012); [12] Welzl (1991).
