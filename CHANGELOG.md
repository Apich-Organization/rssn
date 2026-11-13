# Changelog

All notable changes to rssn are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). The crate version in
`Cargo.toml` is `3.0.0-rc.1`; entries under "Unreleased" are in the working tree
of the `v0.3` branch and not yet part of a tagged release.

## [Unreleased]

### Added

- **Operators in the `numeric` rule set** (`src/rules/numeric.rs`,
  `src/rules/numeric/*.rs`): robust statistics, Hyndman-Fan quantiles and the
  seeded bootstrap interval, windows and spectra, the fractional Fourier
  transform, one-level, multilevel and 2D wavelets, dense and sparse Krylov
  solves, general algebraic multigrid, nested-dissection ordering, total-degree
  homotopy for polynomial systems, complex special functions, rank-revealing
  and multilinear decompositions (`qr_pivoted_*`, `hosvd`, `hooi`, `cp_als`,
  `tt_svd`, `tt_dense`), tensor and sampled-function norms and products,
  `oscillatory_integral`, `anderson_fixed_point`, `newton_cg_minimize`, the
  named-integrator initial-value problems (`ode_solve`, `ode_symplectic`,
  `ode_imex`, `ode_etdrk4`) and the finite-difference solvers `heat_fd`,
  `wave_fd`, `poisson_fd`, `advect_fd` (dimension from the depth of the data).
  They compute from literal numbers only and cap their input size. The list,
  the argument conventions and the kernels deliberately not exposed are in
  [docs/MATH_COVERAGE.md](docs/MATH_COVERAGE.md). Their signatures are in
  `signatures.overlay.json`.
- **3D and 1D PDE kernels** (`kernels::pde`, new `src/kernels/pde/three.rs`):
  `heat_3d`, `wave_3d`, `poisson_3d`, `laplace_3d`, `advection_3d`, `advection_2d`,
  `poisson_1d`, `laplace_1d`, with `Grid3D`/`Field3D` and the matching
  `PdeProblem` variants. Faces are ordered `[x-, x+, y-, y+, z-, z+]`.
- **Dimension-generic simulation kernels** (`sim::dimensional`) and 1D/2D/3D
  forms of Poisson (Jacobi, Gauss-Seidel, SOR, multigrid), advection, diffusion,
  advection-diffusion, viscous and inviscid Burgers, Crank-Nicolson heat and
  Schrödinger, leapfrog wave, structural stiffness assembly
  (`TetrahedronElement3D`), tensor-product interpolation
  (`bilinear_interpolation`, `trilinear_interpolation`, `interpolate_nd`),
  dot and cross products, and slices and streamlines of 3D fields in `io::plot`.
  The family-by-family table is in
  [docs/MATH_COVERAGE.md](docs/MATH_COVERAGE.md).
- **Dimension-free symbolic PDE and transform names**: `solve_heat_equation`,
  `solve_wave_equation`, `solve_laplace_equation`, `solve_poisson_equation`
  (the `_1d`, `_2d`, `_3d` spellings stay as aliases of the same solver), and
  `shear` (two or six arguments; `shear_3d` is new).
- **GPU `f64` by double-single arithmetic** on the wgpu device (`gpu-wgpu`
  feature, `src/backend/device/ds.rs`): `DeviceInfo::emulated_f64` and
  `dtype_epsilon`. About 48 bits and the exponent range of `f32`; see
  [docs/BACKENDS.md](docs/BACKENDS.md) for the measured error and the limits.
- **Launch-fusion planner** (`backend::device::plan`): `Plan` describes a graph
  of kernels, `fuse` merges producers into consumers, `run` executes the result.
- Integration tests `tests/coverage_algebra.rs`, `tests/coverage_api_complex.rs`,
  `tests/coverage_numeric_a.rs`, `tests/coverage_numeric_b.rs`, and
  `dimension_coverage.rs` under `tests/kernels/` and `tests/sim/`.
- Exact rational fast path in the Lie algebra rules (`rules::lie`): products and
  linear solves with literal rational entries are computed directly.
- This changelog and the documentation index [docs/README.md](docs/README.md).

### Changed

- `build.rs` now requires `signatures.witnessed.json`, exactly like
  `signatures.overlay.json`. Before, a missing file silently became an empty
  table and two builds of the same sources could produce different manifests.
  The build error names the file; `tests/build_script.rs` covers it. See
  [docs/BUILDING.md](docs/BUILDING.md) and [docs/KERNEL.md](docs/KERNEL.md).
- Integer-order Bessel `J_n` and `Y_n` (`n >= 2`, up to order 1000) use the
  upward recurrence from `J_0, J_1` and `Y_0, Y_1` (for `J_n` once `x > n`)
  instead of quadrature.
- `graph::dim` documented-bypass list extended with the `numeric` operators whose
  result dimension no rule expresses.
- `lowrank::Gaussian` mixes its seed with the SplitMix64 finalizer; the random
  streams for a given seed differ from earlier versions.

### Fixed

- `rules::pde::extension_tests::heat_on_annuli_and_shells` returned NaN because
  the annulus eigenvalue scans called the quadrature fallback for integer-order
  Bessel functions, which needs seconds per call at arguments in the tens
  (see the Bessel entry above).
- `is_lie_algebra` and `structure_constants` could hang on bases of literal rational
  matrices; fixed by the exact fast path in `rules::lie`.
- The acoustics frozen-bit tests compared bits recorded on one platform's
  `sin`/`exp`; they now compare to round-off of the field amplitude
  (`1e-11` relative), so they pass with other math libraries. The scheme itself
  stays pinned bit for bit by a separate test.
- The DPD wall-decay test runs 2500 steps instead of 1500 so that its 10% bound
  does not depend on the noise realisation.
- `lowrank::Gaussian` seeds `2k` and `2k + 1` used to give the same generator
  (`seed | 1`); now every seed gives its own stream, and seed 0 is valid.
