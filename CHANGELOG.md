# Changelog

All notable changes to rssn are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). The crate version in
`Cargo.toml` is `3.0.0-rc.1`; entries under "Unreleased" are in the working tree
of the `v0.3` branch and not yet part of a tagged release.

## [Unreleased]

### Added

- Three-dimensional shallow water in `sim::euler_nd` and `sim::euler_amr`: `Law::ShallowWater` now runs on 3D grids as the barotropic (`γ = 2`, `p = g h²/2`) system with state `(h, hu, hv, hw)` (HLLC, Roe and Rusanov, all reconstructions, reflective/periodic/fixed boundaries, AMR); 1D/2D results are bit-identical, and a bed or Manning friction are rejected in 3D. `physics_fvm::solve_shallow_water_3d` is the first-order Lax-Friedrichs sibling of `solve_shallow_water_2d`.

- Numeric operators that close gaps in existing families: `heat_fd_bc` and `wave_fd_bc` (finite differences with `dirichlet` or `neumann` walls in 1D/2D/3D; the kernel already supported Neumann, the bridge only exposed Dirichlet), `advect_fd_upwind` (monotone first-order scheme beside Lax-Wendroff), `ifft_abs` (the missing magnitude beside `ifft_re`/`ifft_im`), `tensor_contract` (generalised trace; the kernel `contract` had no operator), `tucker_dense` and `cp_dense` (inverses of `hosvd`/`hooi` and `cp_als`, as `tt_dense` is of `tt_svd`; kernel `tensor::tucker_reconstruct` and `tensor::cp_reconstruct`), and `erf_complex`, `fresnel_s_complex`, `fresnel_c_complex` (the kernels existed without operators).
- Typed, rank-polymorphic entries in `signatures.overlay.json` for 38 operators that had only the blanket derived or witnessed signatures: `heat_fd`, `wave_fd`, `poisson_fd`, `advect_fd` (rank `n` tied across arguments and result), `minres_solve`, `fgmres_solve`, `bicgstab_l_solve`, `randomized_singular_values`, `lobpcg_eigenvalues`, `tensor_*`, `tensordot`, `function_*`, `wavedec`, `waverec`, `nested_dissection_order`, `polynomial_system_roots`, `hosvd`, `hooi`, `cp_als`, `tt_svd`, `tt_dense`, plus every operator added above. A test now requires every numeric table operator to have an overlay signature of its own arity.

### Changed

- The dimension rules of the new siblings are explicit (`heat_fd_bc`, `wave_fd_bc`, `advect_fd_upwind` keep the field's dimension; `tensor_norm`, `tensor_contract`, `tt_dense`, `tucker_dense`, `cp_dense`, `ifft_abs` keep the tensor's; the complex error and Fresnel functions take bare numbers). `hooi` documents that a rank of `0` means the full extent, as `hosvd` does.

### Fixed

- `tensor::tensordot` (and so `tensordot`, `tensor_vec_mul`) panicked on a contracted axis past the rank or named twice; it now returns an error and the operator stays a request.
