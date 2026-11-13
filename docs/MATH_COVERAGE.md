# Mathematical coverage and edge-case policy

rssn aims at the breadth of NumPy + SciPy + SymPy. This page records what the
numerical kernels and the symbolic rules are checked against, the conventions
they follow on degenerate input, and, as honestly as we can, what is still
missing. The regression tests named below live in `tests/` and run in CI.

## Conventions for degenerate input

* **No panics on any input.** Kernels reachable from the public API, the
  dynamic ABI (`Session::call_value`) or the C interface return an error or
  NaN/infinity instead. `tests/abi_fuzz.rs` calls every operator of the
  manifest that has a plain signature with zeros, NaN, infinities, huge
  magnitudes, empty and 1x1 tensors and inconsistent shapes, and fails on a
  panic or a call that runs for more than 20 seconds.
* **Dense linear algebra** (`kernels::dense`, `kernels::matfun`) returns
  `DenseError::NonFinite` for NaN/infinite entries (instead of silently
  propagating them or, for the eigen-solvers' balancing loop, spinning),
  `Shape` for mismatched lengths, `Singular` only for an exactly zero pivot
  or a non-positive-definite Cholesky input. Rank-deficient problems should
  use `pinv`, `lstsq_svd`, `null_space` or `orth` (all SVD based, default
  tolerance `max(m, n) * eps * sigma_max` as in NumPy).
* **Matrix functions.** `sqrtm` and `logm` use the eigen-decomposition for
  symmetric input (exact for singular positive semi-definite `sqrtm`), and
  Denman-Beavers / inverse scaling and squaring otherwise; a matrix without a
  real principal root is an error.
* **Special functions** return NaN outside the domain, `+-inf` at poles and
  never panic. The complex-argument module `kernels::special_complex`
  provides `gamma`, `ln_gamma` (analytic branch), `digamma`, Lambert `W_k`
  and Bessel `J_nu`, `I_nu` (real order).
* **Statistics.** Empty or NaN-containing samples give NaN; `percentile`
  accepts fractional `p` in `[0, 100]`; the p-values of the t, F and chi-square
  tests use survival functions, so tails below `1e-16` are not rounded to 0;
  `mean` does not overflow for finite data whose sum does.
* **Number theory** (`kernels::number_theory`): everything is exact in
  `u128` intermediates. Primality is deterministic for all of `u64`;
  `factorize` is Pollard-Brent rho and factors any `u64` in milliseconds;
  `mod_inverse` accepts negative and extreme operands; `lcm` saturates and
  `checked_lcm` reports overflow.
* **Quadrature.** Adaptive Simpson has a panel budget (a NaN integrand or a
  zero tolerance can no longer make the recursion exponential), Romberg caps
  its depth, non-finite bounds give NaN, and a vanishing interval is detected
  by exact equality of the endpoints rather than `|a - b| < eps`.

## Where each area is tested

| Area | Test file |
|---|---|
| Dense linear algebra, matrix functions, null space | `tests/math_linear_algebra.rs` |
| Complex gamma, digamma, Lambert W, Bessel | `tests/math_special_complex.rs` |
| Real special functions, edge regions and tails | `tests/math_special_edges.rs` |
| Number theory, primality, factorisation, modular roots, discrete log | `tests/math_number_theory.rs` |
| Statistics, quadrature, interpolation edge cases | `tests/math_numeric_edges.rs` |
| ODE / root finder / optimiser robustness | `tests/math_robustness.rs` |
| Panics anywhere behind the dynamic ABI | `tests/abi_fuzz.rs` |
| Symbolic classics (definite integrals, limits, sums, series, solve, number theory) | `tests/math_symbolic_classics.rs` (and the 450-case `tests/symbolic_battery.rs`) |

## Exact and structured algorithms added in 3.0

| Area | Where | Tested in |
|---|---|---|
| Sparse exact linear algebra: Markowitz-style elimination over any domain (rationals, `GF(p)`, symbolic terms), CRT determinants with Hadamard's bound, Dixon `p`-adic lifting with rational reconstruction, Wiedemann minimal polynomial and solve | `kernels::sparse_exact`, rule set `sparse_symbolic` (`sparse(m, n, entries)`, `sp_*`) | unit tests of both modules |
| D-finite functions: annihilators of atoms, closure under `+`, `*`, powers and `D`, coefficient recurrences, exact zero test | `kernels::holonomic`, rule set `holonomic` (`holonomic_*`) | unit tests of both modules |
| Column-pivoted (rank-revealing) QR, randomized SVD, LOBPCG, Anderson acceleration | `kernels::lowrank` | unit tests |

## Known gaps

This list is deliberately blunt; an entry is removed only when a test pins
the new behaviour.

* **Complex special functions** beyond the six above are missing: complex
  Bessel `Y`/`K`, Airy, elliptic integrals, zeta, polylog for `|z| > 1`
  analytic continuation off the real line is only partly covered
  (`polylog_complex`), error function family is covered by `fresnel`.
  `bessel_*_complex` loses digits (about `1e-9`) for orders with `nu^2 > 2|z|`
  near `|z| = 14`.
* **Dense linear algebra is real only.** Complex SVD / QR / LU / Cholesky /
  eigen-solvers exist only in the generic `kernels::matrix::Matrix<T: Field>`
  (determinant, rank, inverse, null space) and not as LAPACK-quality
  routines. Rank-revealing QR is `kernels::lowrank::qr_pivoted` (real only);
  no Schur decomposition is exposed, and `lu_factor` only reports an exactly zero pivot (no condition
  estimate). `Lu::solve` panics on a right-hand side of the wrong length;
  use `dense::solve`, which checks it.
* **Number theory kernels are 64-bit.** Big-integer primality, ECM, quadratic
  sieve and discrete logarithms beyond `2^48` live in the symbolic
  `number_theory` rule set (arbitrary precision), not as numeric kernels.
  `discrete_log` is baby-step giant-step (memory `O(sqrt p)`); there is no
  Pohlig-Hellman or index calculus.
* **Finite fields.** `kernels::finite_field` has `GF(p)` (64-bit) and
  `GF(256)`; `GF(p^k)` and factorisation over finite fields and number
  fields are in the symbolic rule sets (`gf_*`, `gfx_*`, `gfp_factor`) only.
* **Statistics** wrap `statrs` for distributions; there are no non-parametric
  tests (Mann-Whitney, Kolmogorov-Smirnov), no Cornish-Fisher or extreme-tail
  quantile refinements, and no weighted statistics.
* **Integration** of oscillatory and singular integrands relies on Filon,
  tanh-sinh and Gauss-Kronrod with transformations; there is no general
  Levin / Clenshaw-Curtis oscillatory driver, and `adaptive_quadrature`
  (adaptive Simpson) is a low-order method kept for compatibility: prefer
  `integrate::gauss_kronrod_any` or `quadrature::tanh_sinh`.
* **ODE/SDE.** Finite-time blow-up is reported as an error by the adaptive
  solvers but the fixed-step solvers return `inf`/NaN states as errors only
  once the state stops being finite; DAE index above 1, stiff SDEs and
  delay equations with state-dependent delays are not covered.
* **Symbolic.** `tests/math_symbolic_classics.rs` pins the known failures
  with `.gap()` entries; they are the honest list of classical results rssn
  cannot yet derive.
* **Sparse exact.** `sp_rank` over the rationals eliminates with rational
  arithmetic (no modular certificate); there is no sparse Smith or Hermite
  form, no block Wiedemann, and `sp_minpoly_mod` is Monte Carlo (it returns the
  least common multiple of four projections, which is the minimal polynomial
  with probability at least `1 - 8n/p`).
* **D-finite.** One variable, expansion at the origin only, atoms with
  rational Taylor coefficients; no algebraic substitution, no singular
  points (Bessel functions at `0`), no definite integration (creative
  telescoping lives in `calculus/zeilberger` for sums). The operator found
  is the least-order one in the closure's vector space, not always the
  minimal operator of the function.
* **`kernels::lowrank`** is real, dense-vector code: `lobpcg` has no soft
  locking of converged vectors and no generalised (`B`-orthogonal) form.

## Space dimensions and naming

One name per concept, whatever the number of space dimensions. The dimension
comes from the problem (the variables of the equation, the nesting depth of the
data, the length of the argument list); the older `_1d`/`_2d`/`_3d` spellings
stay as aliases of the same operator, and a family now has all three.

| Family | Dimension-free name | Aliases | Kernel |
|---|---|---|---|
| heat / wave / Laplace / Poisson (symbolic) | `solve_heat_equation`, `solve_wave_equation`, `solve_laplace_equation`, `solve_poisson_equation` | `_1d`, `_2d`, `_3d` of each | `rules::pde` |
| heat / wave / Poisson / advection (finite differences, unit box) | `heat_fd`, `wave_fd`, `poisson_fd`, `advect_fd` | none: depth 1, 2, 3 of the data | `kernels::pde` (`heat_1d` .. `heat_3d`, `poisson_1d` .. `poisson_3d`, `advection_1d` .. `advection_3d`, `wave_1d` .. `wave_3d`) |
| affine transforms | `translation`, `scaling`, `shear`, `reflection`, `rotation` | `_2d`, `_3d` (`shear_3d` added) | `rules::discrete::graphics` |

`method_of_lines_1d` (stiff ODE solvers) is the one finite-difference solver
that is still one-dimensional. Robin and periodic faces for heat, wave and
Poisson are not covered.

## Kernels exposed as operators

`rules::numeric` is table driven (`rules/numeric/bridge.rs` decodes literal
numbers, rectangular or ragged nested lists and symbols naming a method, calls
a kernel and encodes the answer). Besides `fftconv` and `interpolate` it
exposes:

- robust statistics (`mad`, `trimmed_mean`, `huber_location`, `ecdf`, `kde`, ...),
  Hyndman-Fan quantiles of every type (`quantile_hf`) and the seeded percentile
  bootstrap (`bootstrap_ci`, statistic named by a symbol);
- windows and spectra (`hann_window`, `fft_abs`, `cross_correlation`, ...) and the
  fractional Fourier transform (`frft_re`, `frft_im`);
- wavelets by name: one level (`dwt_approx`, `dwt_detail`, `idwt`), multilevel
  (`wavedec`, `waverec`) and the matrix transform (`dwt2`, `idwt2`; the sub-bands
  are `list(LL, LH, HL, HH)`);
- dense Krylov solves (`minres_solve`, `fgmres_solve`, `bicgstab_l_solve`), the
  same four Krylov methods on a sparse triplet matrix with the Jacobi
  preconditioner (`sparse_solve`), general algebraic multigrid (`amg_solve`),
  fill-reducing multilevel nested dissection (`nested_dissection_order`) and
  total-degree homotopy continuation for square polynomial systems
  (`polynomial_system_roots`);
- complex special functions, with `list(re, im)` arguments and results
  (`lgamma_complex`, `gamma_complex`, `digamma_complex`, `lambert_w_complex`,
  `bessel_j_complex`, `bessel_i_complex`);
- low-rank and multilinear: `randomized_singular_values`, `lobpcg_eigenvalues`,
  column-pivoted QR (`qr_pivoted_rank`, `qr_pivoted_permutation`,
  `qr_pivoted_diagonal`), `hosvd`, `hooi`, `cp_als`, `tt_svd`, `tt_dense` (the
  decompositions return lists: core, factors, error);
- tensors (`tensor_norm`, `tensor_inner`, `tensor_outer`, `tensordot`) and sampled
  functions (`function_l2_norm`, `function_inner`, ...);
- the finite-difference solvers above.

`rules/numeric/functions.rs` and `rules/numeric/ode.rs` take expressions:
`oscillatory_integral` (Levin), `anderson_fixed_point`, `newton_cg_minimize`
(symbolic gradient and Hessian) and the initial-value problems `ode_solve`
(`dopri5`, `rosenbrock23`, `radau5`, `bdf`, `sdirk3`, `exprb32`),
`ode_symplectic` (Verlet and Yoshida compositions), `ode_imex` and `ode_etdrk4`.

Every one of them computes from literal numbers only; a symbolic argument, a
ragged list where a rectangular one is needed, a rejected parameter or a failed
convergence leaves the call as a request, and every operator caps its input
size so a call stays bounded. Indices (`qr_pivoted_permutation`,
`nested_dissection_order`, the triplets of `sparse_solve`) are 0-based.

Kernels still reachable only from Rust, and why:

- `differential_geometry`, `coordinates`, `indefinite_sum` and
  `functional_analysis::gram_schmidt` have symbolic counterparts that give
  closed forms (`christoffel`, `riemann`, `ricci`, `to_cartesian`,
  `from_cartesian`, `indefinite_sum`, `gram_schmidt`); a numeric operator beside
  them would answer the same question less exactly.
- `ode_contract` only validates problems for the integrators above.
- `ode_ext` variational Lobatto integrators take a mass matrix and a gradient
  closure, and `etdrk4_diag` a spectral diagonal; both are reachable through the
  expression operators only by restating the problem (`ode_symplectic`,
  `ode_etdrk4`), and are not exposed on their own.
- `ode_sdirk` analytic-Jacobian, fixed-step and event variants (`sdirk3`,
  `sdirk3_events`, `sdirk3_fixed*`): `ode_solve` uses the finite-difference
  Jacobian; events return a time and a state pair that has no single-term form.
- `krylov_ext` with arbitrary preconditioner closures and matrix-free operators,
  and `AmgGeneral` as a reusable preconditioner object: a closure or a built
  hierarchy is not a term. `sparse_solve` and `amg_solve` cover the triplet case.
- `homotopy_cert` (Smale alpha, Krawczyk and interval certificates) returns
  certificates with interval balls; there is no term encoding for them, and
  `polynomial_system_roots` returns the uncertified endpoints.
- `ordering_ml` options (`NdOptions`, `vertex_separator_bisection`) beyond the
  default ordering, and `hybrid_nested_dissection_order`, which gives the same
  kind of answer under another heuristic.
- `signal`: `fft`/`ifft` return complex vectors, so they are exposed as the
  `_re`/`_im`/`_abs` projections; `convolve` is `fftconv`.
- `wavelets`: the `*_with` variants with a boundary rule other than periodic and
  the multilevel 2D transforms (`wavedec2`, `waverec2`) return packed
  coefficient arrays whose layout depends on the boundary rule; use `dwt2` level
  by level on the `LL` band.
- `tensor_train` algebra (`add`, `hadamard`, `round`, `inner`) works on a train
  object; `tt_svd`/`tt_dense` are the conversions.
- `special_complex` has no other public function than the six exposed.

## Numerical notes

Integer-order Bessel functions `J_n`, `Y_n` (`n >= 2`) use the upward
recurrence from `J_0, J_1` and `Y_0, Y_1` (for `J_n` once `x > n`); the
quadrature fallback needs seconds for arguments in the tens, which made the
annulus eigenvalue scans (`annulus_root`) time out into NaN.

## Space-dimension coverage of simulation solvers

Dimension-generic kernels live in `sim::dimensional` (flat `x`-fastest arrays,
`dims = [nx]`, `[nx, ny]` or `[nx, ny, nz]`); the legacy `_1d`/`_2d` names are kept
unchanged and the new `_2d`/`_3d` names are thin wrappers around the generic entry
point. Every added function is tested against an analytic solution, a convergence
order, a conservation law, or the lower-dimensional routine (a 3D problem with data
independent of two axes reproduces the 1D result; axis-swap symmetry) in
`tests/sim/dimension_coverage.rs` and `tests/kernels/dimension_coverage.rs`.

| Family | 1D | 2D | 3D | Generic entry point |
|---|---|---|---|---|
| Poisson, Jacobi / Gauss-Seidel / SOR (`physics_cfd`) | `solve_poisson_1d_{jacobi,gauss_seidel,sor}` (new) | `solve_poisson_2d_{jacobi,gauss_seidel,sor}` | `solve_poisson_3d_{jacobi,gauss_seidel,sor}` (new) | `dimensional::solve_poisson_iterative_nd` |
| Poisson SOR (`physics_fdm`) | `solve_poisson_1d` (new) | `solve_poisson_2d` | `solve_poisson_3d` (new) | `physics_fdm::solve_poisson_nd` |
| Poisson, multigrid (`physics_mtm`) | `solve_poisson_1d_multigrid` | `solve_poisson_2d_multigrid` | `solve_poisson_3d_multigrid` (new) | `solve_poisson_multigrid_nd` (any boundary condition) |
| Advection (`physics_cfd`) | `solve_advection_1d` | `solve_advection_2d` (new) | `solve_advection_3d` (new) | `dimensional::solve_transport_nd` with `Transport::Advection` |
| Diffusion (`physics_cfd`) | `solve_diffusion_1d` | `solve_diffusion_2d` (new) | `solve_diffusion_3d` (new) | `solve_transport_nd`, `Transport::Diffusion` |
| Advection-diffusion (`physics_cfd`, `physics_fdm`, `physics_sm`) | `solve_advection_diffusion_1d` (all three) | `solve_advection_diffusion_2d` (cfd, fdm new; sm existing) | `solve_advection_diffusion_3d` (cfd, fdm new; sm existing) | `solve_transport_nd`, `Transport::AdvectionDiffusion`; `advance_transport_nd` for the final state only |
| Viscous Burgers (`physics_cfd`, `physics_fdm`) | `solve_burgers_1d` | `solve_burgers_2d` (new) | `solve_burgers_3d` (new) | `Transport::Burgers` / `Transport::BurgersCentral` |
| Inviscid Burgers, Lax-Friedrichs (`physics_fvm`) | `solve_burgers_1d` | `solve_burgers_2d` (new) | `solve_burgers_3d` (new) | `physics_fvm::solve_burgers_nd` |
| Shallow water, Lax-Friedrichs (`physics_fvm`) | `solve_shallow_water_1d` | `solve_shallow_water_2d` (new) | not meaningful (depth-averaged model is planar) | high-order 1D-3D Euler / shallow water: `solve_conservation_law_nd` (`euler_nd`) |
| Heat, Crank-Nicolson (`physics_cnm`, `classical`) | `solve_heat_equation_1d_cn`, `solve_heat_equation_1d_crank_nicolson` | `solve_heat_equation_2d_cn_adi` (Peaceman-Rachford); `solve_heat_equation_2d_crank_nicolson` (new) | `solve_heat_equation_3d_cn_adi` (new); `solve_heat_equation_3d_crank_nicolson` (new) | `dimensional::solve_heat_adi_nd`; `classical::solve_heat_equation_nd_crank_nicolson`; general boundary conditions: `physics_cnm::solve_heat_equation_nd_cn` |
| Schrodinger, Crank-Nicolson (`physics_cnm`) | `solve_schrodinger_1d_cn` | `solve_schrodinger_2d_cn` (new) | `solve_schrodinger_3d_cn` (new) | `dimensional::solve_schrodinger_cn_nd` |
| Schrodinger, split-step (`classical`) | `solve_1d_schrodinger` | `solve_2d_schrodinger` | `solve_3d_schrodinger` | - |
| Wave, leapfrog (`classical`) | `solve_wave_equation_1d` | `solve_wave_equation_2d` (new) | `solve_wave_equation_3d` (new) | `dimensional::solve_wave_leapfrog_nd` |
| Structural assembly (`physics_fea`) | `assemble_global_stiffness_matrix`, `assemble_1d_stiffness_matrix` (alias, new) | `assemble_2d_stiffness_matrix` | `assemble_3d_stiffness_matrix` + `TetrahedronElement3D` (new) | `assemble_stiffness_matrix` (any DOF map) |
| Tensor-product interpolation (`kernels::interpolate`) | existing 1D interpolants | `bilinear_interpolation` (new) | `trilinear_interpolation` (new) | `multilinear_interpolation`, `tensor_cubic_spline_interpolation`, `interpolate_nd` |
| Dot / cross product (`kernels::computer_graphics`) | - | `dot_product_2d`, `cross_product_2d` (new) | `dot_product`, `dot_product_3d` (new alias), `cross_product` | `dot` (trait `InnerProduct`), `dot_product_nd` |
| Planar slices and streamlines of 3D fields (`io::plot`) | - | `plot_heatmap_2d`, `plot_contour_2d`, `plot_quiver_2d`, `plot_streamlines_2d` | `plot_heatmap_slice_3d`, `plot_contour_slice_3d`, `plot_quiver_slice_3d`, `plot_streamlines_3d` (new) | `slice_scalar_field_3d` + `SliceAxis` |

Notes on scope. The explicit transport schemes are first order (upwind) or
conditionally stable (FTCS); the multi-axis Burgers equation is the scalar
`u_t + u (u_x + u_y + u_z) = nu Delta u` (the same `u` advects along every axis).
Advection is periodic with period `n - 2` per axis (first and last nodes are
identified with `n - 2` and `1`), as in the legacy 1D routine. Families that
already existed in several dimensions (BEM 2D/3D, DG 3D, `physics_fem` Poisson
scenarios 1D-3D, `physics_sm` advection-diffusion 1D-3D, `structured`
heat/wave/Helmholtz 1D-3D, `euler_nd`) were audited and left unchanged. `quadrature` already has dimension-generic `cubature_*` and
`monte_carlo*` routines; the 2D tensor rule there is a private test helper.
