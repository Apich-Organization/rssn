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

## Known gaps

This list is deliberately blunt; an entry is removed only when a test pins
the new behaviour.

* **Complex special functions** beyond the five above are missing: complex
  Bessel `Y`/`K`, Airy, elliptic integrals, zeta, polylog for `|z| > 1`
  analytic continuation off the real line is only partly covered
  (`polylog_complex`), error function family is covered by `fresnel`.
  `bessel_*_complex` loses digits (about `1e-9`) for orders with `nu^2 > 2|z|`
  near `|z| = 14`.
* **Dense linear algebra is real only.** Complex SVD / QR / LU / Cholesky /
  eigen-solvers exist only in the generic `kernels::matrix::Matrix<T: Field>`
  (determinant, rank, inverse, null space) and not as LAPACK-quality
  routines. There is no pivoted rank-revealing QR, no Schur decomposition
  exposed, and `lu_factor` only reports an exactly zero pivot (no condition
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
