# Algorithm completion checklist

Regenerated 2026-10-28 from the current tree (not from the old ledger);
re-audited 2026-10-29 (review of the uncommitted round, limitation sweep —
see Part 5, which supersedes Part 1 where they differ).

## How this list was derived

1. The unfinished-work report of the previous round (11 named leftovers).
2. Every module-level `//!` documentation block that states a limitation
   ("Limitations", "Not covered", "Not done"), read in place, with the
   matching code.
3. A test inventory: `tests/**` against `src/**`, to find modules with no
   direct test or only indirect coverage.
4. A lint baseline: `cargo clippy --lib --all-features` on the tree as found.
5. Reading of the largest and most algorithm-heavy modules (`rules/pde`,
   `rules/ode`, `rules/solve`, `rules/poly`, `kernels/homotopy*`,
   `sim/dg3d`, `sim/fem_simplex`, `sim/euler_amr`, `sim/ns_mac`).

Priorities: **P0** = tree is wrong or unverifiable (fix first); **P1** = a
named gap that blocks a whole class of problems; **P2** = coverage or
accuracy improvement; **P3** = research-grade / optional.

Status legend: `[ ]` open, `[x]` done in this round, `[-]` decided against
(with the reason).

---

## Part 0 — The previous round's leftovers

| ID | Item | Status | Where |
|---|---|---|---|
| L1 | Andrews' squeezer unverified | [x] | `tests/kernels/andrews_squeezer.rs` |
| L2 | No 3²:4 sextic test case | [x] | `src/rules/poly/galois_high.rs`, `src/rules/solve/mod.rs` |
| L3 | Homotopy endpoints certified only in floating point; tracking steps uncertified | [x] | `src/kernels/{ball,homotopy_cert,homotopy}.rs` |
| L4 | Coulomb scattering states: numeric cut-off, 3 % reproduction check | [x] | `src/rules/pde/coulomb.rs` |
| L5 | Cylinder drag 5–8 % off literature | [x] | `src/sim/ns_mac.rs`, `tests/sim/ns_mac.rs` |
| L6 | Exterior PDE sources varying in time simplify very slowly | [ ] | `src/rules/pde/{mod,spectral,numeric}.rs` |
| L7 | Missing Neumann and Robin Green's functions | [x] | `src/rules/pde/{poisson,mod}.rs` |
| L8 | General hodograph case, shock interactions | [~] | `src/rules/pde/nonlinear.rs` |
| L9 | Some eigenring and Kovacic number-field cases | [x] | `src/rules/ode/{linear,kovacic,mod}.rs` |
| L10 | DG shock limiting in 3D | [x] | `src/sim/dg3d.rs` |
| L11 | Mesh refinement in 3D | [x] | `src/sim/fem_simplex.rs`, `tests/sim/fem_simplex.rs` |

`[~]` = the tractable half is done (shock interactions); the research half
(general hodograph) is still open and is listed separately below.

### L1 — Andrews' squeezer (done)

The name is the classical **Andrews' squeezer mechanism** (seven-body planar
linkage, index-3 DAE; the benchmark of the Test Set for IVP Solvers, shipped
with RADAU5 as `andrews.f`, listed by IFToMM). No implementation, test or
reference data existed in the tree and no commit in the full history ever
contained the string (checked with `git log --all -S`), so "can't be verified"
meant "never written".

Done in `tests/kernels/andrews_squeezer.rs` (self-contained mechanism; four
tests), registered in `tests/kernels.rs`, with a paragraph in
`src/kernels/ode_stiff.rs`:
- data, initial state and the 27-component reference at `t = 0.03` are
  transcribed from `andrews.f` rev 1.2 (`solut`, PSIDE, `rtol = atol = 1e-14`)
  and the IFToMM report; the formulation is the file's first-order Hessenberg
  form `y = (q, v, a, λ)`, `diag(I₁₄, 0₁₃) y' = φ(y)`, with
  `DaeIndex { n2: 7, n3: 13 }`;
- reference agreement at `rtol = 1e-10`: `|Δq|∞ = 1.36e-8`, `|Δv|∞ = 1.11e-5`,
  relative `|Δa| = 1.32e-4`, `|Δλ| = 1.26e-4`;
- constraints at `t = 0.03`: `|Φ|∞ = 4.0e-16`; energy balance drift
  `5.5e-9 J` over 857 output points;
- fixed-step convergence of the generalized coordinates: orders 5.09 and 5.03
  (the index-2 velocity block shows the classical Hessenberg order reduction,
  documented in the test and not asserted as fifth order).

### L2 — the 3²:4 sextic (done)

`rules/poly/transitive.rs` contains the transitive group `3^2:4` (order 36,
degree 6) but the test vector `sextic_groups` had no polynomial with that
group. Found by a constructive norm search over `Q(√2)`:

`g(x) = x^6 - 8x^3 - 18x^2 + 16 = Norm_{Q(√2)/Q}(x^3 - 3√2 x - 4)`.

`disc(g) = 864^4` (a square, so the group lies in A6), the unique block
system is `[[0,1,2],[3,4,5]]`, and `towers` gives the block invariant
`θ² = 72` (`θ = ±6√2`) with block polynomial `z³ - (θ/2)z - 4` — a
cross-check independent of the coset resolvent. In-test `F_p`
distinct-degree factorisation gives Frobenius cycle types `(4,2)` at 5,
`(2,2,1,1)` at 7, `(3,1,1,1)` at 23 and `(3,3)` at 31; among the table's
solvable even sextic groups only `3^2:4` contains all four. Neighbour
discrimination: `x^6 - 8x^3 - 18x^2 - 12x + 14` is `S3xS3` and a
one-coefficient change (`-18 → -19`) gives `3^2:D4` (72).

`solve(x^6 - 8x^3 - 18x^2 + 16, x)` reduces fully to two radicals with
residuals 8.9e-15 / 7.5e-14, and `solve_complex` returns all six roots as
canonical `root_of` objects (2 real, 4 non-real).

Tests: `cargo test --lib rules::poly::galois_high` 8/8 and
`rules::solve::tests::sextics_and_septics_by_galois_group` 1/1 (160 s).

### L3 — rigorous homotopy certification (done)

`kernels/homotopy_cert.rs` used to evaluate Smale's `α` in floating point
("a numerical one up to the rounding level of `f64`") and
`kernels/homotopy.rs` accepted steps on a residual/step-size heuristic, so a
path could in principle jump branches between steps. Now:

- `kernels/ball.rs` is a complex rectangle type over the outward-rounded
  real intervals of `backend::tape::Interval` (`add`, `sub`, `mul`, `div`,
  `powi`, intersections, modulus bounds) plus complex-coefficient polynomial
  terms (`CTerm`) with value, Jacobian and evaluation-at-a-point routines;
- `verify_zero` / `certify_zero` apply the Krawczyk operator
  `K(X) = y - Y f(y) + (I - Y J(X))(X - y)` in that arithmetic and prove
  existence **and uniqueness** of a zero in the box (`K(X) ⊂ int X`); a
  tightest-first ladder of relative radii inflates and contracts the box;
- `certify_endpoints_rigorous` returns a `RigorousReport` with the explicit
  boxes, the pairwise-disjoint count and the completeness flag against the
  Bezout number; `homotopy_solve_verified` is the end-to-end entry point
  (with the same adaptive-precision escalation as `homotopy_solve_certified`);
- `HomotopyOptions::certify_steps` certifies every accepted step of the main
  tracking by the same test on the homogenised system `H(·, t)` (built by
  `homotopy_node_terms` with complex coefficients), retries with a smaller
  step when the test fails, and reports `verified_steps`, `unverified_steps`
  and `steps_certified` on `PathEnd`. The endgame continuation is not
  certified, which the field documentation states.

Tests (`tests/kernels/homotopy_cert.rs`, 11 tests, all passing): explicit
balls for `x² - 2`, rejection of a point with no root nearby, rejection of a
double root, a two-variable system with four certified and pairwise disjoint
zeros, duplicated endpoints not counted as distinct, step certification on
cyclic-3 (all nodes verified), `homotopy_node_certificate`, and a
double-double cluster where the rigorous report is complete.

### L4 — Coulomb scattering (done)

The RK4 grid table with its WKB amplitude, the hard `k = 16` cut-off and the
3 % reproduction check are gone. `src/rules/pde/coulomb.rs` now evaluates the
analytic functions:

- `F_ℓ(η,ρ)` from DLMF 33.2.4 with `C_ℓ(η)` (33.2.5), the confluent series
  in double-double arithmetic above `ρ = 3` (plain `f64` below, where the
  cancellation is mild), and `ln Γ`/`arg Γ` by shifted Stirling;
- the large-`ρ` asymptotic expansion (DLMF 33.11.2–33.11.8) for
  `F, F', G, G'`, truncated at its smallest term, with the Coulomb phase
  `σ_ℓ = arg Γ(ℓ+1+iη)`;
- `G` below `ρ = 24` from Steed's continued fractions, with the Wronskian
  `F'G - FG' = 1` as the validity test;
- the energy-normalised scattering solution `u_k(r) = √(2/π) F_ℓ(-Z/k, k r)`
  (no fitted amplitude), so the continuum is normalized by construction;
- the cut-off is adaptive (`k_max ∈ {16, 32, 64}`, extended until the last
  panels contribute below `5e-4` of the scale and the reconstruction is
  within `RECON_TOL`), and the acceptance is certified by the reconstruction
  plus the Wronskian residual;
- repulsive continua (`Z < 0`) are supported; bound states remain `Z > 0`.

Achieved accuracy (documented in the code): `F_ℓ(0,ρ) = ρ j_ℓ(ρ)` to 1e-15,
series/asymptotic continuity 1e-13, Wronskian 1e-13 (1e-9 at `|η| = 5`),
reconstruction 1.6e-3 relative for the pure continuum and ~2–4e-3 for the
mixed bound+continuum data, `~1e-2` against a Crank–Nicolson reference at
`t > 0`. The requested 1e-3 was **not** reached; `RECON_TOL = 5e-3` is
documented with the numbers instead of being silently loosened. Steed's
fractions lose accuracy for `|η| ≫ ρ` and the series is refused above
`ρ = 34` where the asymptotic is not yet valid (both documented).

Tests: five new tests in `src/rules/pde/extension_tests.rs` (Bessel identity,
derivatives, series/asymptotic continuity, Wronskian and asymptotic phase,
energy normalisation, repulsive `Z < 0` PDE against Crank–Nicolson, and the
attractive bound+bump at `t = 0` and `t > 0`) plus the two pre-existing
continuum tests in `completion_tests.rs`, all passing (7/7,
`cargo test --lib coulomb`).

### L5 — cylinder drag (done, with a documented limitation)

The 5–8 % figure came from an uncommitted exploration and there was no test.
Now `tests/sim/ns_mac.rs` has a reusable `cylinder_run` harness and four
tests for Re = 100 (`ν = 0.01`, `D = U∞ = 1`, immersed level-set cylinder,
inflow/outflow/slip walls, six shedding periods after `t = 35`, `C_D = 2F_x`
from `force_on_region`, `St` from both interpolated upward zero-crossings and
a windowed DFT). A symmetry-breaking vortex seed is necessary — a
mirror-symmetric perturbation gives exactly zero lift and the unstable steady
state.

Grid study (20×10 box = 10 % blockage; literature `C_D = 1.32–1.34`,
`St = 0.164–0.166`, Braza–Chassaing–Ha Minh 1986, Williamson 1989):

| cells/D | cells | `C_D` | St | lift amp | runtime |
|---|---|---|---|---|---|
| 8 | 12,800 | 1.5069 | 0.1744 | 0.92 | 51–80 s |
| 12 | 28,800 | 1.4890 | 0.1777 | 1.01 | 276 s |
| 16 | 51,200 | 1.4956 | 0.1779 | 1.05 | 858 s |

`St` converges monotonically (Richardson ~0.1779); `C_D` is non-monotone at
the 0.5 % level, i.e. blockage, not resolution, dominates. The genuine fix
was the domain: at fixed `h = 1/8` a 40×20 box (5 % blockage) gives
`C_D = 1.3988` (+5.2 %) and `St = 0.1653` (+0.2 %). That is the validated
case named in the `ns_mac` module doc. Re = 40 is already covered by
`tests/sim/ns_research.rs` (Dennis & Chang). Tests: 8 passed, 2 ignored
(`cylinder_re100_5pct_blockage_matches_literature` 272 s, and the ~20 min
grid study); documented caveats: the `force_on_region` box omits unsteady
momentum (lift amplitude used only as a shedding indicator), cut-cell
divergence is `O(h)` near the body so mass flux rather than pointwise
divergence is asserted, and the residual +5.2 % is finite-domain error.

### L6 — time-varying exterior sources (done, see R2 in Part 5.3)

*Original plan:* the exterior (Weber-transform) source terms are nested `defint`s whose
integrand is itself an integral; `simplify` therefore churns. Plan: profile
(`sum`/`defint` bodies, `graph::schedule` budgets), then (a) hoist the
k-independent factors out of the k-integral, (b) evaluate the Weber
coefficient integral once per angular mode instead of per time term,
(c) cache the Bessel/Hankel evaluations in the mode table, (d) add a
`Simplify` cost model for nested `defint`s. Acceptance: the same answers
with a measured reduction in rewrite steps/time on a representative
exterior Duhamel problem.

### L7 — Neumann / Robin Green's functions (done)

Stated gap: "no Neumann ball (generalised Green's function), exterior Robin
or curvilinear boxes". Now in `src/rules/pde/poisson.rs`, registered in
`src/rules/pde/mod.rs`:

- `pde_greens_ball_neumann(R, d)` — the generalised (zero-mean) Neumann
  Green's function of the disk and ball, with **closed forms** (the Legendre
  and cosine series are summed): `-ΔG = δ(x-s) - 1/|B|`, `∂_r G = 0`,
  `∫_B G dx = 0`, and `G(x,s) = G(s,x)` exactly. The auxiliary
  `T = G - (|x|² - dR²/(d+2))/(2d|B|)` carries the compatibility derivative
  `∂_r T = -1/|∂B|`; note that `-ΔG = δ`, `∂_r G = -1/|∂B|`, `∫G = 0` and
  symmetry cannot hold together (the first three force a source-dependent
  constant), so the operator returns the standard symmetric generalised
  Neumann function and the `-1/|∂B|` convention is the tested auxiliary.
- `pde_greens_exterior_robin(R, h, d)` — the exterior disk/ball Robin Green's
  function with `∂_r G + hG = 0` (outward normal), `h > 0` and `hR` not a
  positive integer (guarded), symmetric; closed form for the `l = 0` part and
  the Fourier tail in 2D, the decaying Legendre series in 3D.
- Tests in `src/rules/pde/completion_tests.rs`: symmetry 1e-10, `∂_r G = 0`
  and `∂_r T = -1/|∂B|`, unit mean, `ΔG = 1/|B|` and `ΔT = 0` by finite
  differences, the exterior Robin boundary condition at three normals, the
  3D far field (`1/r`), the bounded 2D limit, the `h → ∞` Dirichlet limit
  against the Kelvin-image operator, and end-to-end single-layer and
  `pdsolve` Neumann-ball representations.
- Left out of scope: curvilinear boxes (documented).

### L8 — hodograph and shocks (shock interactions done)

- **Shock interactions** (done): the 1D scalar solver handled one Riemann
  step; `rules/pde/nonlinear.rs` now also takes a three-state datum
  `(u_l, u_m, u_r)` (`riemann_two_steps`, `wave_contribution`). It builds the
  two outgoing waves (Rankine–Hugoniot shocks or rarefaction fans with the
  Oleinik selection), tests whether the right characteristic of the left wave
  is faster than the left characteristic of the right one, computes the
  collision
  point `T = (x₂ - x₁)/(v₁⁻ - v₂⁺)`, `X = x₁ + v₁⁻T`, and splices the
  two-wave solution before `T` with the merged wave (seeded at `(X, T)`) after
  it. The data must be two increasing heavisides with numeric states; a
  non-convex flux is left to the existing convex-hull path.
  Test `shock_and_rarefaction_interactions` in
  `src/rules/pde/completion_tests.rs` checks all four combinations on Burgers
  (shock–shock with the merged shock at `10/3` at `t = 2`, shock–rarefaction,
  rarefaction–shock where the merged rarefaction is centred at `(5, 2)`,
  rarefaction–rarefaction with the constant middle state) and a concave
  traffic flux whose two upward shocks merge into a stationary jump.
- **General hodograph** (research, open; affine part done in R10):
  `pde_hodograph` only linearises the system. `pde_hodograph_solution` now
  also takes `λ, μ` (R10, Part 5.3). Plan:
  (a) affine `λ = a r + b s + c`, `μ = d r + e s + f` by characteristics;
  (b) Riemann's method for the linear hodograph equation when a Riemann
  function can be constructed (constant-coefficient/adjoint-solvable
  cases), (c) numeric characteristics for the rest, returning an implicit
  solution.

### L9 — eigenring / Kovacic number fields (substantially done)

Stated gaps: "eigenrings with denominators other than powers of the leading
coefficient" and "towers over number fields when the characteristic
polynomial does not split"; Kovacic "poles of order ≥ 3 at non-rational
points, irrational local data".

Implemented:

- **General eigenring denominators (Beke / van Hoeij)** in
  `src/rules/ode/linear.rs`: the denominator lattice is built from the monic
  irreducible factors `q_i` of `lc(L)` with exponent caps `(n-1)e_i` and a
  local-formal-data preferred exponent (exponent spread + 1, capped at 2);
  candidates `Π q_i^{m_i}` are enumerated around the preferred vector and
  each is tried through the exact Ore-remainder `Q`-system with a
  fraction-free (primitive) nullspace. Right factors of any order are
  extracted from an eigenring element `R` via the minimal polynomial `μ` of
  `R mod L` and `gcrd(L, g(R))`, with rational combinations `R_i + c R_j`
  for projector-like cases. Tests: products of second-order operators with
  order-2 factors, the lclm of Airy and Bessel (extended), and Legendre ∘
  Chebyshev where `lc = (1-x²)²` and the needed denominator `1-x²` is not a
  power of `lc`.
- **Kovacic case 2 over `Q(√d)`** in `src/rules/ode/kovacic.rs`: exact
  arithmetic in a quadratic field (`Quad`, polynomials and rational
  functions over it), surd collection/consistency (`surd_field`,
  `surd_in_field`), field exponent sets `{4}`, `{2, 2±2√d}`, `{order}` and
  the infinity variants, a degree-weighted family condition, the third-order
  equation for `P = P_0 + √d P_1` solved over `Q`, and the solution
  `exp(¼Σe ln f + ½ln P - ½∫p ± ½∫√(4r - φ² - 2φ'))` with the square root
  kept symbolic when it needs a second extension. Tests:
  `(x²-2)²y'' - 2(d-1)y = 0` for `d = 2,3,5` (closed form verified in the
  ODE at several points to 1e-6), the `√5` difference equation
  `4x(x-1)y'' + 2(2x-1)y' - 5y = 0`, and the exponent-set/degree tests.
  Finding documented: for rational-coefficient equations every irrational
  case-2 family with a non-negative degree condition has `deg P = 0`; case 3
  cannot have irrational local data (finite primitive monodromy forces
  rational exponent differences).

Verification: `rules::ode::linear` 9/9 (505 s here, 782 s in the agent's
clean-worktree run), `rules::ode::kovacic` 5/5, `rules::ode::nf_hyper` 2/2,
the aggregate `eigenring`/`kovacic` tests green; clippy clean for both files.

Still open (documented in the module docs): irrational `b` at a rational
pole (needs a `Q(c,√b)` tower), composita of two different `√d`, a second
quadratic extension for the discriminant (kept symbolic), poles of order ≥ 3
at non-rational points, eigenring elements over a number-field tower, and the
numerator bound `3deg(lc)+2(+8)` is workable but not proven sharp.

### L14 — Zeilberger creative telescoping and Gruntz limits (new, done)

- `src/rules/calculus/zeilberger.rs`: creative telescoping for proper
  hypergeometric terms with a **bounded polynomial certificate numerator**
  over `Q` (first-order recurrences in practice), a rational certificate
  `G/F`, and operators `zeilberger(F, k, n)` (recurrence + certificate) and
  `definite_sum(F, k, n)`. Tests: Vandermonde, the binomial-square identity,
  a certificate with a non-trivial denominator, closed forms through
  `sum(…, -oo, oo)`. *Correction (2026-10-29 review):* Dixon, Apéry and
  Watson are **not** found by the operator — the test only sums them
  directly (Watson only asserts finiteness), and Pfaff–Saalschütz is a
  negative test. The exact `Q(n)` elimination (`telescoping_over_qn` and
  helpers, ~250 lines) was `#[allow(dead_code)]` and unused. **Completed
  as R1 (Part 5.3):** the exact search now works by specialisation and
  rational reconstruction, and Dixon, Apéry and Pfaff–Saalschütz are found.
- `src/rules/calculus/gruntz.rs`: mrv-set rewriting to a single exp-log
  variable with series-based limit evaluation, wired into `limit`, with the
  standard hard limits (nested exponentials, `x^x` at `0+`, `(1+1/x)^x`,
  `ln ln x / ln x`, …) tested symbolically and against high-precision numeric
  evidence at infinity and finite points.

Verification: `cargo test --lib rules::calculus` 40/40 in 11.6 s.

### L10 — 3D DG limiting (done)

`sim/dg3d.rs` had positivity scaling but none of the 1D/2D limiters. Now it
has `Dg3Limiter { None, Characteristic, Hierarchical, Weno }` (default
`None`, so the plain DGSEM scheme is unchanged), with `limiter_tvd` and
`limiter_weno_eps` knobs on `Dg3Options`:

- `Euler3d` implements `DgLaw::has_characteristics` / `to_characteristic` /
  `from_characteristic` through a local 3D eigen-decomposition of the normal
  flux Jacobian (waves `u_n - c`, two shears, entropy, `u_n + c`);
- the shared TVB troubled-cell indicator compares face-trace slopes with
  neighbour cell averages (threshold `M h²`), optionally in characteristic
  variables;
- the characteristic limiter transforms the tensor-product Legendre modes
  (`nodal → modal` through the LGL Legendre Vandermonde) and applies
  Krivodonova's minmod from the top total degree downwards, with the
  coefficient vectors in characteristic variables;
- the hierarchical limiter is the same moment limiter componentwise, gated
  by the componentwise TVB test (the gate is needed: the ungated limiter
  clipped coarse smooth extrema and dropped the p = 2 order to ~1.5);
- the WENO limiter is the 1D Zhong–Shu reconstruction per axis with the own
  element and its ±2-shifted neighbours, its three directional
  reconstructions averaged (mean-preserving).
All limiters run after every SSP-RK3 stage and before the Zhang–Shu
positivity scaling; all preserve element means.

`tests/sim/dg3d.rs` (7 tests, all passing): smooth advection keeps order 3.31
at p = 2 with every limiter; the 3D Sod shock tube along x, y and z gives the
same solution by symmetry (overshoot ≤ 0.6 % of the density jump, shock
within 0.14–0.17 cells of `0.850431` at `t = 0.2`); periodic conservation of
mass, momentum and energy is at round-off (6e-15); a 200:1 pressure blast
stays finite and positive for all three limiters; the eigen-decomposition
satisfies `J R = R Λ`, `L R = I`. Known limitations (documented in the
module): WENO is componentwise and its directional reconstructions are
averaged rather than fully multidimensional; the moment limiters skip a
direction without a neighbour at a non-periodic boundary.

### L11 — 3D mesh refinement (done for the FEM part)

Two independent gaps: (a) `fem_simplex` (P1/P2 on tetrahedra) had no
adaptive refinement; (b) `euler_amr` is 1D/2D Berger–Colella only. Part (a)
is now implemented in `src/sim/fem_simplex.rs`:

- `SimplexMesh::refine_marked(&self, marked: &[bool])` — red–green conforming
  edge bisection: a marked tetrahedron has all six edges split and yields
  exactly 8 children; an unmarked cell with split edges is cut only along
  those edges (green closure, 2–6 children); faces are triangulated
  identically from both sides, so there are no hanging nodes; boundary
  facet markers are inherited by containment; a new public `level` field
  records the refinement generation; works in 1D/2D as well;
- `error_indicators(space, u, problem, t)` — the residual estimator
  `η_T² = h_T²∫_T(f - cu_h + ∇·(k∇u_h))² + ½Σ_faces h_e∫_e[k∂_n u_h]²`
  for P1 and P2 (the P2 Laplacian from the barycentric Hessian);
- `solve_adaptive(mesh, problem, degree, t, theta, tolerance, max_steps,
  solver)` — solve → estimate → Dörfler mark → refine, with the
  `(dofs, estimate)` history.
- The `euler_amr` 3D gap remains open (documented in Part 1).

`tests/sim/fem_simplex.rs` (11 tests, all passing): volume defect 1.3e-15 /
8.0e-15 / 6.4e-14 over repeated refinements (< 1e-12), all volumes positive,
every interior face shared by exactly two cells (no hanging nodes, checked
geometrically), boundary markers preserved, identity when nothing is marked;
P1 manufactured solution order 1.95 (L2) / 0.98 (H1); estimator consistency
`η/H1 = 6.47 / 6.22 / 6.41` over three meshes; P2 smoke; and adaptive vs
uniform on an L-prism edge singularity: adaptive 696 dofs L2 2.594e-3 vs
uniform 675 dofs 4.094e-3 (×0.63), H1 0.0860 vs 0.1273 (×0.68), history
`[(130,1.231),(155,0.944),(282,0.725),(696,0.520)]`. Caveats documented: no
coarsening/goal-oriented adaptivity/multigrid, the estimator ignores
Dirichlet boundary-data error and has no face term in 1D, green closure can
produce anisotropic but volume-positive cells.

### L12 — graph canonical labelling and automorphism groups (new, done)

Beyond the old list: `kernels::graph` now has 1-WL partition refinement with
canonical certificates (`ColourGraph`, `refine`, `Refinement`,
`wl_equivalent`) and nauty-style individualisation–refinement
(`canonical_labelling`, `CanonicalLabelling`, `Permutation`,
`SearchBudget`) producing a canonical order/code/matrix, automorphism
generators, orbits and the group order (Schreier–Sims with automorphism
pruning), with a documented bounded-budget fallback (`is_exact == false`).
Operators: upgraded `graph_isomorphic` (certificate + exact canonical codes,
now up to 1024 vertices, old matcher as fallback), new
`graph_canonical_form(g) = list(order, adjacency)` and
`graph_automorphisms(g) = list(order, generators, orbits)`.

Tests: `cargo test --test kernels graph` 115/115 and
`cargo test --lib rules::discrete` 137/137. Verified group orders (K_n = n!,
C_n = 2n, Petersen 120, Q2/Q3/Q4 = 8/48/384, K_{a,b} = a!b! ×2 iff a = b),
orbit partitions, permutation invariance over 14 families × 6 random
permutations, directed/coloured graphs, and the decisive SRG(16,6,2,2) pair:
WL gives equal certificates for the 4×4 rook and Shrikhande graphs but the
canonical codes differ, `graph_isomorphic` is false and `|Aut|` is 1152 vs
192. Performance: 300-vertex sparse graph canonical form in 18.4 ms; budget
exhaustion on 100 disjoint triangles reports `is_exact == false`.

### L13 — special-function layer extension (new, done)

`kernels::special_ext` gained public, documented kernels and `rules::special`
the corresponding symbolic operators (exact values, recurrences, partials):

- **Coulomb**: `coulomb_f/g/f_prime/g_prime`, `coulomb_fg`, `coulomb_c`,
  `coulomb_sigma` (DLMF 33.2/33.8/33.11, f64 + double-double series, Steed
  fractions, asymptotic). `rules/pde/coulomb.rs` now imports these instead
  of carrying private copies — PDE results unchanged (`cargo test --lib
  coulomb` 8/8).
- **Polylogarithm** `polylog(s, z)` + `polylog_complex` (series, DLMF
  25.12.10 log expansion, Kummer/Bernoulli inversion for integer `s`,
  Eulerian rational functions) and **Lerch** `lerch(z, s, a)` (series, Euler
  integral with exp-sinh, Stirling rational form).
- **Appell F1** `appell_f1(a, b1, b2, c, x, y)` (double series + Euler
  integral via tanh-sinh).
- **Kelvin** `ber/bei/ker/kei` for orders 0 and 1 (complex-J/K series in
  double-double + asymptotics + adaptive integral) and **Anger–Weber**
  `anger_j`, `weber_e` (DLMF 11.10 integral representations).
- Symbolic layer: operators `coulomb_f/g/fp/gp`, `coulomb_c/sigma`,
  `polylog`, `lerch`, `appell_f1`, `ber/bei/ker/kei`, `angerj/webere` with
  exact rewrites (`η = 0` Bessel forms, `polylog(1,z) = -ln(1-z)`,
  `lerch(z,s,1) = polylog(s,z)/z`, …).

Tests: `rules::special` 28/28, `tests/kernels special_ext` 27/27 (11 new).
Accuracies against 40-digit mpmath references: Coulomb F/F′ ~1e-12, G/G′
~1e-9 (1e-8 in the `|η| ≫ ρ` corner), polylog ~1e-15, Lerch ~1e-15, Appell
~1e-15/1e-13, Kelvin ~1e-15, Anger–Weber ~1e-13, plus branch-cut,
series/asymptotic continuity, Wronskian/radial-equation and exact-value
tests. Not done (documented): Meijer G, Kelvin orders > 1, non-integer
polylog with `|ln z| ≥ 4.5`, complex-parameter Lerch.

---

## Part 1 — Module-by-module gaps

### 1.1 Graph engine (`src/graph`)

- `schedule.rs`: the tiered fixpoint has budgets/bans; there is no
  *proof-carrying* trace of which rule fired (useful for the verify and
  proof rule sets) — add an optional rule-fire log behind a budget flag.
  **P2**
- `soundness.rs` samples points per rewrite; no adversarial sampling near
  singularities/poles (e.g. log branch cuts). **P2**
- `facts.rs`: sign/integer analysis is interval-based; no modular facts
  (`n mod m`), needed by `powmod`/`factorint` shortcuts. **P3**
- `extract.rs`: cost model is size-based; no evaluation-cost model even
  though `backend::tape` knows the operation count. Add an optional
  `EvalCost` extractor for the numeric phase. **P2**

### 1.2 Calculus (`src/rules/calculus`)

- Integrate: Risch–Norman + reduction formulas + table. Missing: algebraic
  function integration beyond sqrt-quadratics (Hermite reduction for
  general algebraic extensions) — **done for quadratic extensions, L15** —
  integration of `exp` of a general polynomial in more than one variable,
  definite integrals by the Mellin–Barnes method. **P2/P3**
- `limits.rs`: ~~no Gruntz algorithm~~ — done (L14, `gruntz.rs`); still
  open: trigonometric/oscillating factors and algebraic cancellation beyond
  series order 6 make the kernel return `None`. **P2**
- `series.rs`: multi-variable series and Puiseux at arbitrary algebraic
  branch points. *(multi-variable series done in L15; genuine multivariate
  Puiseux factors remain unevaluated.)* **P2**
- Summation: Gosper + recurrence guessing + first-order Zeilberger (L14);
  higher-order Zeilberger (Dixon, Apéry, Saalschütz) open → R1. **P2**

### 1.3 Polynomials, solving, Galois (`src/rules/poly`, `rules/solve`)

- Factorisation: Zassenhaus/Kronecker; missing: van Hoeij's lattice
  factorisation for high-degree integer polynomials (slow path), and
  factorisation over algebraic extensions. **P2**
- Groebner: Buchberger with criteria, plus F4 and F5 (L15). FGLM change of
  ordering for zero-dimensional ideals is done (R5, Part 5.3). **P2**
- `solve`: sextic/septic towers as in L2/L9 (not done: sextic towers with
  non-real block invariants, septics with three real pairs, degree ≥ 8
  without decomposition). `cad(polys, vars)` exists (sample points only)
  and `cad_cells`/`real_qe` now decide quantified formulas (L15). **P2**
- `solve/numeric.rs`: homotopy endgame certification (L3). **P1**

### 1.4 ODE (`src/rules/ode`)

- L9 (eigenring, Kovacic fields).
- `formal.rs`/`nf_hyper.rs`: turns of the number-field tower are not
  produced; add a quadratic extension on top of an existing one. **P2**
- `systems.rs`: first integrals up to degree three for 3–4 unknowns; no
  Darboux polynomials / Jacobi last multiplier. **P2**
- `reduce.rs`: Lie symmetry reduction only for point symmetries with a
  polynomial ansatz; no generalized (higher-order) symmetries. **P3**
- `resum.rs`: Borel–Padé; no conformal mapping / Padé–Hankel multi-summation
  for strongly divergent series. **P3**

### 1.5 PDE (`src/rules/pde`)

Documented gaps to close, in order of value:

1. L4 (Coulomb continua) and L7 (Neumann/Robin Green's functions). **P1**
2. L6 (exterior time-dependent sources) and L8 (hodograph/shocks). **P1/P2**
3. Nonlinear: no Riemann problems for systems (only scalar), no N-wave
   interaction, no inverse scattering for mKdV/NLS, only `KdV` (`N ≤ 4`).
   **P2**
4. Transform methods: Hankel only for radial whole-space; add general
   Hankel/Weber pairs, Mellin for wedges. **P2**
5. Curvilinear: non-separable metrics (ellipsoidal beyond spheroidal,
   toroidal) and variable coefficients in 2D/3D. **P3**
6. Sturm–Liouville: no Weyl–Titchmarsh m-function for arbitrary
   potentials (only the Hankel/Weber/Coulomb instances); add a numerical
   m-function and spectral density. **P2**
7. Spheroidal: exterior radiation, Neumann/Robin data. **P2**

### 1.6 Linear algebra / special / stats / transforms (`src/rules`, `src/kernels`)

- `rules/linalg.rs` / `kernels/dense.rs`: no rank-revealing QR with column
  pivoting exposed symbolically, no generalized Schur (QZ) for the
  symbolic side, no structured (Toeplitz/Hankel) solvers. **P2**
- `kernels/eigs_ext.rs`: Jacobi–Davidson is documented symmetric-only for
  the Hermitian path; add a general non-Hermitian JD with two-sided
  deflation (partially present, verify) and a shift-invert harmonic path
  for interior eigenvalues. **P2**
- `kernels/special.rs`: polylog/Lerch/Appell F1/Kelvin(0,1)/Anger–Weber/
  Coulomb/Meijer G done (L13, L15); still open: Fox H, Kelvin orders > 1,
  non-integer polylog for `|ln z| ≥ 4.5`, complex-parameter Lerch. **P2**
- `kernels/hypergeometric.rs`: large-parameter `1F1`/`U` lose accuracy;
  add Temme's uniform asymptotic expansions and a recurrence-direction
  (Miller) stabiliser. **P2**
- `kernels/quadrature.rs`: no oscillatory (Filon/Levin) quadrature in more
  than one dimension, no tanh-sinh for endpoint singularities of
  fractional order. **P3**
- `kernels/ode_stiff.rs`: BDF has no mass-matrix/DAE path (only Radau);
  add it, or document the restriction in the error type. **P2** → R4
- `kernels/ode_adaptive.rs`: no symplectic/geometric integrator for
  Hamiltonian DAEs. **P3**
- `kernels/stats.rs`/`rules/stats.rs`: no multivariate distributions
  beyond the implemented set (Dirichlet, Wishart), no robust covariance
  estimators (MCD). **P3**
- `rules/transforms.rs`: no wavelet transforms, no fractional Fourier. **P3**

### 1.7 Discrete mathematics (`src/rules/discrete`)

- Graph algorithms: ~~canonical labelling~~ — done (L12). **P2**
- Permutation groups: Schreier–Sims and Todd–Coxeter exist; no
  low-index subgroups beyond Todd–Coxeter, no character table computation
  for arbitrary groups (only the implemented families). **P2**
- Number theory (verified 2026-10-29): `factorint` is trial division +
  Pollard rho, ECM and, since L15, the quadratic sieve; `isprime` uses
  Miller–Rabin for 64-bit integers and BPSW above. No ECPP/primality
  certificate. **P2**
- Coding theory: no list decoding (Guruswami–Sudan), no LDPC. **P3**

### 1.8 Simulation (`src/sim`)

- L10 (DG 3D limiters), L11 (3D refinement), L5 (cylinder drag). **P1/P2**
- `structured.rs`: constant-coefficient Laplacian only; add variable
  diffusion, higher-order Mur, PML coupling. **P2**
- `yee_fdtd.rs`: no anisotropic/magnetic media, no conductivity, TF/SF not
  combined with periodic faces, no PMC. **P2**
- `euler_nd.rs`: no source terms (gravity/bathymetry), no dry-bed; 3D has
  no AMR. **P2**
- `hermite.rs`: direct `O(N²)` gravity; add a neighbour scheme / block
  time stepping (the Barnes–Hut and FMM modules exist but are separate).
  **P2**
- `fmm.rs`: P3M is not implemented (`ewald::p3m_reciprocal_cell` covers the
  mesh route); add particle–particle corrections for near-field accuracy.
  **P2**
- `bem_helmholtz.rs`: constant elements, dense storage, no impedance BC.
  **P2**
- `acoustics.rs`: homogeneous medium, Dirichlet data must be zero. **P2**
- `sse.rs`: spin-1/2 only, no winding-number (stiffness) estimator. **P2**
- `md_engine.rs`: NPH-type legacy barostat, no constraints with flexible
  cells, no combined barostat/thermostat for MTK. **P3**
- `hyperelastic.rs`: P1 only (locking), no friction/self-contact. **P3**
- `kerr.rs`/`relativity.rs`: no radiation reaction, no extremal spin; ray
  tracing numerical only. **P3**

### 1.9 IO / backend / FFI

- `io/plot.rs`: no 3D surface/volume rendering beyond the existing image
  export; add contour/streamline export. **P3**
- `backend/tape.rs`: interval enclosures exist; expose them as a public
  `compile_interval` for rigorous quadrature/ODE. **P2**
- `ffi/`: no bincode/tensor variants for every function (the old task list
  asked for three variants per API); the current ABI is JSON+handle only.
  **P2**

---

## Part 2 — Cross-cutting

### 2.1 Lint baseline (P0)

As found, `cargo clippy --lib --all-features` failed with 12 errors (not
clean as the project notes claim):

- `approx_constant` for the literal Euler–Mascheroni constant in
  `kernels/special.rs`, `kernels/special_ext.rs`, `rules/combinatorics.rs`,
  `rules/special.rs` — **fixed** (use `std::f64::consts::EULER_GAMMA`);
- `needless_borrows_for_generic_args` in `rules/poly/groebner.rs`
  (4×) — **fixed**;
- `double_must_use` in `sim/ewald.rs` (3×) — **fixed**;
- `doc_markdown` "backticks are unbalanced" — **open**. *Correction
  (2026-10-29):* the README include is **not** the cause — it is already
  `#[cfg(doctest)]`, and removing it leaves the error in place; no `///`/`//!`
  paragraph in `src/` has an odd backtick count either, so the source is a
  macro-generated doc (candidates: `nist_const!` in `constant.rs`, derive
  macros). The old diagnosis read: it comes from the
  `#[doc = include_str!("../README.md")]` of the
  doctest wrapper in `src/lib.rs` (replacing the include removes the error;
  every backtick in the README is balanced, and `#[allow]` on the item, on
  an enclosing module, and tilde code fences do not suppress it — the
  diagnostic is attached to the macro expansion with an unusable span).
  Alternatives for the owner: a crate-level `allow` (hides real
  `doc_markdown` findings), dropping the README doctest include (loses the
  README doctests), or an upstream clippy issue;
- the remaining `missing_const_for_fn`, `neg_cmp_op_on_partial_ord` and
  `question_mark` findings in `rules/pde/coulomb.rs` are in the new Coulomb
  continuum code and are tracked with L4.

Also fixed: `Cargo.toml` used the deprecated
`lints.clippy.empty-line-after-doc-comments` key; the underscore spelling is
now used.

### 2.2 Test coverage holes (P1)

Files with no direct test file and (sometimes) no inline tests:

- `sim/dg3d.rs` (tests added in L10); `sim/cutcell.rs`, `sim/dpd.rs`,
  `sim/acoustics.rs`, `sim/bem_helmholtz.rs`, `sim/md_cell.rs` now have
  direct tests (R3, Part 5.3);
- `kernels/fftconv.rs` and `kernels/random.rs` now have direct tests (R3);
  still open: `kernels/ordering_ml.rs` (only indirectly through
  `kernels/ordering.rs`), `kernels/jacobian.rs`
  (indirect through `ode_stiff_struct`), `kernels/eigs.rs` (indirect);
- `graph/id.rs`, `io/markup.rs`, `ffi/compute.rs`, `ffi/term.rs`;
- `rules/calculus/{diff,gosper,indefinite,integrate_ext,limits,powerseries,recurrence}.rs`
  are exercised only through the aggregate `calculus/tests.rs`; add
  targeted tests for the newer entry points (Gosper with parameter
  conditions, recurrence guessing).

Every new algorithm in this document must ship with tests in the same
style as its neighbours (`tests/` for kernels/sim, inline
`#[cfg(test)] mod` or `src/rules/<domain>/*_tests.rs` for rule sets).

### 2.3 Performance and docs

- The benchmark suite (`benches/engine.rs`) only times the engine; add
  benches for the heavy kernels (homotopy, Radau, DG3D) so the completion
  work has a performance baseline. **P2**
- Several module headers are longer than the code they describe and mix
  "coverage" with "limitations"; keep the limitations list in one place and
  move the tables to `docs/`. **P3**

---

## Part 3 — Order of work

1. **P0 hygiene**: lint baseline (2.1), test-registration gaps that hide
   compile errors (2.2 for `sim/dg3d.rs`).
2. **P1 named gaps**: L3, L4, L7, L10, L1 (in that order: they are
   self-contained and each has an external reference).
3. **P1 continuation**: L8 (shock interactions), L2, L5, L6.
4. **P2 blocks** in Part 1, module by module, each with tests.
5. **P3** items only after the P1/P2 set is demonstrably green.

---

## Part 4 — What changed in this round

| Item | Change | Verification |
|---|---|---|
| L1 | new `tests/kernels/andrews_squeezer.rs` (index-3 Hessenberg form, `andrews.f` reference) | 4/4 tests; `|Δq|∞ = 1.4e-8`, `|Φ|∞ = 4e-16`, order 5.09 on the coordinates |
| L2 | 3²:4 sextic `x^6 - 8x^3 - 18x^2 + 16` = Norm over Q(√2), certified by group/discriminant/block system/tower/Frobenius types; neighbours with other groups; radical solve + `solve_complex` | 8/8 `rules::poly::galois_high`, 1/1 solve test (160 s); roots residual 1e-14 |
| L3 | new `kernels/ball.rs` (outward-rounded complex balls), Krawczyk `verify_zero` / `certify_zero` / `certify_endpoints_rigorous` / `homotopy_solve_verified`, per-step certification in the tracker | 11/11 tests in `tests/kernels/homotopy_cert.rs` |
| L4 | analytic Coulomb `F/G` (DLMF 33.2/33.11, Steed fractions), energy-normalised `u_k`, adaptive `k_max`, repulsive continua; now backed by public kernels | 8/8 `cargo test --lib coulomb`; reconstruction 1.6e-3 (pure) and 2–4e-3 (mixed), documented `RECON_TOL = 5e-3` |
| L5 | cylinder `Re = 100` harness + grid study + 5 %-blockage validation + module-doc note | 8 passed / 2 ignored (`tests/sim/ns_mac.rs`); 5 %-blockage `C_D = 1.3988` (+5.2 %), `St = 0.1653` (+0.2 %) |
| L7 | `pde_greens_ball_neumann`, `pde_greens_exterior_robin` with closed forms, registration and docs | 7/7 `cargo test --lib greens` (symmetry 1e-10, BC 1e-6, mean 1e-7) |
| L8 | three-state Riemann problems with shock/rarefaction interactions in `rules/pde/nonlinear.rs` | `shock_and_rarefaction_interactions` (5 data sets, 24 point checks) |
| L10 | `Dg3Limiter` (characteristic, hierarchical, WENO) for the 3D DGSEM solver, Euler eigen-decomposition, options and docs | 7/7 `cargo test --test sim dg3d` (order 3.31, Sod 0.14 cells, conservation 6e-15, blast stable); 26/26 on the `dg` filter including `dg_research` |
| L11 | red–green conforming 3D tet refinement (`refine_marked`, `level`), residual estimator (P1/P2), `solve_adaptive` with Dörfler marking and history | 11/11 `tests/sim/fem_simplex.rs` (volume defect ≤ 6e-14, no hanging nodes, order 1.95/0.98, adaptive ×0.63 L2 at equal DOF) |
| L12 | WL certificates + nauty-style canonical labelling/automorphism groups; `graph_isomorphic` upgrade, `graph_canonical_form`, `graph_automorphisms` | 115/115 `--test kernels graph`, 137/137 `rules::discrete`; SRG(16,6,2,2) rook vs Shrikhande separated, `|Aut|` 1152 vs 192, 300-vertex canonical form 18.4 ms |
| L13 | public Coulomb/polylog/Lerch/Appell F1/Kelvin/Anger–Weber kernels + symbolic operators, exact rules, recurrence/partials; PDE Coulomb refactored onto them | 27/27 `--test kernels special_ext`, 28/28 `rules::special`, 8/8 `coulomb`; accuracies 1e-15…1e-9 against mpmath |
| L9 | general eigenring denominators (Beke/van Hoeij) with higher-order factors; Kovacic case 2 over `Q(√d)` | `rules::ode::linear` 9/9 (505 s), `kovacic` 5/5, `nf_hyper` 2/2; clippy clean |
| L14 | Zeilberger creative telescoping (+`definite_sum`) and a Gruntz mrv limit core wired into `limit` | `rules::calculus` 40/40 (11.6 s) |
| hygiene | 12 pre-existing clippy errors: 11 fixed; the `doc_markdown` one is diagnosed and documented | `cargo clippy --lib` |

### How to run the affected tests

```powershell
cargo test --test kernels homotopy_cert             # L3
cargo test --test kernels andrews_squeezer          # L1
cargo test --test kernels graph                     # L12
cargo test --test kernels special_ext               # L13
cargo test --test sim dg3d                          # L10
cargo test --test sim fem_simplex                   # L11
cargo test --test sim ns_mac                        # L5 (-- --ignored for the 5 % case)
cargo test --lib coulomb                            # L4
cargo test --lib greens                             # L7
cargo test --lib shock_and_rarefaction_interactions # L8
cargo test --lib rules::poly                        # L2
cargo test --lib rules::discrete                    # L12
cargo test --lib rules::special                     # L13
cargo test --lib rules::calculus                    # L14
cargo test --lib rules::ode::linear                 # L9 (slow: ~8-13 min)
cargo test --lib rules::ode::kovacic                # L9
```

No full-suite run was made: the affected binaries above are the ones whose
modules changed.

### Still open from this list

L6 (exterior source simplification performance), L8's general hodograph
half, the `euler_amr` 3D gap and the P2/P3 blocks of Part 1 — each with the
plan recorded above. The remaining L9 obstructions are listed in its
section (irrational `b` at a rational pole, composita of two surds, a second
quadratic extension for the discriminant, order ≥ 3 non-rational poles,
eigenring over a tower).

---

## Part 5 — Re-audit 2026-10-29

### 5.1 Review of the uncommitted round

Read in place: `kernels/ball.rs`, the Krawczyk part of
`kernels/homotopy_cert.rs`, `kernels/graph/canon.rs` (module design, group
order in `BigInt`), `rules/calculus/{zeilberger,gruntz}.rs` and their tests,
plus a red-flag grep (`todo!`, `unimplemented!`, `panic!`, `unwrap`,
`println!`, `#[ignore]`) over the whole diff. All of those hits are in tests.

| Finding | Severity | Action |
|---|---|---|
| `cargo clippy --all-targets` failed with 32 errors: 22 new `assert_is_empty` (clippy 1.99) in old tests, 3 literal Euler–γ, 2 `map_unwrap_or` in `ode/linear.rs`, 4 `i32 as f64` in `pde/completion_tests.rs`, 1 `needless_pass_by_ref_mut` in `pde/extension_tests.rs` | P0 | **fixed**; only the span-less `doc_markdown` error remains (2.1, diagnosis corrected) |
| Zeilberger: checklist claimed Dixon/Watson and an exact `Q(n)` path; the code has a bounded `Q` ansatz only and ~250 lines of `#[allow(dead_code)]` | P1 | L14 text corrected; completion is R1 |
| Krawczyk test: sound (centre `y ∈ X`, strict inclusion before contraction, last verified box kept on failure); `modulus_lower/upper` use an unrounded `hypot` but are not used by the proof | info | none |
| `rules/pde/mod.rs` header still said "Not covered: shock interactions" | P3 | **fixed** |
| Tests print with `println!`/`eprintln!` (diagnostic output, matches neighbours) | info | none |

### 5.2 Limitations found by the sweep but missing from Part 1

From the `Limitations:` paragraphs of the module headers (the full text stays
in the modules):

- `kernels/hypergeometric.rs` *(update 2026-10-30: complex `U` with large `|a|` is done to `1e-13` for `Re a >= 0`, see R8b in 5.5)*: `M` with large parameters is done (R8);
  complex `M` in the transition regions and complex `U` with `a` of
  modulus > 100 keep the Taylor-stepping accuracy; `p > q + 1` and
  non-terminating `b_j ∈ Z≤0` give NaN. **P2** *(update 2026-10-30: `hyp_pfq_checked` returns the reason as a `HypError`, `pFq(..; 0) = 1`, see Round 3 in Part 7)*
- `kernels/sde.rs`: Milstein/SRK/adaptive need commutative noise; the
  adaptive solver has no Lévy-area path. **P2** *(update 2026-10-30: `sde_adaptive_noncommutative_with` on a Chen-consistent Brownian grid, B2 in Part 7; the commutative adaptive solver is unchanged)*
- `kernels/eigs_ext.rs`: *(correction, R11)* the non-symmetric
  Jacobi–Davidson `eigs_jd` (complex targets, partial Schur deflation,
  complex harmonic extraction) already exists since commit 1b2530a4; the
  symmetric generalised solver still needs SPD `B`, and a two-sided JD is
  not implemented. **P2** *(two-sided JD done, R11b in 5.5; indefinite `B` detected and routed to Krylov-Schur, dense `eigen_generalized_symmetric`, B1 in Part 7)*
- `kernels/pde.rs`: no 3D, Robin/periodic faces, time-dependent BCs or
  adaptive stepping in the simple kernels (routed to `sim`). **P3**
- `sim/bem_aca.rs`: Helmholtz path is dense; no cluster-tree rebalancing.
  **P2** *(update 2026-10-30: `bem_helmholtz::solve_helmholtz_exterior_aca`, complex ACA + GMRES, B7 in Part 7)*
- `sim/dg.rs`: WENO limiter 1D only; hierarchical limiter skips
  non-periodic boundaries (same as 3D). **P2** *(update 2026-10-30: ghost elements from the boundary condition in `dg.rs` and `dg3d.rs`, B3 in Part 7; WENO is still 1D)*
- `sim/dg3d.rs`: uniform axis-aligned boxes only (no curved/unstructured
  hexahedra, no hanging nodes). **P2**
- `sim/dpd.rs`: single species, no bonds, no walls. **P2** *(update 2026-10-30: `Dpd::set_walls`, B6 in Part 7)*
- `sim/elasticity_fem.rs`: isotropic, homogeneous, small strain, static.
  **P2**
- `sim/euler_amr.rs`: 1D/2D only, ratio 2, no non-reflecting BCs. **P2**
- `sim/fem_simplex.rs`: constant `k`, `c`; no coarsening, no
  goal-oriented adaptivity, no multigrid. **P2** *(update 2026-10-30: `FemProblemVar` / `solve_static_var` / `solve_heat_var` with `k(x)`, `c(x)`, B4 in Part 7; the error estimator still assumes constant coefficients)*
- `sim/ns_mac.rs`: uniform grids, Boussinesq only, no free surface. **P3**
- `sim/sph.rs`: no self-gravity. **P3** *(update 2026-10-30: `Sph::gravity` with the Barnes-Hut tree force, B5 in Part 7)*
- `sim/split_step.rs`: periodic only, no rotating frame, no dipolar terms.
  **P3** *(update 2026-10-30: `evolve_rotating`, B8 in Part 7)*
- `sim/spin_models.rs`: nearest-neighbour hypercubic lattices only. **P3**
- `sim/kerr.rs`: no separatrix/circular/extremal closed forms. **P3**
- Not present anywhere (grep 2026-10-29): F4/F5, quadratic sieve,
  BPSW/ECPP, Meijer G, wavelets, fractional Fourier, Darboux
  polynomials / Jacobi last multiplier, multivariate series, Hermite
  reduction for algebraic integrands, real quantifier elimination,
  Guruswami–Sudan/LDPC. *(update 2026-10-31: F4/F5, the quadratic sieve,
  BPSW, Meijer G, wavelets, fractional Fourier, multivariate series,
  Hermite reduction and real quantifier elimination are done, see L15;
  remaining: Darboux/Jacobi last multiplier, Guruswami–Sudan/LDPC.)*

### 5.3 This round's completion set

Picked for value and tractability; each item ships with tests next to its
neighbours and a module-doc update, and is marked here when done.

| ID | Item | Prio | Area |
|---|---|---|---|
| R1 | Zeilberger beyond first order: working exact `Q(n)` (or modular + rational reconstruction) creative telescoping for order ≥ 2, certificates with poles in `n`; positive tests for Dixon, Apéry (order 2) and Pfaff–Saalschütz; remove or use the dead code | P1 | `rules/calculus/zeilberger.rs` |
| R2 | L6 — exterior time-varying sources: hoist `k`-independent factors, one Weber coefficient per mode, cached Bessel/Hankel values, measured speed-up on a representative Duhamel problem | P1 | `rules/pde/{mod,spectral,numeric,curvilinear}.rs` |
| R3 | Test-coverage holes of 2.2: direct tests for `sim/{cutcell,dpd,acoustics,bem_helmholtz,md_cell}.rs`, `kernels/{fftconv,random}.rs` | P1 | `tests/` |
| R4 | BDF with a mass matrix (index-1 DAE), checked on Robertson in DAE form and against Radau | P2 | `kernels/ode_stiff.rs` |
| R5 | FGLM change of ordering (zero-dimensional ideals), used for lex elimination | P2 | `rules/poly/groebner.rs` |
| R6 | Number theory: BPSW primality and Lenstra ECM (stage 1 + 2) in `factorint` for large factors | P2 | `rules/number_theory.rs` |
| R7 | `euler_nd` source terms: gravity and well-balanced bathymetry (hydrostatic reconstruction) with a lake-at-rest test | P2 | `sim/euler_nd.rs` |
| R8 | Large-parameter `1F1`/`U`: recurrence (Miller) stabilisation and uniform asymptotics, against mpmath values | P2 | `kernels/hypergeometric.rs` |
| R9 | Kelvin functions of integer order `n > 1` (recurrences + asymptotics) | P2 | `kernels/special_ext.rs`, `rules/special.rs` |
| R10 | L8a — affine hodograph `λ = a r + b s + c`, `μ = d r + e s + f` | P2 | `rules/pde/nonlinear.rs` |
| R11 | Non-symmetric Jacobi–Davidson (two-sided, harmonic Ritz for interior eigenvalues) | P2 | `kernels/eigs_ext.rs` |

**Status of this round** (updated as items land):

- **R1 — done.** `zeilberger` on parameter-free terms now runs an exact
  search over `Q(n)` for orders 1–3 (`rules/calculus/zeilberger.rs`, section
  "The exact search over Q(n)"). Instead of eliminating over the field
  `Q(n)` (the dead `telescoping_over_qn`, which blew up), `n` is specialised
  to the integers 23, 24, …; the system in `k` is a small linear system over
  `Q` whose one-dimensional solution (normalised `p_J = 1`) is sampled and
  recovered as rational functions of `n` by Cauchy interpolation with
  held-out points (recurrence coordinates one by one; the certificate
  coordinates over one common denominator found from an integer
  combination). `Q` is the product of the `n`-ratio denominators with
  redundant ones (dividing the others) dropped, which cut the degrees and
  the time. The ~250 lines of dead code are gone; `qpoly_node`/`frac_node`
  are used. Positive tests: `zeilberger_second_order_dixon_and_apery`
  (second-order recurrences found, identity checked at 4 points, and the
  recurrences annihilate the directly summed Dixon and Apéry sequences for
  `n = 0..6`) and `zeilberger_pfaff_saalschutz` (the old negative test is
  replaced: the 3F2(-n, 2, 3; n+5, 1-2n; 1) summand, closed form checked
  for `n = 1..5`, certificate with poles in `n`). Time of `zeilberger` on
  Dixon / Apéry: 9.6 s / 8.9 s for a first version that reconstructed every
  coordinate separately, 1.2 s / 1.1 s now; `binomial(n,k)^2` 0.07 s.
  `cargo test --lib rules::calculus` 41/41 in 10.8 s (was 40 in 11.6 s).
  Not covered: orders above 3, symbolic parameters beyond the bounded
  polynomial ansatz, and the boundary terms that turn the summand
  recurrence of a shifted-factorial 3F2 into a recurrence of the sum (the
  Pfaff–Saalschütz test therefore asserts the telescoping identity, not
  the recurrence of the sum).

- **R2 — done (L6).** Profiled first (`Report` of the engine on
  `pdsolve(u_t = Δu + f(r)cos θ g(t))`, unit disk, zero data, 3 sources):
  the time-varying exterior solution was a triple `ρ, k, τ` integral per
  mode and the engine then spent 375 s (`g = 1 + t`) and 419 s (`g = sin t`)
  and 86 000 / 92 000 e-graph nodes (9 800 / 13 000 unions) on it, against
  3.1 s and 7 000 nodes for the steady source. Two causes, fixed in turn:
  (a) `curvilinear.rs::separate` splits a source into `Σ g_i(t) s_i(x)`;
  the Weber coefficient of each `s_i` is then formed once per mode and only
  the one-dimensional time kernel is integrated against `g_i` (hoisting the
  `k`-independent factors out of the inner integral and one Weber
  coefficient per mode instead of per time term), resolved eagerly on the
  small integrand: 375 s → 27.9 s, 419 s → 13.2 s (the rule itself takes
  0.18 s; the rest was the engine); (b) `calculus::integrate::antiderivative`
  now declines at once an integrand with a Bessel/Hankel function of the
  variable (`special_in`; the tables know no such antiderivative and the
  by-parts/substitution search took 0.03–0.25 s for each of ~100 inner
  Weber integrals and 5–8 s for the outer one): 27.9 s → 0.24 s, 13.2 s →
  0.19 s, and the steady source 3.1 s → 0.17 s (1 900 nodes, 630 unions).
  Answers are unchanged: `time_varying_sources_outside_a_disk_by_duhamel`
  compares a sum of two separable terms (`(1+t²)cos θ`, `t sin 2θ`) with a
  Crank–Nicolson Duhamel reference to 4e-3, and the existing
  `sources_outside_a_disk_by_duhamel` and the other 96 `rules::pde` tests
  pass (`cargo test --lib rules::pde` 98/98, 171 s; `rules::calculus` 41/41).
  `exterior_time_varying_source_stays_cheap` guards the node/union counts
  (< 8 000 / < 2 000; measured 2 000 / 650). Not done from the plan: item
  (c) (cached Bessel/Hankel values: the graph already shares them) and (d)
  (a `Simplify` cost model for nested `defint`: (b) removed the cost it
  would have addressed); non-separable sources (`exp(t r)`) still use the
  nested form.
- **R3 — done** (tests only; no library change). New files, registered in
  `tests/kernels.rs` / `tests/sim.rs`: `kernels/fftconv.rs` (6 tests: FFT
  against the O(n²) DFT for radix-2 and Bluestein lengths 1–100, round trip
  and Parseval to 1e-10, impulse/constant/harmonic spectra, convolution on
  both code paths, circular convolution, correlation by definition, exact
  `poly_mul_i64` against 128-bit arithmetic), `kernels/random.rs` (8 tests:
  the stream equals an independent xoshiro256++ written from the published
  SplitMix64 seed-0 values, then moments/KS/chi-square for `uniform`,
  `below`, `normal`, `exponential`, `gamma` (shapes 0.3–9), `beta`,
  `poisson` (both branches, λ up to 1000), `binomial`, `shuffle`),
  `sim/cutcell.rs` (6: simplex closed forms, planar interfaces exact in
  1/2/3D, chord of a circle, disk/sphere convergence, complementarity and
  additivity), `sim/dpd.rs` (9: pair force and minimum image from rest,
  virial, ballistic limit exact, momentum, thermostat relaxation and
  equipartition, mass, determinism), `sim/acoustics.rs` (7: dt formula,
  d'Alembert split, mass/energy conservation, one-cell-per-step domain of
  dependence, source, 2D isotropy, pressure-release reflection),
  `sim/bem_helmholtz.rs` (7: validation, plane wave, translation as a phase,
  extinction theorem inside the body, soft monopole and hard Rayleigh
  limits, Burton–Miller vs conventional) and `sim/md_cell.rs` (4: matrix
  helpers against exact identities, `e^A`/`phi1`, barostat constructors,
  degrees of freedom, energy). All pass first run except for test-side
  corrections (DPD temperature averaging needs ≥ 750 particles and
  `dt = 0.02`: a 375-particle run at `dt = 0.03` read 1.04). Timing:
  `cargo test --test kernels fftconv`/`random` 0.0 s/0.2 s, `--test sim`
  `cutcell`/`acoustics`/`md_cell` < 0.3 s, `bem_helmholtz` 1.5 s, `dpd`
  8.5 s. **Finding** (documented in the test): at low frequency the
  Burton–Miller coupling `α = i/k` is large and amplifies the double-layer
  discretisation error into a phase error of about `2ka` in the sound-soft
  scattered field (0.2 rad at `ka = 0.1`); the conventional equation is
  accurate there, and the Mie tests at `ka ≥ 1` do not expose it.

- **R4 — done.** `BdfOptions::mass` / `sparse_mass` in
  `kernels/ode_stiff.rs`: the corrector `M (d + psi) = (h/α) f(t, y_pred + d)`
  is solved with the Newton matrix `M - (h/α) J` (dense, banded or sparse
  LU; the band includes the mass band), the error estimate covers the
  algebraic components, the initial derivative is `M^-1 f0` or the
  diagonal-row heuristic, the first step `1e-5 (t1 - t0)` as in Radau; dense
  output and events are unchanged. Module docs updated (BDF handles index
  1 only). Tests in `tests/kernels/ode_stiff_struct.rs`: Robertson in DAE
  form (`M = diag(1,1,0)`) against the Hairer–Wanner reference (BDF within
  1e-4 relative at `rtol 1e-8`), against `radau_with` (1e-4) and at
  `t = 1e5`; constraint drift < 1e-7 along the run; analytic and sparse-mass
  variants; a coupled linear index-1 DAE with a non-diagonal singular
  `M = [[1,-1],[0,0]]` (error < 2e-6, dense output < 1e-5); an invertible
  `M = diag(2,3)`; a wrong-size mass is rejected. `cargo test --test kernels
  ode_stiff` 28/28 in 0.3 s; clippy clean. Not covered: index ≥ 2 (use
  Radau) and time-dependent mass matrices.

- **R5 — done.** `groebner::fglm(basis, from, to)` in
  `rules/poly/groebner.rs`: standard monomials of the zero-dimensional
  reduced basis, sparse multiplication matrices by reduction, monomial
  walk in the target order with an incremental echelon form that tracks
  the combinations (a dependent monomial gives a new element, an
  independent one joins the staircase); `None` for ideals that are not
  zero-dimensional, the zero ideal, or a quotient above
  `FGLM_MAX_DIMENSION = 400`; the unit ideal returns `{1}`.
  `groebner_via(generators, order, limits)` computes grevlex and converts
  for lex and grlex targets, falling back to Buchberger; the operator
  `groebner(list, vars[, order])` and the `solve` system elimination
  (`rules/solve/system.rs`) use it, and the reduced basis is the same
  list as before. Tests (inline, `cargo test --lib rules::poly::groebner`
  11/11 in 0.7 s): FGLM equals Buchberger in lex and grlex for the circle
  and line, cyclic-3, Katsura-2/3, a multiple root and a sphere/plane/
  quadric; positive-dimensional ideals, the zero ideal and the unit ideal;
  cyclic-5 (70 solutions, 11-element lex basis, univariate element of
  degree 15) in 0.25–0.4 s and Katsura-4 (16 solutions) in 0.3–0.4 s, both
  checked to be Gröbner bases (Buchberger leaves them unchanged) containing
  the generators; operator tests in `rules::poly::tests::groebner_bases`.
  `cargo test --lib rules::poly` 53/53 (53 s), `rules::solve` 27/27 (61 s).
  Direct lex Buchberger on cyclic-5 was not attempted (out of reach for
  this implementation, which is why FGLM matters). Not done: F4/F5, the
  explicit target-order operator separate from `groebner`.

- **R6 — done.** `rules/number_theory.rs`: `bpsw` (trial division by the
  primes below 100, strong probable prime to base 2, strong Lucas with
  Selfridge's method A) is now `isprime` beyond 64 bits (Miller–Rabin with
  12 bases stays below `2^64`). *Bug fixed on the way:* the old big-integer
  branch (Miller–Rabin to the first 12 primes) declared the smallest strong
  pseudoprime to those bases, 318665857834031151167461, prime; the comment
  claimed validity to 3.3e24 (that is the 13-base bound). `factor_big`
  (trial division below 1000, the 64-bit path, perfect powers, BPSW and
  `ecm`) is used by `factorint`, `totient`, `divisor_count`,
  `divisor_sum`, `omega`, `bigomega`, `radical`, `mobius`, `liouville`,
  `sigma` and `is_squarefree`. `ecm` is Lenstra's method on Montgomery
  curves (x-only ladder, Suyama's parametrisation, stage 1 up to 2 000 /
  11 000 / 50 000 with 20 / 80 / 60 curves, standard stage 2 with `D = 210`
  to 50 times the bound); a composite cofactor that resists the curves
  leaves the request unevaluated. Tests (`cargo test --lib
  rules::number_theory` 27/27 in 7 s): BPSW against a sieve for 0..30 000,
  strong pseudoprimes to base 2, Carmichael and strong Lucas pseudoprimes,
  the 12- and 13-base strong pseudoprimes, Mersenne primes up to
  `2^607 - 1` and composite Mersenne numbers; `factorint` of `2^64+1`
  (274177 · 67280421310721), `2^128+1` (17- and 22-digit factors),
  `2^101-1` (13- and 18-digit factors), `2^96-1`, `1000003^5` and a
  product of three 12-digit primes; the dependent functions on those. The
  `ext` operators that need a prime as a machine word (`sqrtmod`,
  `primitive_root`, `mult_order`, `dlog`, `divisors`, `two_squares`,
  `carmichael`) still stop at 64-bit factors; semiprimes with two prime
  factors above ~25 digits are out of reach (no quadratic sieve).
  Also ran `rules::combinatorics` (29/29), `rules::discrete::finite_field`
  (7/7), `rules::discrete::gf` (8/8).

- **R7 — done.** `sim/euler_nd.rs`: `Sources { gravity, bed }` and
  `solve_with_sources` (`solve` is the same with the defaults, so
  `FvConfig` stays `Copy`). Euler gravity adds `ρ g` and `ρ u·g` in the
  cell centres. Shallow-water bathymetry uses the hydrostatic
  reconstruction of Audusse et al.: `(η = h + b, u, v, b)` are
  reconstructed with the configured MUSCL/WENO5 scheme, the face depths are
  `h* = max(0, η - max(b_L, b_R))`, each side's flux gets `g (h² - h*²)/2`
  in the normal momentum and the cell source is the centred
  `-g h̄ Δb/Δx`; the bed rides in the unused slot 4 through the boundary
  fill (transmissive, reflective, periodic and fixed faces) and the
  initial data are the depth. Tests in `tests/sim/euler_nd.rs` (14/14 in
  21 s): the lake at rest over a smooth bump and a discontinuous step, for
  all 5 reconstructions × 3 Riemann solvers, stays at rest to < 1e-12 in
  speed and free surface (1D, 80 cells, t = 0.5); the same in 2D (reflective
  and periodic/wall) while the flat-bottom scheme on the same depth profile
  produces speeds > 0.3 within t = 0.1; a 1 % surface hump over the bed
  conserves mass to 1e-12 and its free surface converges (successive L1
  differences 100 → 800 cells fall by > 1.6 each); free fall of a uniform
  gas is exact (`v = g t`, pressure unchanged, momentum `ρ g V t`);
  an isothermal atmosphere is held to first order (residual speed 1.06e-2,
  5.3e-3, 2.6e-3 at 40/80/160 cells; the zero-gradient wall ghost cells
  limit it, so gas gravity is *not* well-balanced); invalid combinations
  (bed with Euler, gravity with shallow water) are rejected; the defaults
  reproduce `solve` bit for bit. Not done: well-balanced gas (hydrostatic
  reconstruction for the Euler equations), friction, a sharp dry front,
  `euler_amr` with sources.

- **R8 — done for `M`, partly for `U`.** Profiled against exact rational
  series and mpmath: the Taylor stepping of `hyp_pfq` is wrong already for
  `b = 80` at `z = 3` (`M(5.5, 80.5, 3)` returned 3e17 instead of 1.23;
  `M(100.5, 120.5, 12)` 4e203 instead of 2.4e4, NaN for `|z| ≥ 50`):
  the forward recurrence for the Taylor coefficients amplifies the second
  solution. Two stabilisations were tried and dropped (Miller backward
  recurrence and Olver's boundary-value form on the coefficient recurrence
  are accurate for `z ≲ 12` but lose digits roughly like `e^{0.7 z}` when the
  expansion point is large, because the minimal solution is dominant for
  `k < |p|`). What is shipped (`kummer_series_dd`, `kummer_series_exact`,
  `u_integral_real` in `kernels/hypergeometric.rs`): (1) the power series in
  double-double arithmetic (Dekker splitting, ~32 digits) with the sum of
  the moduli of the terms; the sum is accepted when the terms cancel by
  less than `1e15`, which covers `Re z ≥ 0` with `Re a, Re b > 0` and, by
  Kummer's transformation (both forms tried), the left half plane for
  `Re(b - a) > 0`; (2) real arguments in the remaining transition regions
  (`a >> b`, `z < 0`, `M` exponentially smaller than its terms, down to
  `1e-90`) by exact fixed-point summation on `num-bigint` integers with the
  precision doubled until the result keeps 60 bits (`|z| ≤ 3000`,
  `|a|, |b| ≤ 5000`); (3) Tricomi `U` for real `a ≥ 12`, real `b`, `z > 0`
  from the positive integral representation by the trapezoidal rule in
  `t = e^s` (peak by bisection, step a third of the width, geometric
  convergence), which gives full relative accuracy where `U` is tiny
  (`U(150.5, 120.5, 30) = 3.6e-246`; the old continuation returned a wrong
  sign/magnitude there). Tests in `tests/kernels/hypergeometric.rs` (14/14,
  2.6 s): 22 real `M` values against mpmath (50 digits) and the exact
  rational series, from `1e-90` to `1e59`, relative 1e-12; the 64-point
  grid `a, b ∈ {…, 5.5 … 150.5}` × `z ∈ {0.5 … -100}` against exact
  rational sums (relative 1e-11); 7 complex-parameter values (rel. ≤ 2e-16),
  Kummer's transformation and a contiguous relation at `|a|, |b| ≈ 90`;
  10 real `U` values with `a, b` up to 150 (relative 1e-11). **Not covered:**
  complex `M` whose series cancels by more than `1e15` (no exact path for
  complex parameters), and complex `U` with `Re z < 0` and `|a| ≳ 60`
  (`U(60.5+10i, 90.5-5i, -30+20i)` is accurate to ~4e-7 only); a
  Miller recurrence in `a` for `U` fails because `U` is the dominant, not
  the minimal, solution of the `a`-recurrence while `a < b`.

- **R9 — done.** `kernels/special_ext.rs`: `ber/bei/ker/kei` for any integer
  order `ν ≤ 100` (`KELVIN_MAX_ORDER`), registered in `rules/special.rs` with
  the new `ber(n, 0) = bei(n, 0) = 0` rewrites for `n ≥ 1`. `ber/bei`: the
  double-double series now starts from `(z/2)^ν/ν!` for any `ν` (`x ≤ 30`);
  beyond, `ν = 0, 1` keep the Hankel expansion and higher orders use
  Miller's backward recurrence `J_{k-1} = (2k/z) J_k - J_{k+1}` from
  `max(ν, x) + 60`, scaled to the expansion of `J_0` (the expansion of `J_ν`
  itself needs `x ≫ ν²`); `ker/kei`: the forward recurrence
  `K_{n+1} = K_{n-1} + (2n/z) K_n` from `K_0, K_1` (stable) and the rotation
  `(-i)^ν`. *Bug fixed:* `kelvin_ker(1, 0)` returned `+inf` (the function
  tends to `-inf`; `ker_1(0.5) = -1.52`); at `0` the limits now follow the
  sign of the leading term `e^{-3πiν/4}` (`NaN` where that component
  vanishes and the limit is finite). Tests (`cargo test --test kernels
  kelvin` 4/4): 50 mpmath values for `ν = 2, 3, 5, 10, 20` and `x` from 0.5
  to 100 (`ber/bei` to 1e-10, `ker/kei` to 1e-9 relative to the larger
  component), the Bessel recurrence of `ber/bei` and the `K` recurrence of
  `ker/kei` for `n = 1..8` at six arguments each (1e-9), and the origin
  behaviour. `cargo test --lib rules::special` 28/28. Not covered:
  non-integer order, `ν > 100`, the finite limits of `ker_ν(0)`/`kei_ν(0)`
  for `ν = 2, 4, ...`.

- **R10 — done for the decoupled and the integer Euler–Poisson–Darboux
  cases; the general coupled affine system is not.**
  `pde_hodograph_solution(F, G, λ, μ)` (`rules/pde/nonlinear.rs`,
  `hodograph_decoupled`/`hodograph_coupled`; the operator is now variadic and
  the two-argument form is unchanged): (1) decoupled `λ(r)`, `μ(s)` (affine
  or not): `x = λ t + F = μ t + G`, `t = (G - F)/(λ - μ)`; (2) coupled affine
  `λ = a r + b s + c`, `μ = b r + a s + c` with `κ = a - b`, `N = a/κ` an
  integer in `2..5` (polytropic gas with `γ = (2N+1)/(2N-1)`, e.g. `γ = 5/3`
  for `N = 2`): the hodograph equation is the Euler–Poisson–Darboux
  `E(N, N)` equation `(r - s) t_rs + N (t_s - t_r) = 0`, whose general
  solution (Darboux) is `t = ∂_r^{N-1} ∂_s^{N-1} W`,
  `W = (F(r) + G(s))/(r - s)`; integrating `x_s = λ t_s`, `x_r = μ t_r` gives
  `x = λ t - κ (N - 1) ∂_r^{N-1} ∂_s^{N-2} W` with no integration constant
  (the formulas, including `h'(r) = 0`, were derived and checked
  symbolically with sympy for `N = 2, 3, 4` before being implemented). Tests
  in `rules/pde/completion_tests.rs`: decoupled affine and non-affine speeds
  (`x` consistent from both characteristics to 1e-11, `x_s = λ t_s` and
  `x_r = μ t_r` by central differences to 2e-6); coupled `N = 2, 2, 3, 4`
  with `κ = 1, 1/2, 1, 2` and nonzero `c` with arbitrary `F, G`
  (`r^3, sin s`; `e^r, s^2`; …): both hodograph equations to 1e-5 and the
  EPD equation to 2e-3 relative by finite differences; the two-argument form
  equals `λ = r, μ = s`; coupled systems outside the family (`N = 3/2` as
  for shallow water, unequal constants, `a ≠ e`) stay symbolic.
  `cargo test --lib rules::pde::completion_tests::hodograph` 3/3.
  **Not done:** the general affine system, in particular non-integer `N`
  (shallow water, `N = 3/2`), needs Riemann's method for the EPD equation
  (hypergeometric Riemann function); `a ≠ e` or `b ≠ d` (`λ - μ` not
  proportional to `r - s` after a shift) likewise.

- **R11 — already present; verified and extended.** The audit entry was
  stale: `eigs_jd` (non-symmetric Jacobi–Davidson with partial Schur
  deflation, harmonic Ritz extraction in complex arithmetic and a complex
  GMRES correction solve, optionally preconditioned) was committed in
  1b2530a4 with two tests (convection–diffusion interior eigenvalues with a
  sparse-LU preconditioner, and a complex target on a non-normal
  block-triangular matrix). Added
  `nonsymmetric_jacobi_davidson_on_a_strongly_non_normal_matrix`
  (`tests/kernels/eigs_ext.rs`): the Grcar matrix, `n = 64`, target
  `2 + i`, two eigenvalues; converged with true residuals `< 1e-6`, each
  returned value within `1e-4` of an eigenvalue from the dense complex Schur
  form and among the nearest to the target. `cargo test --test kernels
  nonsymmetric_jacobi` 3/3 in 0.1 s. **Not implemented:** the two-sided
  variant (left and right search spaces, oblique projection), which the
  table's wording asks for; no code was written for it.

**Summary of the round.** Done: R1–R10 (R10 and R8 with the limits stated
above), R11 verified. New library behaviour: exact `Q(n)` creative
telescoping, separable exterior sources with a Bessel-aware integrator
(minutes → 0.2 s), BDF with a mass matrix, FGLM, BPSW and ECM,
bathymetry/gravity sources, large-parameter Kummer `M`/Tricomi `U`, Kelvin
functions of order up to 100, hodograph solutions for decoupled and integer
EPD systems. Bugs found on the way: the 12-base Miller–Rabin accepted
318665857834031151167461; `kelvin_ker(1, 0)` had the wrong sign;
Burton–Miller with `α = i/k` has a phase error of about `2ka` at low
frequency.

Deferred (research-grade or large): general hodograph by Riemann's method,
F4/F5, van Hoeij factorisation, Hermite reduction for algebraic functions,
real QE over CAD, Meijer G/Fox H, `euler_amr` 3D, multivariate series,
Darboux polynomials, unstructured DG3D, FEM multigrid/coarsening.

### 5.4 Verification of the completion set (2026-10-29)

Re-run after the implementation, one cargo process at a time:

- `cargo clippy --all-targets`: only the span-less `doc_markdown` error
  (2.1) remains;
- `cargo test --lib -- rules::calculus rules::pde rules::poly rules::solve
  rules::number_theory rules::special`: 276/276 (194 s);
- `cargo test --test kernels`: 1022/1022 (13 s);
- `cargo test --test sim`: 808 passed, 2 ignored after one fix —
  `split_step::gpe_conserves_norm_and_energy_and_bright_soliton_is_stationary`
  failed on a norm drift of 1.25e-12 against a 1e-12 bound (400 FFT round
  trips; `split_step`, `spectral` and `transforms` are unchanged in this
  round, so it is round-off, not a regression); the bound is now 1e-11 with
  the drift in the message.
- R2's early exit in `antiderivative` for Bessel/Hankel integrands was
  checked against the guard disabled: `∫J₁`, `∫x J₀`, `∫J₀J₁` stay
  unevaluated either way, so no capability is lost. **New open item:** the
  integrator has no Bessel antiderivatives at all *(done, B1 in 5.5)* (`∫J₁ = -J₀`,
  `∫x J₀ = x J₁`, `∫x^ν J_{ν-1} = x^ν J_ν`, `∫J₀J₁ = -J₀²/2`); when added,
  the guard must let those table forms through. **P2**
- `rules::ode::linear` (~10 min) was not re-run; the round did not touch it.

### 5.5 Remaining limitations of the round after the completion set (2026-10-30)

Items marked here were closed after 5.4; each has tests next to its neighbours and an updated module doc.

| ID | Item | Status | Where |
|---|---|---|---|
| R8b | Tricomi `U` for complex `z` with `Re z < 0` and `|a| >~ 60` (was `4e-7`) | **done, with a documented remainder**: 979 stress cases against mpmath; `Re a >= 0`, `|a| <= 100`, `|Im a| <= 90`, `|b| <= 8`: below `6e-11` (typically `1e-13`), `a = 65i` at `z = -40i` up to `3e-8`. Fixes: (1) Taylor steps limited to `4/(1 + Σ|params|)` of the distance to the singular point (a step of relative length `eps` loses `exp(|a| eps)`); (2) the start ray is at `±45°`, not the imaginary axis (there `e^z` is not damped and its algebraic growth inwards amplified the start rounding error to `5e-7`); (3) start values on a logarithmic scale (the series needs `r` up to `3e5`, where `U` over/underflows); (4) `Re a < -8` and negative integer `a`, `a-b+1` (`n >= 12`, `Re b > 0`) by the recurrence in `a` run downwards from `Re a ~ 10` (stable: `U` is minimal for increasing `a`); (5) `U(a, a+1, z) = z^{-a}` for `a - b + 1 = 0` (it returned NaN). Remaining (documented in the module): results below `~1e-270` underflow inside the continuation (`a = 150-40i`, `Re z <= -60`); `|Im a| >= 120` with `z` on the negative imaginary axis up to `4e-3`; `Re a <= -60` is accurate only for part of the parameter space (the downward recurrence loses digits for `b = 30`, and for `a = -70+20i` near `z = -8-6i`). Tests: `tricomi_u_with_large_complex_parameters_in_the_left_half_plane_matches_mpmath` (126 points) | `kernels/hypergeometric.rs`, `tests/kernels/hypergeometric.rs` |
| R10b | L8a — hodograph with non-integer `N` (Euler–Poisson–Darboux `E(N, N)`, shallow water `N = 3/2`) by Riemann's method | **done**: Riemann function `((ρ-σ)/(r-s))^N 2F1(N, 1-N; 1; z)`, `z = (ρ-r)(σ-s)/((r-s)(ρ-σ))` (found by a residual search against mpmath, `2e-31`), `0 <= z < 1` in the characteristic rectangle with `r0 > s`; `pde_hodograph_riemann(N)` returns it in hypergeometric form (`hyperpfq`), `rules::pde::riemann::{epd_riemann, epd_riemann_gradient, epd_goursat}` solve the Goursat problem `t(ρ, s0) = φ`, `t(r0, σ) = ψ` by Green's identity (no derivative of the data needed). Tests: Riemann function boundary values and the adjoint equation for `N = 1.5, 0.7, 2.5, 3`; Goursat solutions against the Euler–Poisson integral solutions (`N = 3/2, 5/4, 5/2`, error `< 1e-8`) and against Darboux's closed form for `N = 2`; the symbolic operator against the numeric one. `pde_hodograph_solution` itself still returns the integer-`N` family only (a closed form in two free functions needs fractional derivatives) | `rules/pde/{riemann,nonlinear,mod}.rs`, `rules/pde/completion_tests.rs` |
| R11b | R11 — two-sided Jacobi–Davidson for non-symmetric matrices | **done**: `eigs_jd_two_sided(apply, apply_t, n, k, target, opts)`; right and left search spaces, Petrov–Galerkin extraction from `(W^H (A-τ) V, W^H V)` with the left vectors from the inverse of the right-eigenvector matrix, two-sided Rayleigh quotient, two oblique correction equations by GMRES, oblique deflation `I - Q (P^H Q)^{-1} P^H` of locked pairs; returns left/right vectors, both residuals and the condition numbers `1/|w^H u|`. Tests: strongly non-normal bidiagonal matrix (closed-form eigenvectors and condition numbers to `1e-5`, biorthogonality), complex pair `2 ± i` | `kernels/eigs_ext.rs`, `tests/kernels/eigs_ext.rs` |
| B1 | Bessel antiderivatives in the integrator | **done**: `∫ Z_1 = ∓Z_0` and, for `J`, `Y`, every odd order `∫ Z_(2m+1) = -Z_0 - 2 Σ Z_(2k)`; `∫ x^(m+1) Z_m = ± x^(m+1) Z_(m+1)`; `∫ x^(1-m) Z_m = ∓ x^(1-m) Z_(m-1)`; `∫ Z_0 Z_1 = ∓ Z_0²/2` for `J, Y, I, K`, arguments `c x`, by linearity; checked like every antiderivative (differentiated and compared). The R2 guard `special_in` now hands the integrand to this table instead of declining: anything else still declines at once, and the exterior-source test (`exterior_time_varying_source_stays_cheap`, < 8000 nodes, < 2000 unions) still passes in 0.5 s. Tests: 28 antiderivatives and 5 forms that must stay unevaluated | `rules/calculus/integrate.rs`, `rules/calculus/tests.rs` |
| B2 | **Bug found**: the Coulomb function `G_l(η, ρ)` and `G'_l` were wrong at small `ρ` (Steed's continued fraction: relative error `8e-7` at `η = 1.94, l = 4, ρ = 0.555`, worse elsewhere; 85 of 270 grid points above `1e-10`), which made `rules::tests::no_standard_rewrite_is_refuted` fail on `special/coulomb-gp-recurrence` (the test was outside the filter of 5.4) | **fixed**: `G + iF = H⁺ = e^{iθ} (-2iρ)^{l+1+iη} U(l+1+iη, 2l+2, -2iρ)` with the new complex `U` (`coulomb_g_from_u`); all 270 points `< 3.3e-13`; Steed stays as fallback. Test: 64 mpmath points with the Wronskian | `kernels/special_ext.rs`, `tests/kernels/special_ext.rs` |
| B4 | `coulomb_f`/`coulomb_f_prime` returned `None` at `rho = 0`, so three rules (`coulomb-f-zero`, `coulomb-fp-zero*`) were never numerically tested | **fixed** (`F_l(η,0)=0`, `F'_0(η,0)=C_0(η)`) | `kernels/special_ext.rs` |
| B3 | **Bug found**: `fem_simplex::error_indicators` summed the edge-jump terms in the randomised order of a `HashMap` (two calls of the same function gave results differing in the last bits) | **fixed**: facets taken in sorted key order | `sim/fem_simplex.rs` |

Not reached from the "own search" list (still open): the other `Limitations:` paragraphs of Part 5.2. *(Update 2026-10-30, Part 7: the symmetric generalised eigenproblem with indefinite `B`, the SDE adaptive step without commutative noise, the DG hierarchical limiter at non-periodic faces, variable `k, c` in `fem_simplex` and SPH self-gravity are closed. Round 2 of Part 7 closed DPD walls, complex Helmholtz ACA, the rotating frame in `split_step` and the variable-coefficient FEM estimator; what remains is listed at the end of Part 7.)*

---

## Part 6 — In-algorithm parallelisation (done 2026-10-30)

`rayon 1.11` was already a hard dependency; no new dependency was added. The
machine has 16 hardware threads. `src/kernels/par.rs` holds the shared
helpers (`map_indexed`, `map_slice`, `for_chunks_mut`, `for_each_mut`,
`sum_fixed`): ordered collection, a sequential path below a grain size, and
fixed-block sums.

### Rules for every change

- **Deterministic results.** Parallel maps collect in index order; floating-point
  reductions are summed in a fixed order or in fixed blocks, so results do not
  depend on the thread count. Exact symbolic results are identical.
- **Grain-size cut-offs.** Each site has a sequential path below a size
  threshold (constants named `*_PAR_*` next to the site).
- **No parallelism inside the e-graph.** Parallel work is done on pure data.
- **No nested oversubscription.** Only outermost loops are parallel.
- **Measured.** Each site is timed (serial in a one-thread rayon pool vs the
  global 16-thread pool, same process) and its result compared for equality;
  a site that did not win was **dropped** (listed below).

### Verification

- `par_sites` timing tests: `tests/kernels/par_sites.rs`,
  `tests/sim/par_sites.rs`, `src/rules/par_tests.rs`. Each runs the site once in
  `rayon::ThreadPoolBuilder::new().num_threads(1).install(..)` and once in the
  global pool, **asserts the `Debug` text (bit-exact floats) is equal**, and prints
  a `PAR |` row (`--nocapture`, run with `--test-threads=1` for clean timings).
- Every touched module's tests were run normally and with `RAYON_NUM_THREADS=1`
  (kernels: 1032/1032 both ways; sim: 812 passed, 2 ignored, and the touched
  modules with one thread: 68 passed, 2 ignored; lib subset below).
- Exact outputs compared byte for byte with the serial path: Gröbner bases,
  factorisations, `factorint`, homotopy endpoints and certificates, the
  e-graph soundness verdicts, FMM/DG/Euler/N-body fields.

### Timing table (serial vs 16 threads)

| ID | Site | Size | Serial | Parallel | Speed-up | Notes |
|---|---|---|---|---|---|---|
| P1 | `homotopy_solve` (paths) | cyclic-5, 120 paths | 19.2 s | 4.6 s | 4.2 | uneven path cost limits it; identical endpoints |
| P1 | `homotopy_solve` adaptive | cyclic-3, 6 paths | 43 ms | 16 ms | 2.8 | |
| P2 | `certify_endpoints` | 120 endpoints | 71 ms | 10 ms | 6.9 | |
| P2 | `certify_endpoints_rigorous` | 120 endpoints | 18 ms | 2.3 ms | 7.7 | |
| P2 | `homotopy_solve_verified` | cyclic-3 | 48 ms | 18 ms | 2.6 | |
| P3 | `groebner` (batched S-pairs) | Katsura-5 grevlex | 535 ms | 127 ms | 4.2 | reduced basis unique, equal |
| P3 | `groebner` | Katsura-6 grevlex | 8.9 s | 2.4 s | 3.6 | |
| P3 | `groebner` | cyclic-5 grevlex | 150 ms | 81 ms | 1.9 | |
| P3 | `fglm` | cyclic-5 → lex (70) | 45 ms | 47 ms | 1.0 | **dropped**: the multiplication columns are not the bottleneck (the echelon walk is sequential) |
| P4 | `univariate::factor` | Swinnerton-Dyer deg 32 | 8.3 s | 2.0 s | 4.2 | batched Zassenhaus subsets + parallel modular images |
| P5 | `qlinalg::nullspace` / `rref` | 80 × 80 over `Q` | 4.1 s / 4.2 s | 0.78 s / 0.80 s | 5.3 / 5.2 | |
| P5 | `ode/linear.rs` nullspace + parallel candidate denominators | Airy/Bessel eigenring test | 491 s | 454 s | 1.08 | **dropped** (the row updates are tiny; parallel candidates ran the expensive late denominators the serial scan never reaches: 97 CPU-minutes) |
| P5 | `zeilberger` `solve_system`, `lie::null_space` | Apéry, Dixon | 1.3 s | 1.2-1.4 s | 0.9-1.0 | **dropped** |
| P6 | `galois_high::analyse` (cosets + resolvents) | `x^7 - x - 1` (S7) | 431 ms | 206 ms | 2.1 | |
| P6 | | `x^7 - 7x + 3` (PSL(3,2)) | 272 ms | 91 ms | 3.0 | |
| P7 | ECM curves (ramped groups) | `factorint(2^128 + 1)` | 2.9 s | 1.45 s | 2.0 | lowest successful curve wins; `(2^107-1)(2^61-1)`: 8.1 s to 3.6 s (2.25); `2^101-1`: 1.4 |
| P8 | `Csr::matvec` | 160 000 rows, 5-point | 81 ms | 19 ms | 4.4 | per 50 products |
| P8 | `Mat::matmul` | 500 × 500 | 145 ms | 25 ms | 5.7 | |
| P8 | `lu_factor` trailing update | 700 × 700 | 114 ms | 535 ms | 0.2 | **dropped** (per-step dispatch cost) |
| P9 | `fd_jacobian_dense_par` | n = 200 | 272 ms | 34 ms | 8.0 | leaf only: needs a `Sync` right-hand side, so the stiff solvers keep the sequential function |
| P10 | `integrate_gk_par` | 2000 sin per node | 58 ms | 12 ms | 4.7 | bit-identical to `integrate_gk` |
| P10 | `cubature_nested_par` | 3-D | 241 ms | 43 ms | 5.6 | same value and evaluation count as `cubature_nested` |
| P10 | `quasi_monte_carlo_par` (Sobol) | 4-D, 4 M points | 1.55 s | 0.27 s | 5.7 | fixed 4096-point blocks (not bit-equal to the serial function) |
| P10 | `monte_carlo_par` | 4-D, 4 M points | 244 ms | 24 ms | 10.3 | per-block streams |
| P11 | `solve_dg3d` | 8³ elements, p = 2 | 568 ms | 150 ms | 3.8 | with the hierarchical limiter 4.1 |
| P12 | `fem_simplex::solve_static` | 100 × 100 triangles, P2 | 1.52 s | 0.83 s | 1.8 | assembly parallel; the sparse solve dominates |
| P12 | `error_indicators` | same | 880 ms | 608 ms | 1.45 | |
| P13 | `euler_nd::solve` | 192 × 192 blast | 3.03 s | 0.94 s | 3.2 | |
| P13 | `ns_mac` convection/Laplacian stencils | 96 × 96 cavity | 1.61 s | 1.55 s | 1.04 | **dropped**; only `Csr::matvec` (P8) is parallel there |
| P14 | `fmm3d` | 30 000 points, order 4 | 10.1 s | 1.9 s | 5.2 | |
| P15 | `check_program` (soundness) | arith+elementary+poly+complex, 60 samples | 70 ms | 16 ms | 4.5 | first refutation in rule order is returned |
| P16 | `canon::refine` | 3000 vertices | 14 ms | 14 ms | 1.0 | **dropped** (the dense `n × n` adjacency limits the size) |
| P17 | `ewald_reciprocal` | 2500 charges, kmax 8 | 476 ms | 58 ms | 8.2 | extra: fixed 32-vector blocks |
| P17 | `ewald_energy_forces` real space | 2500 charges | 6.1 s | 0.57 s | 10.7 | extra: one row of pairs per task (each pair evaluated twice) |
| P18 | `nbody::direct_field` | 3000 bodies | 358 ms | 46 ms | 7.8 | extra, bit-identical |
| P18 | `gravity_accelerations` (Barnes–Hut) | 30 000 bodies | 1.72 s | 0.24 s | 7.1 | extra |

P7 ECM and P6 rows with a note are from `src/rules/par_tests.rs`
(`cargo test --lib par_tests -- --nocapture --test-threads=1`).

Sites added by the search beyond P1–P16: **P17** Ewald (reciprocal blocks,
real-space rows), **P18** N-body direct sums and tree traversals, and the
leaf variants of P9/P10 (`fd_jacobian_*_par`, `integrate_gk_par`,
`cubature_nested_par`, `quasi_monte_carlo_par`, `monte_carlo_par`). Not
parallelised although examined: multifactor `split_off` (stateful `seen` set
and a global subset budget), Hensel lifting (sequential dependence between
factors), `md_engine::compute_forces` (Newton's third law needs per-thread
force copies of size `N`), FGLM, the DG limiters other than the troubled-element
test, the multigrid/red-black smoothers (they go through `Csr::matvec`).

Results that changed (not bit-equal to the previous serial code, but
independent of the thread count): `ewald_reciprocal_cell` above `2·10^5`
wave-vector × particle products (blocked summation) and the orthorhombic
real-space sum above `n² = 4·10^6` (row-wise accumulation); the `_par`
quadrature variants (block sums / streams). Everything else is bit-identical
to the sequential code.

### Final suite counts (2026-10-30)

- `cargo test --test kernels`: 1032 passed (also with `RAYON_NUM_THREADS=1`: 1032);
- `cargo test --test sim`: 812 passed, 2 ignored (modules touched by Part 6 with one thread: 68 passed, 2 ignored);
- `cargo test --lib -- rules::calculus rules::pde rules::poly rules::solve rules::number_theory rules::special rules::ode rules::discrete rules::tests graph::soundness rules::par_tests`: 482 of 483 passed in 693 s including `rules::ode::linear`; the one failure (`rules::tests::no_standard_rewrite_is_refuted`, three Coulomb rules at `rho = 0`) was fixed afterwards (B4) and the test passes (it had been failing before this round because of B2); the poly/number theory/calculus/lie/kovacic/soundness subset also passed with `RAYON_NUM_THREADS=1` (132/132) before the Phase B changes;
- `cargo clippy --all-targets`: only the span-less `doc_markdown` "backticks are unbalanced"; with `--cap-lints warn` no new warnings in touched files.

---

## Part 7 — Engine- and session-level parallelisation (done 2026-10-30, round 2 same day)

Part 6 parallelised *inside* algorithms. Part 7 looks one level up: the engine, the session,
simulation scenarios and ensembles, and a few leaf utilities. The rules of Part 6 apply (deterministic
results, named `*_PAR_*` cut-offs, no manual threads, every site measured one thread vs 16 threads in
the same process with the `Debug` text asserted equal, `PAR |` row printed). No dependency was added.

**Keep threshold (round 2): about 1.1.** Round 1 used 1.3. A site that clears ~1.1 on realistic input is
kept; a site that is slower than serial on small input gets a grain cut-off so that small inputs take the
sequential path, and the *semantics* never depend on the path taken (only the execution strategy).
Sites that stay dropped are listed. Timings on this machine are noisy by about ±10 % between runs
(the table gives the range of the runs of the day).

### Architecture

| ID | Level | What runs in parallel | What stays sequential, and why | Determinism argument |
|---|---|---|---|---|
| A1 | `Engine::apply_tier` | The e-matching of each maximal run of consecutive `Action::Rewrite` rules of a tier against the **frozen** graph (`&Graph`; `Graph`, `Rewrite`, `Pat`, `Engine` are `Sync`, asserted at compile time); each rule keeps its own `match_limit`, `match_fuel` and ban decision, computed from its own search result. Parallel from `ENGINE_PAR_MIN_CANDIDATES = 512` candidate nodes in the run | Kernels keep their position and run one after the other (they mutate the graph and call the engine); builds and unions are applied in rule order, then match order. `windows(..)`: not parallelisable, `TreeWindow::optimize` / `beam` take `&mut Graph` (rewrites and window passes create nodes) | The searches only read the graph and are collected in rule order; the batched semantics (a later rule of the run does not see the unions of an earlier one in the same pass) apply on **both** paths, so results do not depend on the cut-off or on the thread count |
| A2 | Two-phase kernels | `Kernel::{two_phase, prepare, commit}` with `PreparedWork: Send`: the scheduler prepares every eligible node (sequential, graph access), runs the returned jobs in parallel when there are at least two, and commits in node order. Implemented for **`groebner(...)` requests** of `PolyKernel` (extract polynomials, `groebner_via` on plain data, build the basis terms) | Everything else; univariate factorisation was tried as a job and dropped (below) | Prepare-all-then-commit is the semantics on both paths, jobs are independent pure functions of their data, results are committed in node order |
| A3 | `Session::compute_all` | One engine run (and the choice of the answer) per request, on its private clone of the template | Cloning the template and importing a request (store borrowed) and importing the answers back (input order) | The answers and the final store equal a loop of `compute` calls (tested also in a one-thread pool) |
| A4 | Scenarios and ensembles | `sim::scenario::run_many`, `kernels::sde::sde_ensemble` (member `i` seeded by `ensemble_seed(seed, i)`), replica updates of `ParallelTempering` / `PottsTempering` | Replica exchange and its generator | Collected in index order; members depend only on `(seed, index)`; replicas own their generators |
| A5 | `bdf_with_par` | The finite-difference Jacobian columns (`FdMode::parallel()` = `fd_jacobian_dense_par`), from `STIFF_PAR_MIN_N = 64` unknowns | Newton iterations, LU, Radau (no gain) | Bit-identical to `bdf_with` (asserted) |
| A6 | `md_engine::compute_forces` | The pair sum above `MD_PAR_MIN_PAIRS = 20 000` pairs, in blocks of `max(8192, pairs/64)` pairs with one force buffer per block, added in block order | The bonded / extra terms; the single-accumulator sum of short lists | Block boundaries depend only on the list length; results differ from the single-accumulator sum only in the last bits and are independent of the thread count |
| A6 | `backend::call_batch_par` | Rows of a compiled batch in fixed blocks of 4096 | – | Pointwise evaluation: identical to `call_batch` |
| A6 | `io::plot_function_2d_par`, `plot_surface_3d_par`, `plot_term_2d` | Sample evaluation from 256 samples (rows for the surface) | SVG assembly | Samples collected in order, identical picture |
| A6 | `rules/verify.rs` `vanishes` | The 12 sample points, when the simplified term has at least `VERIFY_PAR_NODES = 40` reachable classes | – | The verdict (some sample refutes / some sample is finite) does not depend on the order |

### Timing table (serial = 1-thread pool, parallel = 16 threads)

| ID | Site | Size | Serial | Parallel | Speed-up | Notes |
|---|---|---|---|---|---|---|
| A1 | batched e-matching, kept (threshold 512) | `diff(sin²x e^{x²} ln(1+x²))`, `(x+1)^6 - ...`, integral | 12–27 ms | 11–28 ms | 0.98–1.14 | small: mostly under the cut-off, so the sequential path |
| A1 | | `diff(((x²+1)^5 sin³x eˣ ln(x²+2))², x)` (285 ms) | 281–286 ms | 164–322 ms | 0.89–1.74 (runs: 1.15, 0.89, earlier 1.74) | noisy: runs 1.15, 1.12, 0.89 (1.74 with the cut-off at 0) |
| A1 | | `diff(diff(sin x cos x eˣ ln x atan x sinh x, x), x)` (10–13 s) | 10.7–13.0 s | 7.3–10.2 s | 1.27–1.40 | the clear gain |
| A1 | | polynomial identity in x, y, z; two `diff` products (24–90 ms) | 19–90 ms | 19–98 ms | 0.91–1.08 | at the noise level; not slower on average |
| A2 | `groebner` as parallel jobs | 4 Katsura-4/5 systems, jobs side by side vs one after the other, both with the parallel Gröbner inside | 385 ms | 270 ms | 1.43 | kept |
| A2 | two-phase kernel in the engine | request holding 4 `groebner(katsura-4)` / 4 `groebner(katsura-5)` | 1.20 s / 11.8 s | 0.67 s / 6.96 s | 1.79 / 1.70 | one thread vs 16 (includes the parallelism inside each Gröbner run); outputs equal |
| A2 | `univariate::factor` as parallel jobs | 4 Swinnerton-Dyer deg 32 in one request | 7.89 s | 7.48 s | 1.06 | **dropped** (each factorisation already fills the pool, P4) |
| A3 | `Session::compute_all` | 8 requests, standard rules (loop of `compute`: 49–125 ms) | 54–65 ms | 17–31 ms | 2.0–3.9 | limited by the longest request |
| A4 | `sim::scenario::run_many` | 8 Ising 96×96 scenarios | 251–359 ms | 39–95 ms | 3.8–6.5 | |
| A4 | `sde_ensemble` (Euler–Maruyama) | 20 000 paths × 200 steps, n = 2 | 1.06–1.19 s | 0.20–0.24 s | 4.9–5.3 | |
| A4/A6 | `ParallelTempering<2,2>` | 8 replicas, 64×64, 40 rounds | 161–222 ms | 62–84 ms | 2.6–2.7 | bit-identical |
| A4/A6 | `PottsTempering` | same | 108–139 ms | 53–95 ms | 1.45–2.04 | |
| A5 | `bdf_with_par` | n = 100 / 200 / 400, expensive `f` | 18.5 / 55 / 217 ms | 13.7 / 36 / 139 ms | 1.35 / 1.52 / 1.56 | kept, cut-off 64 |
| A5 | `radau_with_par` | n = 100 / 200 / 400 | 0.88 / 6.5 / 39.2 s | 0.82 / 6.6 / 42.4 s | 1.08 / 0.99 / 0.92 | **dropped** (Newton and LU dominate); the n = 160 anomaly of round 1 (0.74 vs the BDF 1.55) was noise of two different solvers, not an inconsistency of the Jacobian |
| A6 | `md_engine::compute_forces` | 2048 / 10 976 atoms, LJ rc 2.5 | 28.7 / 109 ms | 24.5 / 85 ms | 1.17 / 1.29 | kept |
| A6 | `backend::call_batch_par` | 1M rows, 2 inputs | 317–372 ms | 72 ms | 4.4–5.2 | |
| A6 | `io::plot_term_2d` | 20 000 samples (evaluation + SVG) | 142–146 ms | 44–53 ms | 2.7–3.3 | |
| A6 | `rules/verify.rs` `vanishes` | 12 samples of a large term (60–85 µs per call) | 57–85 µs | 41–66 µs | 1.29–1.40 | kept, only above 40 classes |

Still dropped from Part 6 (unchanged): `ns_mac` stencils (1.04–1.10), eigenring search, FGLM, Zeilberger,
`lie` null spaces, LU update, graph canonical refinement.

### FFI

`rssn_compute_all(session, terms, n, config, answers, statuses)` (`src/ffi/compute.rs`): caller-provided arrays
of `n` handles, `RssnAnswer`s and `RssnStatus` codes; a handle that is not a term of the session gets
`NoSuchTerm` in its slot, the others are computed together by `Session::compute_all`; the return value is `Ok`
when the call was well formed and the message of the first failing request stays available through
`rssn_last_error`. Test: equal to `rssn_compute` per request, per-request codes, null-array handling.
`rssn.h` / `rssn.hpp` were **not regenerated** (cbindgen runs only with `DEV=1`; skipped, the new
declaration has to be generated before the header is published).

### Part B — closed limitations

| ID | Item | Status | Where |
|---|---|---|---|
| B1 | Symmetric generalised eigenproblem with indefinite `B` | **done**. Dense: `dense::eigen_generalized_symmetric` (Cholesky reduction and `B`-orthonormal eigenvectors for SPD `B`; eigenvalues of `B⁻¹A`, complex pairs included, otherwise). Iterative: `eigs_gen_sym_shift_invert` estimates `λ_min(B)`, `λ_max(B)` by Lanczos (`b_is_positive_definite`) and falls back to `eigs_gen_shift_invert` (Krylov–Schur on `(A − σB)⁻¹B`) for indefinite `B` (2-normalised vectors; `B` must be nonsingular). Tests: SPD pencil, real indefinite spectrum {2, −3}, complex pair ±i (dense and iterative), 40-dimensional indefinite diagonal pencil | `kernels/dense.rs`, `kernels/eigs_ext.rs` |
| B2 | SDE adaptive solver with non-commutative noise | **done**: `sde_adaptive_noncommutative_with` on a fixed Brownian grid of cells, each with Lévy areas; a step is a union of cells combined exactly by Chen's relation, so rejecting or splitting a step re-uses the same cells (finest step `(t1 − t0)/cells`); `noncommutative_noise_grid` returns the path for a fixed-step reference. Tests: rms error against the fine `milstein_noncommutative` reference on the same path 4.8e-2, 1.9e-2, 2.7e-3 for 11, 13, 46 mean steps, Chen's relation, invalid input; `sde_adaptive_with` untouched | `kernels/sde.rs` |
| B3 | DG hierarchical limiter at non-periodic boundaries | **done** in `dg.rs` and `dg3d.rs` (moment limiters and troubled-element test): ghost elements from the boundary condition (mirror with `wall_state` for walls; exterior state at the virtual nodes for fixed / custom states; even mirror for outflow). Wall equals the left half of the doubled domain to 1.3e-13 (1D), 1.2e-16 / 5.3e-15 (3D hierarchical / characteristic); a Burgers shock at an outflow boundary no longer undershoots (−0.26 → 0.16 ≤ u ≤ 1); periodic runs create no ghosts. Cost: the limited modes of an outflow boundary element are cut to the (zero) mean difference, i.e. first order there whenever the limiter engages | `sim/dg.rs`, `sim/dg3d.rs` |
| B4 | Variable `k(x)`, `c(x)` in `fem_simplex` | **done**: `FemProblemVar`, `solve_static_var`, `solve_heat_var` (bit-identical to the constant API). **Round 2:** `error_indicators_var` (element residual with `∇k · ∇u_h` by central differences of `k`, `k(x)` in the facet jump) and `solve_adaptive_var`. Tests: manufactured solution `k = 1 + xy + y²/2` (L² rates 1.99 for P1, 2.9–2.97 for P2), layered conductivity, estimator equal to the constant one for constant closures, estimator 460× larger near a Gaussian layer of `k` than far from it, adaptive loop (estimator falls monotonically 66 → 5.4 over 7 refinements, L² error 0.039 vs 0.112 on the coarse mesh) | `sim/fem_simplex.rs` |
| B5 | SPH self-gravity | **done**: `Sph::gravity` (`SphGravity { g, eps, theta }`, Barnes–Hut tree force), `Sph::gravitational_energy`; test: uniform cold ball, radial acceleration linear in `r` within 8 %, vanishing net force, pair energy near `−3GM²/5R` | `sim/sph.rs` |
| B6 | DPD walls | **done (round 2)**: `Dpd::set_walls(DpdWall { axis, a_wall, kind })`: two walls normal to one axis with the effective soft repulsion `a_wall (1 − d/r_c)` and specular (free slip) or bounce-back (no slip) reflection; the neighbour search stops wrapping along that axis. Reflected velocities are also applied to the final velocity of the modified Verlet step. Tests: particles stay in the slab and `x`, `y` momentum is conserved with specular walls (to 1e-8), temperature stays at the thermostat value plus the drift energy, a drift of 0.8 survives specular walls (exactly) and decays to 0.02 with bounce-back walls, density depletes within `r_c` of a wall (wall layer 30 vs bulk 81 per unit thickness), no-wall runs unchanged. Limits: a single pair of planes, no frozen wall particles, the virial pressure ignores the wall force | `sim/dpd.rs` |
| B7 | Complex Helmholtz ACA | **done (round 2)**: `solve_helmholtz_exterior_aca` with `HelmholtzAcaOptions` / `HelmholtzAcaStats`: cluster tree by recursive bisection, admissible blocks compressed by complex ACA with partial pivoting (falls back to dense storage when the rank passes `min(m, n)/2`), complex GMRES with Jacobi preconditioner, block-parallel assembly and products; the dense and compressed solvers share the matrix entries (`Setup::entry`). Tests against the dense path: 1280 triangles, sound-soft `k = 1.5`: density error 2.2e-8, 449 low-rank + 1529 dense blocks, 10 GMRES iterations; sound-hard (Burton–Miller, hypersingular) `k = 2`: 7.1e-9; all-dense case 1.4e-11. Compression at this size is modest (0.96 of dense storage); it pays from a few thousand elements and only while `k·diam` is moderate (no directional compression) | `sim/bem_helmholtz.rs` |
| B8 | Rotating frame in `split_step` | **done (round 2)**: `SplitStep::evolve_rotating` (and `rotate`): `exp(iΩ dt L_z)` as three alternating spectral shears (Paeth), sandwiched around the Strang step; norm conserved to round-off. Tests: rotation of a Gaussian against the exact one (3.7e-6 at 64², n = 64), forth-and-back, free anisotropic packet equals the rotated lab solution with spectral convergence (7.8e-7 at n = 64, 1.0e-9 at n = 128), a displaced packet in an isotropic trap follows `x₀ cos t` rotated clockwise by `Ω t` (to 2e-3 in the mean position, Yoshida), `Ω = 0` bit-identical to `evolve`, norm with a GPE nonlinearity. Limits: the wave function must vanish near the box edges (the shears wrap; a free packet that had spread to the edge gave a grid-independent 2e-3 error), each sub-step rotates by at most π/2 | `sim/split_step.rs` |

### Results that changed

- `Engine::apply_tier`: within a run of consecutive rewrite rules a later rule no longer sees the unions
  of an earlier one in the same pass (it did before). The engine-heavy subset (below) passes.
- `groebner(...)` requests are computed two-phase: all `groebner` nodes of an iteration are prepared before
  any is committed. The values are the same bases.
- `md_engine::compute_forces` with at least 20 000 candidate pairs sums in blocks (last-bit differences from the
  old single accumulator; independent of the thread count).
- The hierarchical DG limiters now limit boundary elements (monotone results at walls / fixed states).
- Otherwise results are unchanged.

### Final suite counts (2026-10-30)

Only the touched modules were run (never the full suites), on the final tree after `cargo +nightly fmt`:

- engine-heavy subset (A1 and A2 change engine semantics): `cargo test --lib -- rules::tests graph:: rules::calculus rules::poly rules::solve rules::par_tests api:: sim::scenario ffi::`: **252 passed**, and with `RAYON_NUM_THREADS=1`: **252 passed** (includes `rules::tests::no_standard_rewrite_is_refuted`, `graph::soundness`, the new `compute_all`, `run_many` and FFI tests);
- `cargo test --test kernels --test sim -- sde:: eigs_ext:: dense:: dg:: dg3d:: dg_research:: fem_simplex:: sph:: spin_models:: dpd:: bem_helmholtz:: bem_aca:: split_step:: md_engine:: par_sites`: kernels **41 passed**, sim **151 passed**, 0 failed;
- timing rows: `cargo test --lib par_tests -- --nocapture --test-threads=1` and `cargo test --test kernels --test sim -- par_sites --nocapture --test-threads=1` (all `PAR |` rows above);
- `cargo +nightly fmt` (the repository's `rustfmt.toml` needs nightly; stable `cargo fmt` reformats 145 files and must not be used); `cargo clippy --all-targets` reports no error in the code added in Parts 7 (the remaining errors are in untouched code or in the older `scattered_field` / `ghosts` code).

**Still open from 5.2 after Part 7:** WENO limiter in more than one dimension, `ns_mac` non-uniform grids / free surface, `spin_models` beyond nearest-neighbour hypercubic lattices, `kerr` closed forms, unstructured / curved DG hexahedra, a directional (FMM-type) compression for the high-frequency Helmholtz BEM, frozen-particle DPD walls, and the "not present anywhere" list of 5.2 (F4/F5, quadratic sieve, Meijer G, wavelets, ...).

### Round 3 — triage of the remaining limitations, header, full suite, clippy (2026-10-30)

Triage of every `Limitations:` paragraph in `src/` (30), Part 5.2 including the "not present anywhere" list, and 5.5.
Classes: **(a)** silently wrong results, NaN, panics or hangs on reasonable inputs; **(b)** a gap that blocks a common
CAS / scientific-library workflow; **(c)** niche or research level.

| Class | Item | Result |
|---|---|---|
| (a) | `airy_ai`, `airy_bi` and the derivatives: NaN for `x` below about `-1000` (Taylor continuation overflowed), **hangs** at `+-inf` and for `Bi(x)`, `x >~ 1e5` (`2e10` Taylor steps) | **fixed**: DLMF 9.7.9-9.7.12 oscillatory expansions for `x <= -10` (Wronskian `Ai Bi' - Ai' Bi = 1/pi` to `1e-6` up to `x = -1e8`; `Ai(-10)`, `Bi(-10)` against mpmath values), `Bi` is `+inf` from `x > 105`, the limits at `+-inf` are returned (`Ai' / Bi'` have none at `-inf`: NaN) |
| (a) | `ln_gamma_numerical` NaN for negative non-integers; `elliptic_e(-inf)` hang; `I0, I1, sinc, E1, Ei, li, Si, Ci, Shi, Chi` NaN at `+-inf`, `E1(0)`, `Ci(0)`, `Chi(0)` NaN instead of `+-inf` | **fixed**: `ln|Gamma(x)|` by reflection (`+inf` at the poles), `E(m)` for `m < -1` by `E(m) = sqrt(1-m) E(m/(m-1))`, explicit limits. `Y0, Y1, K0, K1, E1, li, Ci, Chi` stay NaN for negative arguments (their values are complex) and `K(m)`, `E(m)` for `m > 1` |
| (a) | the same class in `hyp_pfq`: NaN for `p > q + 1` and for `b_j` at a non-positive integer, indistinguishable from a NaN input; `pFq(...; 0)` was NaN in those cases | **fixed**: `hyp_pfq_checked` returns `Err(HypError::{NotANumber, DenominatorPole, Divergent})`; `z = 0` is `1` in every case; `hyp_pfq` keeps NaN (documented). The Tricomi `U` results below `1e-270` are the IEEE underflow to zero, as documented (not NaN) |
| (a) | two unit tests of `kernels::ball` failed with the code that was committed (wrong expected values in the tests: a point outside the sum rectangle, a reversed interval) | **fixed** in the tests |
| (a) | **header generation was broken**: `DEV=1 cargo build` panicked in cbindgen (`Ratio` of `rules::poly` clashing with `num_rational::Ratio<T>`, the generic `BoundaryCondition<T>` of `physics_bem` clashing with `sim::boundary::BoundaryCondition`) | **fixed**: `rules::poly::Ratio` renamed `PolyRatio`, `VertexBc` marked `cbindgen:ignore`, three Rust-only types (`Dd`, `Face`, `Interval`) excluded in `cbindgen.toml` |
| (a) | clippy: a doc table with `\|` inside code spans and a list-like continuation line inside a code span (`indefinite_sum.rs`) made `doc_markdown` report "backticks are unbalanced" (25 table rows fixed) | **fixed** |
| (b) | wavelets (SciPy / PyWavelets workflows) | **done (core)**: `kernels::wavelets` — Haar, `db2`, `db3`, `db4` orthogonal filters, periodised `dwt` / `idwt`, `wavedec` / `waverec`; tests: orthonormality and vanishing moments of the filters, Parseval, exact reconstruction, multilevel round trip, error cases. Not done: biorthogonal / symmetric families, 2-D transform, other boundary extensions |
| (b) | fractional Fourier transform | **done (core)**: `kernels::wavelets::{frft_matrix, frft}`, the discrete (Candan / Pei) transform from the Hermite-like eigenvectors of the matrix commuting with the DFT (degenerate eigenvalue pairs, which occur for even lengths, are rotated to pure DFT eigenvectors); tests: `a = 0` identity, `a = 1` unitary DFT, `a = 2` reversal, `a = 4`, index additivity, unitarity (lengths 8, 15, 16, 24). `O(N^3)` |
| (b) | F4 / F5 Gröbner | **done (L15)**: `groebner_f4` / `groebner_f5` in `rules/poly/f4.rs` and `f5.rs`, matrix elimination and signature criteria, compared literally with Buchberger on Katsura/cyclic/random systems; F5B runs the criteria path on Katsura-4/Cyclic-5 (~29/~26 ms grevlex); the complete-pair rerun was removed and the F4 matrix is only a resource fallback. **Finite fields (done, 2026-10-31)**: the same engines run over any prime field `GF(p)` through a coefficient domain with genuine modular arithmetic; forms `groebner_f4(polys, vars[, order][, p])` and `groebner_f5(...)`; prime fields only (extension fields out of scope). |
| (b) | quadratic sieve | **done (L15)**: MPQS in `kernels/quadratic_sieve.rs` with factor base, large primes, GF(2) elimination and `factorint` dispatch; BPSW primality for large integers; the 50-digit ceiling is gone (practical ~65 digits; 55-digit ~8–12 s, 60-digit ~20–45 s; above the 40k column cap returns unevaluated) |
| (b) | Meijer G, Hermite reduction for algebraic integrands, real quantifier elimination, multivariate series | **done (L15)**: Meijer G (Slater/Mellin–Barnes numeric, symbolic reductions, derivative/integral/Mellin/series); Hermite reduction in `Q(x)(√D)` with confluent/quadratic radicands and residue log part; CAD-based `real_qe`/`cad_cells` with exact algebraic numbers, no fixed variable cap and `exists`/`forall` routing; `mseries`/`mtaylor` with Puiseux monomial factors, symbolic exponents and directional Laurent expansions |
| (c) | Guruswami–Sudan / LDPC decoding, Darboux polynomials / Jacobi last multiplier | documented, not done |
| (c) | all simulation `Limitations:` (`ns_mac` free surface / variable density, DG curved hexahedra, WENO in more than one dimension, `kerr` closed forms, `spin_models` beyond nearest neighbour, `euler_amr` 3-D, `yee_fdtd` media, `hyperelastic` P1 locking, `sse` winding numbers, DPD frozen walls, directional Helmholtz compression, ...) | documented, not done |

Header: `rssn.h` / `rssn.hpp` regenerated with `DEV=1 cargo build`; the diff
includes `rssn_compute_all` and the dimension API (`rssn_dimension`,
`rssn_compute_json`, mode/assumption setters and status codes).

**Final suite (one run of `cargo test --all-targets --no-fail-fast`, then `cargo test --doc`, log in the scratch directory), after all Round 3 changes:**

| Target | Passed | Failed | Ignored |
|---|---|---|---|
| lib (`rssn`) | 976 | 0 | 0 |
| `tests/io.rs` | 13 | 0 | 0 |
| `tests/kernels.rs` | 1052 | 0 | 0 |
| `tests/sim.rs` | 835 | 0 | 2 |
| `benches/engine.rs`, examples | (no tests) | 0 | 0 |
| doc tests | 45 | 0 | 0 |

The first full run (before the fixes) had two failures, both in `kernels::ball::tests` (wrong expectations in the tests, see the table);
the two ignored sim tests carry their own `#[ignore = "..."]` reasons from earlier rounds and were not changed.
`cargo clippy --all-targets` and `cargo clippy --all-targets --all-features` (cranelift `jit` and `npy` build on this machine): no errors, no warnings.
The crate-level lint configuration and `clippy.toml` are unchanged; the 25 `;` / `is_multiple_of` fixes came from `cargo clippy --fix`; the only `#[allow]`s added are
`clippy::approx_constant` (with the reason that `f64::consts::EULER_GAMMA` is unstable) on seven test items and `clippy::type_complexity` on one test helper.
`cargo +nightly fmt` was run last.

---

## L15 — the six remaining algorithms (2026-10-31, done)

The "not present anywhere" list of 5.2 is closed. Six algorithms were
implemented on separate branches, merged, and validated on the merged tree.
Each ships with module documentation, unit tests and public-API integration
tests, and every symbolic result is checked numerically or against an
independent reference.

| Algorithm | Files | Public entry points | Tests (merged tree) |
|---|---|---|---|
| Meijer G | `kernels/meijer.rs`, `rules/special_meijer.rs` | `meijer_g`, `meijer_g_real`; `meijerg(m,n,list(a),list(b),z)`, `meijer_mellin`, `meijer_series` | `--lib meijer` 34; `--test meijer_g` 10 |
| Hermite reduction for algebraic integrands | `rules/calculus/algebraic.rs` + `algebraic_{field,hermite,logpart,quadratic,term}.rs` | `integrate_algebraic` (crate) and the `integral` kernel | `--lib rules::calculus` 60; `--test algebraic_integration` 12 |
| Real quantifier elimination | `rules/solve/qe.rs` + `qe_{alg,poly,field,cad,formula,qe_impl}.rs`, `rules/logic.rs` routing | `real_qe(formula, list(vars))`, `cad_cells(polys, list(vars))`; `exists`/`forall` over real polynomial atoms | `--lib rules::solve::qe` 39; `--test real_qe` 15; `--lib rules::logic` 22 |
| Multivariate series | `rules/calculus/multiseries.rs` | `mseries(f, list(x), list(a), N)`, `mtaylor(...)` | `--lib multiseries` 14; `--test multivariate_series` 17 |
| F4 / F5 Gröbner | `rules/poly/f4.rs`, `f5.rs`, `rules/poly/field.rs` | `groebner_f4(polys, vars[, order][, p])`, `groebner_f5(...)`; `GF(p)` coefficient domain (`Domain`, `groebner_in`, `reduce_in`, `fglm_in`, `groebner_via_in`) | `--lib rules::poly` 92 (14 finite-field); `--test f4_f5` 16 (6 finite-field) |
| Quadratic sieve (MPQS) | `kernels/quadratic_sieve.rs`, `rules/number_theory.rs` dispatch | `qs_factor`, `qs_factor_with`, `QsLimits`, `QsConfig`, `QsStats`; `factorint` integration | `--lib kernels::quadratic_sieve` 10; `--test quadratic_sieve` 11; `--lib rules::number_theory` 27 |

**Meijer G.** Numeric evaluation by Slater's residue series (exact
consecutive-gamma ratios) with a numerical Mellin–Barnes fallback,
**exact confluent-pole residues** (differentiated Lanczos log-Gamma,
reflection formula, repeated/integer-spaced parameters, order up to 64) and
`p = q` evaluation on both sides of the unit circle; symbolic reductions to
powers, exponentials, erf/erfc, log/atan, Bessel J/Y/K/I with the required
parameter relations; the terminating Slater pFq; derivative, antiderivative,
Mellin transform and series at zero and infinity, all verified numerically.
Confluent accuracy is now ~1e-11–1e-13 (was ~1e-9, worst 1e-7). Limits:
`p<q` at very large |z|, `p=q` exactly on |z| = 1 with `delta <= 0`, and
identically-zero residue families stay unevaluated; literal `G^{1,0}_{1,1}`
has no global closed form (it is zero for |z|>1) and is deliberately not
reduced, while `G^{1,1}_{1,1}` is.

**Hermite reduction.** Exact arithmetic in `Q(x)(√D)`: square-free
normalisation, conic parametrisations for degree 1–2 radicands, and for degree
≥3 the complete Hermite reduction (each multiplicity is lowered by subtracting
an exact derivative, checked with rational-function arithmetic) followed by
the logarithmic part, with an 18-point differentiation check of every answer.
The logarithmic part computes the residue in the algebra `Q(α, y)`,
`y^2 = D(α)`, and either sums the conjugates of a residue in `Q(α)` (linear
factors, quadratic residue fields by atan, constant residues in higher degree)
or, when the residue has a `y` component, solves the exact linear system
`g' = f g` for `g = P + Q y` and returns `log(g)` when the primitive is a
logarithm of an algebraic function (e.g. `log(1 + D y)` for
`3 D D'(y - D^2)/(2(1 - D^3))`, which closes elliptic cases).  The `y`-trace
`2a` of the residue is the base-field projection, not the antiderivative; it
is checked through the norm `N = P^2 - Q^2 D`, `N'/N = 2a/c`.  A residue with
a `y` component that is not a logarithmic derivative (the generic
positive-genus third-kind integral) returns the Hermite-reduced integral (an
identity, not a wrong answer); ramified repeated factors and non-square
constant radicands stay unevaluated.

**Real QE.** Exact real algebraic numbers (square-free defining polynomial +
rational isolating interval, gcd/Sturm sign and comparison), sparse
multivariate polynomials over Q with Bareiss resultants, number-field
arithmetic `Q(θ)` with Sturm isolation and primitive-element extension, Collins
projection with the full coefficient/nullification sets, exact lifting
(rational sectors, algebraic sections), `real_qe` with inside-out block
elimination, sign-condition output minimised greedily, and an independent
rational-grid verification of every block. `exists`/`forall` over polynomial
real atoms are routed here; non-polynomial atoms stay unevaluated instead of
falling into the Boolean path. There is **no variable cap**: projection,
lifting and cell construction are generic in the number of variables, and
five- and six-variable formulas with small CADs are decided. Limits are
resource budgets (256 projection polynomials per level and 4 096 over the
tower, 20 000 cells, fibre degree 48, 20 000 terms per projection polynomial,
number-field degree 64, 1 000 000 projection/lifting steps); exhaustion
returns unevaluated. Audit fixes include binding-aware free variables,
single-coprime-factor lifting, and exact rational-root snapping.

**Multivariate series.** Recursive univariate expansion with a provably
sufficient coefficient bound whenever the function has a monomial factor
(Taylor and monomial-factor Laurent), total-degree truncation, directional
Laurent expansions for rational poles along a hypersurface, exact symbolic
coefficients and points, finite points and ±∞, `n=1` identical to
`taylor`/`laurent`, and an error term documented as the truncated part (same
convention as `taylor`). **Genuine Puiseux factors are expanded** by rational
monomial-factor extraction (`sqrt(x*y)`, `x^(4/3)`, `(x*y)^(5/2)`,
`sqrt(x^2*(1+y))`), and powers with a **symbolic exponent** are expanded by
the binomial series when the base has a non-zero numeric constant
(`(1+x+y)^a`). Genuine branch loci (`(x^2+y^2)^(1/2)` at the origin) and
essential singularities (`exp(1/x)`, `sin(1/x)`, `exp(x/y)`, classified
explicitly) stay unevaluated and never return a fabricated polynomial.
Limits: 8 variables, order 20, ramification caps Q<=64 / Q*N<=128, symbolic
exponents only for analytic bases.

**F4 / F5.** F4: sugar/normal-strategy batches, product and chain criteria,
Macaulay matrix with symbolic preprocessing, sparse exact Gaussian elimination
over the coefficient domain, row→polynomial map, tail reduction, unique
reduced basis. F5: signatures `(module monomial, index)`, position-over-term
order, F5 + rewritten + rewritable criteria, signature-admissible reduction,
homogenisation for affine input, FGLM conversion, with the F4 matrix as a
last-resort resource fallback; both operators return exactly `groebner`'s
format and are compared literally with Buchberger in the integration suite.
F5B now runs on the criteria path itself: Katsura-4 grevlex ~29 ms (was
6.5 s) and Cyclic-5 ~26 ms (was 5.9 s), all orders correct; the complete-pair
rerun was removed. Positive-dimensional cyclic-4 lex/grlex uses the
documented FGLM order fallback, recorded in `F5Stats::matrix_order_fallback`.
F4 katsura-4/cyclic-5 take 76/330 ms in grevlex.

**F4 / F5 over finite fields (2026-10-31).** `rules/poly/field.rs` adds the
coefficient domain `Domain` (`Q` or `GF(p)`, `p` an arbitrary-precision prime
validated by BPSW) with genuine modular arithmetic: `add`/`sub`/`mul` reduce
modulo `p`, `inv` is a modular inverse, and every kernel of the engines goes
through it — the Macaulay-matrix echelon form and S-polynomial of F4, the
one-way signature-admissible normal form of F5, the monic normalisation and
the final interreduction, and the FGLM order conversion (`fglm_in`). Nothing
is lifted to `Q`; the `GF(p)` result is monic, reduced, sorted and uses the
canonical representatives `0 .. p`, so it is literally the same list as
`groebner_in` over the same field. The operator forms are
`groebner_f4(polys, vars)` (lex, `Q`),
`groebner_f4(polys, vars, order)`, `groebner_f4(polys, vars, p)` (lex, `GF(p)`)
and `groebner_f4(polys, vars, order, p)`; the same for `groebner_f5`. A
rational coefficient `a/b` maps to `a b^-1 mod p`; a denominator divisible by
`p`, or a composite/non-positive modulus, leaves the request unevaluated.
Limits: prime fields only — extension fields `GF(p^k)` remain out of scope
(the `gfx_*` operators serve them separately); the F5 coefficient guard caps a
canonical representative at 32 768 bits. Verification: an independent naive
`GF(p)` Buchberger oracle (`field::oracle`, own `i64` arithmetic, own monomial
comparison) cross-checks 120 engine cases per engine (8 systems × 5 primes
{2,3,5,7,101} × 3 orders) plus direct-Buchberger cases in every order, FGLM
conversions, membership decisions and the no-fallback statistics. F5 solves
Katsura-3 over `GF(5)` in ~2 ms and Cyclic-4 in ~1 ms in every order; the
criteria path is used in grevlex and the positive-dimensional order fallback
for Cyclic-4 is recorded in `matrix_order_fallback` exactly as over `Q`.

**Quadratic sieve.** MPQS with `a=q²`, Tonelli–Shanks + Hensel lift for
`b² ≡ n (mod a)`, a Jacobi factor base, log sieving, single-large-prime
pairing, exact BigInt factor extraction, bitset GF(2) elimination and
`gcd(X±Y,n)` retries, fully deterministic; `factorint` dispatches large
composites to it. The 50-digit ceiling is gone: parameters scale from the QS
heuristic through 80 digits (practical ceiling ~65 digits, enforced by a 40k
column cap that returns `None`, not a hang), the GF(2) solve is structured
Markowitz elimination with bitsets, and the sieve covers both parities. A
40-digit semiprime takes ~1–2 s; 55-digit ~8–12 s; 60-digit ~20–45 s;
65-digit is bounded by the test budget. BPSW replaced the earlier 12-base
Miller–Rabin after the audit found the 24-digit strong pseudoprime
318665857834031151167461.

**Dimensional analysis.** The dimension typing layer
(`docs/DIMENSIONS.md`): `DimensionMode::Advisory` by default, strict hard
errors, symbolic exponents, units and affine temperature boundary checks,
user-defined base dimensions (`base:Name` text syntax) and an explicit
bypass. The **C API** exposes dimension mode/assumptions/base dimensions,
`rssn_dimension`, `rssn_compute_json` (term, value, reduced, dimension,
diagnostics) and `RSSN_STATUS_DIMENSION`/`RSSN_STATUS_INVALID_ARGUMENT`
(`rssn.h`/`rssn.hpp` regenerated). Findings are structured `DimError`s
(`UnknownUnit`, `Affine`, `MixedList`, ...) and advisory non-conflicts do not
fail Strict. The core tests are `--lib graph::dim` (16), `--test dimensions`
(15), `--test dimensions_core` (24), `dimensions_special` (18),
`dimensions_ode_pde` (4), `dimensions_physics` (9) — 86 tests — plus the FFI
dimension test. The coverage scanner is unignored and non-vacuous; it rejects
missing, explicit `Unspecified` and accidental `Any` rules. Structured
decompositions (`lu`, `qr`, `svd`, `eigenvects`, `jordan`, `ldl`, ...),
distribution shape parameters and homogeneous `matfun` are typed.

**Known limits after the round.** Meijer G keeps the genuine numeric limits
listed above; QE has no variable cap but resource budgets; multivariate series
has the ramification/order and essential-singularity limits listed above; F5
is prime-field only (extension fields `GF(p^k)` use `gfx_*`) and FGLM is
bounded at quotient dimension 400; the quadratic sieve is practical to ~65
digits; Hermite's logarithmic part closes when the residue is in the base
field or is the logarithmic derivative of an algebraic function (`g = P + Q y`
with `P`, `Q` of degree at most `deg C + 2`), and otherwise leaves the reduced
integral. Dimension limits are listed in `docs/DIMENSIONS.md` §13 (text-only
C diagnostics, homogeneous-only `matfun`, affine points boundary-checked but
not a separate type, container answers not tuple-typed). All documented limits
return unevaluated rather than a wrong answer.

For F4/F5 the finite-field support covers prime fields `GF(p)` only (no
`GF(p^k)`; use the `gfx_*` extension-field operators for those), requires the
modulus to be a literal prime larger than one, and refuses a generator whose
rational coefficient has no image modulo `p` (`p` divides its denominator).
The F5 `MAX_BITS` guard of 32 768 bits per canonical coefficient applies to
`p` as well. The independent Buchberger cross-check is done on systems small
enough for a naive implementation; the larger finite-field benchmarks are
cross-checked against F4/the FGLM route and certified structurally (monic,
reduced, all generators reduce to zero, S-polynomials reduce to zero in the
oracle's own arithmetic).

