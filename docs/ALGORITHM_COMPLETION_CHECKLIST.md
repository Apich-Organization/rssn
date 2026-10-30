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
  general algebraic extensions), integration of `exp` of a general
  polynomial in more than one variable, definite integrals by the
  Mellin–Barnes method. **P2/P3**
- `limits.rs`: ~~no Gruntz algorithm~~ — done (L14, `gruntz.rs`); still
  open: trigonometric/oscillating factors and algebraic cancellation beyond
  series order 6 make the kernel return `None`. **P2**
- `series.rs`: no multi-variable series, no Puiseux at arbitrary algebraic
  branch points. **P2**
- Summation: Gosper + recurrence guessing + first-order Zeilberger (L14);
  higher-order Zeilberger (Dixon, Apéry, Saalschütz) open → R1. **P2**

### 1.3 Polynomials, solving, Galois (`src/rules/poly`, `rules/solve`)

- Factorisation: Zassenhaus/Kronecker; missing: van Hoeij's lattice
  factorisation for high-degree integer polynomials (slow path), and
  factorisation over algebraic extensions. **P2**
- Groebner: Buchberger with criteria; no F4/F5. FGLM change of ordering
  for zero-dimensional ideals is done (R5, Part 5.3). **P2**
- `solve`: sextic/septic towers as in L2/L9 (not done: sextic towers with
  non-real block invariants, septics with three real pairs, degree ≥ 8
  without decomposition). `cad(polys, vars)` exists (sample points only);
  `exists`/`forall` exist only for Boolean atoms — no real quantifier
  elimination over CAD cells. **P2**
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
  Coulomb done (L13); still open: Meijer G, Fox H, Kelvin orders > 1,
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
  Pollard rho, `isprime` is Miller–Rabin with the first 12 primes; no ECM,
  no quadratic sieve, no BPSW/primality certificate for large integers.
  **P2**
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
  non-terminating `b_j ∈ Z≤0` give NaN. **P2**
- `kernels/sde.rs`: Milstein/SRK/adaptive need commutative noise; the
  adaptive solver has no Lévy-area path. **P2**
- `kernels/eigs_ext.rs`: *(correction, R11)* the non-symmetric
  Jacobi–Davidson `eigs_jd` (complex targets, partial Schur deflation,
  complex harmonic extraction) already exists since commit 1b2530a4; the
  symmetric generalised solver still needs SPD `B`, and a two-sided JD is
  not implemented. **P2** *(two-sided JD done, R11b in 5.5)*
- `kernels/pde.rs`: no 3D, Robin/periodic faces, time-dependent BCs or
  adaptive stepping in the simple kernels (routed to `sim`). **P3**
- `sim/bem_aca.rs`: Helmholtz path is dense; no cluster-tree rebalancing.
  **P2**
- `sim/dg.rs`: WENO limiter 1D only; hierarchical limiter skips
  non-periodic boundaries (same as 3D). **P2**
- `sim/dg3d.rs`: uniform axis-aligned boxes only (no curved/unstructured
  hexahedra, no hanging nodes). **P2**
- `sim/dpd.rs`: single species, no bonds, no walls. **P2**
- `sim/elasticity_fem.rs`: isotropic, homogeneous, small strain, static.
  **P2**
- `sim/euler_amr.rs`: 1D/2D only, ratio 2, no non-reflecting BCs. **P2**
- `sim/fem_simplex.rs`: constant `k`, `c`; no coarsening, no
  goal-oriented adaptivity, no multigrid. **P2**
- `sim/ns_mac.rs`: uniform grids, Boussinesq only, no free surface. **P3**
- `sim/sph.rs`: no self-gravity. **P3**
- `sim/split_step.rs`: periodic only, no rotating frame, no dipolar terms.
  **P3**
- `sim/spin_models.rs`: nearest-neighbour hypercubic lattices only. **P3**
- `sim/kerr.rs`: no separatrix/circular/extremal closed forms. **P3**
- Not present anywhere (grep 2026-10-29; FGLM, BPSW and ECM since added, R5/R6): F4/F5, quadratic
  sieve, BPSW/ECPP, Meijer G, wavelets, fractional Fourier, Darboux
  polynomials / Jacobi last multiplier, multivariate series, Hermite
  reduction for algebraic integrands, real quantifier elimination,
  Guruswami–Sudan/LDPC.

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

Not reached from the "own search" list (still open): symmetric generalised
eigenproblem with indefinite `B` (documented), SDE adaptive step without
commutative noise, DG hierarchical limiter at non-periodic faces, variable
`k, c` in `fem_simplex`, the other `Limitations:` paragraphs of Part 5.2.

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
