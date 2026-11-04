# Algorithm completion status

**This file is the current, prioritised status of what is not finished. It is not
a work log.** The detailed narrative of past rounds lives in
[`docs/archive/COMPLETION_LOG_2026-11-03.md`](archive/COMPLETION_LOG_2026-11-03.md).

**Ground truth.** Each module's `//!` documentation block states its own current
limitations. This file summarises and prioritises them; where the two disagree,
the module header wins and should be updated first. A limitation that is not in
a module header and not listed here is undocumented — treat that as a defect in
the documentation, not as a supported feature.

## How to read this

Priorities: **P0** = the tree is wrong or unverifiable (fix first); **P1** = a
named gap that blocks a whole class of problems; **P2** = coverage or accuracy
improvement; **P3** = research-grade or optional.

Definition of done for any item: implementation + module-header limitation text
reduced to what remains + tests next to the module's neighbours (inline
`#[cfg(test)]` or `tests/<area>/…` in the same style) + `cargo check --lib` green
+ the item's targeted test command green. There is no requirement to run a full
suite: project experience is that a local algorithm change is pinned by its own
tests and nothing else (a full-suite run has never caught a regression caused by
such a change). Run only the tests of the module you touched.

Resource discipline (one shared cargo target dir): prefer `cargo check --lib`
over building; run one cargo process at a time (builds serialise on the target
lock); never `cargo build`, `cargo clean`, `cargo update`, release builds, or
`cargo test` with no filter; `cargo fmt` needs nightly (`cargo +nightly fmt`),
stable must not be used.

## Open items, by area

### Simulation (`src/sim`)

| ID | Item | Prio | Where |
|---|---|---|---|
| SM-01 | Free surface, variable density derived from an equation of state, and a thermal field together with variable density | P3 | `ns_mac.rs` |
| SM-02 | Curved/unstructured hexahedra, hanging nodes, h/p adaptivity; WENO is per-axis averaged (not a genuinely multidimensional stencil) | P2/P3 | `dg.rs`, `dg3d.rs` |
| SM-03 | 3-D Berger–Colella refinement; characteristic non-reflecting patch boundaries; shallow water in 3-D; cheaper regridding | P2 | `euler_amr.rs` |
| SM-04 | Sharp dry-bed front; per-cell Manning roughness; AMR with sources; hydrostatic gas closure is locally isothermal, axis-aligned gravity, reflective walls only | P2 | `euler_nd.rs` |
| SM-05 | Anisotropic/tensorial ε, μ, σ; saturable/nonlinear media; PMC with Dey–Mittra and TF/SF; non-vacuum incident medium | P2/P3 | `yee_fdtd.rs` |
| SM-06 | P1 volumetric locking in near-incompressible bending; friction and self-contact; damping beyond generalised-α | P3 | `hyperelastic.rs` |
| SM-07 | Extremal spin, radiation reaction/self-force, separatrix and circular Kerr closed forms | P3 | `kerr.rs`, `relativity.rs` |
| SM-08 | P3M; >1PN terms; neighbour scheme for the direct N-body; SPH periodic images | P3 | `fmm.rs`, `hermite.rs`, `nbody.rs`, `sph.rs` |
| SM-09 | SSE winding-number (stiffness) estimator and general imaginary-time measurements; non-cubic (triangular/hexagonal/fcc) spin lattices | P2 | `sse.rs`, `spin_models.rs` |
| SM-10 | Impedance boundaries and directional/FMM-type compression in BEM; multiple species/polymers/multi-axis walls in DPD | P2/P3 | `bem_helmholtz.rs`, `bem_aca.rs`, `dpd.rs` |
| SM-11 | Anisotropy/tensor materials, plasticity/damage, elastodynamics in the FEM elasticity layer; tensorial or time-dependent diffusivity | P3/P2 | `elasticity_fem.rs`, `structured.rs` |

### Kernels (`src/kernels`)

| ID | Item | Prio | Where |
|---|---|---|---|
| KN-01 | Rader/prime-length FFT, odd-length real packing, threaded transforms | P2/P3 | `transforms.rs` |
| KN-02 | SDIRK DAE/mass-matrix path, banded/sparse Jacobian and factorisation reuse; dense-output order at fixed time | P2 | `ode_sdirk.rs` |
| KN-03 | Wavelet families beyond db10/sym6/bior4.4; complex/packets/lifting/3-D transforms; multilevel on non-dyadic shapes | P2/P3 | `wavelets.rs` |
| KN-04 | ILU(k)/threshold dropping/reordering; mesh-independent iteration counts; ILU/IC without pivoting | P2 | `krylov.rs` |
| KN-05 | SQP second-order check is a dense FD Hessian; no SOCP/SDP/integer constraints; augmented-Lagrangian multipliers asymptotic only | P3 | `optim.rs` |
| KN-06 | `1F1`/`U` in the turning region `|z| > Re b` (uniform Airy asymptotics); complex corners of `U`; `|ln z| ≥ 4.5` polylog | P2/P3 | `hypergeometric.rs`, `special.rs` |
| KN-07 | Adaptive rules on simplices/balls; no rigorous enclosure; 2-D oscillatory quadrature is Filon for linear phases with an adaptive Genz–Malik general-phase fallback (no 2-D Levin) | P2/P3 | `quadrature.rs`, `convergence.rs` |
| KN-08 | Multivariate distributions beyond Dirichlet/Wishart, robust covariance (MCD); Tucker/HOSVD/tensor train; higher-order or C²-monotone splines | P2/P3 | `stats.rs`, `tensor.rs`, `interpolate.rs` |
| KN-10 | No reverse-mode/Hessian-vector products; third and higher derivatives; interval rules for opaque calls | P2/P3 | `backend/tape.rs` |

### Rules (`src/rules`)

| ID | Item | Prio | Where |
|---|---|---|---|
| RU-01 | Zeilberger: orders above 6, symbolic parameters beyond the bounded polynomial ansatz, certificates with degree above the capped search, non-integer Gamma shifts; boundary conversion only for finite support | P2/P3 | `calculus/zeilberger.rs` |
| RU-02 | Limits: boundedness analysis is syntactic; series fallback stops at the extended order; no oscillatory-endpoint theory | P3 | `calculus/{limits,gruntz}.rs` |
| RU-03 | Definite integrals by Mellin–Barnes; integration of `exp` of a general multivariate polynomial; algebraic integration beyond the implemented logarithmic part | P2/P3 | `calculus/{integrate,algebraic}.rs` |
| RU-04 | van Hoeij lattice factorisation for high degree; factorisation over algebraic extensions; sextic towers with non-real block invariants; degree ≥ 8 without decomposition | P2/P3 | `poly`, `solve` |
| RU-05 | Number-field towers in formal/`nf_hyper`; Darboux search degree ≤ 6 and no symbolic parameters; generalised (higher-order) symmetries; conformal-mapping/Padé–Hankel resummation | P2/P3 | `ode/{formal,nf_hyper,nonlinear,reduce,resum}.rs` |
| RU-06 | Riemann problems for systems other than flat-bed shallow water (Euler systems); more than two interacting waves; Weyl–Titchmarsh for half-line/continuous spectra (regular finite intervals are done); spheroidal exterior radiation; non-separable curvilinear metrics | P2/P3 | `pde/{nonlinear,sl,spheroidal,curvilinear}.rs` |
| RU-07 | Character tables for arbitrary groups; Guruswami–Sudan list decoding and LDPC | P3 | `discrete`, `number_theory.rs` |

### Kernel core, dimensions, I/O, FFI

| ID | Item | Prio | Where |
|---|---|---|---|
| CO-01 | E-graph: witness-vs-certificate ball semantics are hardened (rigorous `Ball` ops, `approx_interval`, `EvalCost`, singular-adjacent soundness sampling); remaining: interval rules for `tan`/`atan2`, reverse-mode/Hessian-vector products | P2 | `src/graph/mod.rs`, `src/backend/tape.rs` |
| DM-01 | Dimension layer: codes exist but are not yet surfaced through the C API; `matfun` homogeneity; affine point typing; container typing; strict-input typing (the per-rewrite scan and complex-plane checks are closed) | P2/P3 | `docs/DIMENSIONS.md` §13 (authoritative), `src/dim`, `src/graph/dim*.rs` |
| IO-01 | Plot coverage beyond SVG contours/streamlines/quiver and PGM/PPM rasters (surface and 3-D views, PNG/JPEG encoder would need a dependency); further symbol coverage | P3 | `src/io` (module header lists the gaps) |
| FF-01 | C ABI remains JSON+bincode single-request handles (no per-function tensor ABI); session/handle lifetime; headers now regenerated for the bincode entry points | P2 | `src/ffi`, `src/api` |

### Cross-cutting

- Benchmarks only time the engine (`benches/engine.rs`); heavy kernels (homotopy,
  Radau, DG3D, FMM) have no committed performance baseline. **P2**
- `rules::ode::linear` / `kovacic` tests take 8–13 minutes; keep them out of
  routine verification filters. **P2**
- `cargo clippy --lib --tests` still reports 18 pre-existing diagnostics at HEAD
  under rustc 1.99.1 (toolchain drift in modules this status page does not
  touch: `kernel/*`, `kernels/{sparse,tensor}`, `rules/poly/{apart,groebner}`,
  `rules/calculus/tensor`; CI runs clippy with `continue-on-error`). The work
  of 2026-11-03/04 added none: every diagnostic that appeared on a line added
  in those rounds was fixed. **P2**

## Recently completed

**Round 7 — deep completion (2026-11-04).** Nine items closing documented
limitations in the modules that had a limitation record, plus four follow-ups in
the newly audited peripheral modules. Each ships with a module-header update and
targeted tests.

| Item | Where | Verification / measured evidence |
|---|---|---|
| 2-D WENO limiter (per-axis Zhong–Shu, direction averaged) | `sim/dg.rs` | 39 passed; 2-D advection order 3.19 and limited == unlimited bitwise; 2-D Sod positive; element-mean drift 8e-15; 1-D path bit-identical |
| Well-balanced hydrostatic gas gravity | `sim/euler_nd.rs` | 21 passed; isothermal atmosphere max |u| 5.4e-16/8.2e-16/1.2e-15 at 40/80/160 cells (was ~1e-2); pressure error 1.7e-15; free fall still exact |
| 2-D oscillatory quadrature (Filon; general-phase Genz–Malik fallback) | `kernels/quadrature.rs` | 24 passed; cos(1000x)cos(1000y) error 5.3e-22 with 68 evals (vs ~263k for a tensor rule); fixed a pre-existing 1-D Filon β-series sign error |
| Peripheral kernels audit + closures | `kernels/{stats,signal,coordinates,vector,ordering,random,robust,differential_geometry,geometric_algebra}.rs` | 9 limitations sections added; `NormalDist::new`, continuous inverse CDFs, `signal::ifft`, `blackman_window`, `jacobian_determinant`; 88 passed |
| Weyl–Titchmarsh m-function + spectral density | `rules/pde/sl.rs` | 14 passed with the Riemann work; `-p k cot(kL)` to 1e-8; atom weights; eigenvalue peaks vs `sl_eigen` ≤1e-6 |
| Exact shallow-water Riemann solver (Toro two-shock) | `rules/pde/nonlinear.rs` | Stoker reference to 1e-6; RH/fan identities; symmetry; agrees with `sim::euler_nd` (mass 1e-10, star 2–3 %, shock ≤3 cells) |
| `EvalCost` extractor + singular-adjacent soundness sampling | `src/graph/{extract,soundness}.rs` | `graph::` 147 passed; weights documented; scan of 261 standard rewrites 2.86 s, no refutations; discontinuity no longer falsely refuted |
| Per-rewrite dimensional scan + complex-plane checks | `src/graph/dim*.rs`, `tests/dimensions*.rs` | 261/261 rewrites exercised, 882 instantiations, 0 refutations; fixed `power_branch` typing and spurious contour conflicts; dimension suites 21+75+76 |
| Filled contours + PGM/PPM raster + quiver; batch bincode + regenerated headers | `src/io`, `src/ffi` | `--test io` 37 passed; `rssn_compute_all_bincode` tested against JSON; `rssn.h`/`rssn.hpp` +126/−6 each from `DEV=1 cargo build` |

Pre-existing bugs found and fixed along the way: the 1-D Filon series sign
(`+4θ⁴/105` → `−4θ⁴/105`), the `power_branch`/branch-wrapper dimension
inconsistency (found by the new scan), spurious dimension conflicts on
`count_zeros_poles`/`winding_number`/`poles`/`zeros_of`/`singularity`/
`branch_cut`, soundness false refutations at floor/ceil discontinuities, and the
stale checked-in C headers.

**Round 6 — peripheral-module limitation audits (2026-11-04).** Four areas that
had no limitations record were audited in place; each module header now carries
the prioritised list, and the tractable gaps were closed:

| Area | Closed |
|---|---|
| E-graph (`src/graph`) | rigorous outward-rounded `Ball` ops + `Graph::approx_interval` over tape interval enclosures; sign facts only from exact/rigorous evidence; witness radii propagated with ball arithmetic; extraction 2.25×/3.0× faster (class-indexed tables, measured) — `--lib graph::` 137 passed |
| Dimensions (`src/dim`, `src/graph/dim*.rs`, `docs/DIMENSIONS.md`) | `DimError::code()` (8 stable codes); `matexp` needs only a dimensionless product; dimensionless `z` required for contour/Cauchy operators; §13 renumbered — 20+74 lib and 92 integration dimension tests passed |
| I/O (`src/io`) | `plot_contour_2d`/`contour_lines_2d` (marching squares) and `plot_streamlines_2d`/`streamline_paths_2d` (RK4); LaTeX/Typst escaping; pretty-printer Greek completion — `--test io` 26 passed |
| C API (`src/ffi`, `src/api`) | status/error matrix + 4 silent-gap fixes; additive bincode compute path (`rssn_compute_bincode`); session-unit resolution in `value_of` — `--lib ffi::` 11 and `--lib api::` 36 passed; header regen follow-up documented |

**Round 5 — limitation completions across sim / kernels / rules (2026-11-03).**
Twelve items, each with module-doc update and targeted tests:

| Item | Where | Tests |
|---|---|---|
| ohmic conduction + PMC walls | `sim/yee_fdtd.rs` | 14 (`--test sim -- yee_fdtd`) |
| variable-coefficient diffusion | `sim/structured.rs` | 14 (`-- structured`) |
| Manning bed friction | `sim/euler_nd.rs` | 20 (`-- euler_nd`) |
| heterogeneous elastic materials | `sim/elasticity_fem.rs` | 22 (`-- elasticity`) |
| real FFT pair (`rfft`/`irfft`) | `kernels/transforms.rs` | 22 (`--test kernels -- transforms`) |
| SDIRK3 dense output + events | `kernels/ode_sdirk.rs` | 7 (`-- sdirk`) |
| multilevel 2-D wavelets | `kernels/wavelets.rs` | 18 (`-- wavelets`) |
| Krylov auto-preconditioner + SQP second-order check | `kernels/krylov.rs`, `kernels/optim.rs` | 58 (`-- krylov optim`) |
| Zeilberger order 6 + boundary terms | `rules/calculus/zeilberger.rs` | 14 fast batch (`--lib`), 13+1 ignored full module |
| limits with bounded oscillatory factors + deeper series | `rules/calculus/{limits,gruntz}.rs` | 14 fast batch |
| Darboux polynomials + Jacobi last multiplier | `rules/ode/{nonlinear,mod}.rs` | 14 fast batch |

Also fixed: a pre-existing broken wavelets moment test (absolute bound
unsatisfiable for db9/db10 — now scale-relative) and a pre-existing wrong limit
(`(sin x − x + x³/6)/x⁵` returned `oo`, now `1/120`). Details, numbers and the
remaining-limits list are in the archive, Part 11.

Earlier rounds (all closed unless listed above): L1–L15 (index-3 DAE benchmark,
sextic Galois test case, rigorous homotopy certification, Coulomb scattering,
cylinder drag, exterior time-varying sources, Neumann/Robin Green's functions,
shocks, eigenring/Kovacic fields, 3-D DG limiting, 3-D FEM refinement, graph
canonical labelling, special-function layer, Zeilberger+Gruntz, the six remaining
algorithms), R1–R11, B1–B8, the Part 6/7 parallelisation rounds, and Round 3/4
hygiene and correctness fixes. The archive has one section per round.

## Verification recipes

Run only the group you changed; each line is one command.

```powershell
cargo check --lib                 # type gate for any change
cargo check --lib --tests         # also checks every integration-test target

cargo test --test kernels -- transforms sdirk wavelets krylov optim
cargo test --test kernels -- quadrature stats signal coordinates
cargo test --test sim -- structured yee_fdtd euler_nd elasticity
cargo test --test sim -- dg
cargo test --test dimensions      # dimensions_core / _special / _ode_pde / _physics likewise
cargo test --test io
cargo test --lib -- graph::store graph::extract graph::soundness graph::facts
cargo test --lib -- graph::dim ffi:: api:: weyl_titchmarsh riemann
cargo test --lib -- hypergeometric
cargo test --test kernels -- homotopy_cert
cargo test --test sim -- ns_mac dg3d fem_simplex
```

Full-suite counts, when a round records them, are in the archive. Do not run
`cargo test` unfiltered on this machine.

## Starting a new round

1. Read the module headers in the area you intend to touch (this list is a
   summary, not a substitute).
2. Pick items that are tractable and testable; a documented unevaluated result
   is better than a plausible wrong one. Research-grade gaps stay documented.
3. Implement, update the module's `Limitations` text, add tests next to the
   module's neighbours, run the type gate and the targeted tests.
4. Append one row (and only the numbers that matter) to the archive log; add or
   remove the item here. Keep this file short — it is a status page.

## Archive

- [`docs/archive/COMPLETION_LOG_2026-11-03.md`](archive/COMPLETION_LOG_2026-11-03.md)
  — the full narrative log through Round 5 (Parts 0–11): item-by-item plans,
  implementation notes, measured accuracy numbers, discoveries and failures.
