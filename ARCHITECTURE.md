# RSSN Implementation Reference

[English](ARCHITECTURE.md) | [简体中文](ARCHITECTURE.zh-CN.md)

The [system architecture manual](https://aevia.apich.org/architecture/03-rssn-engine?lang=en) covers RSSN in the context of Aevia and the registry. This reference documents the kernel internals.

**rssn** is built around one idea: **every symbolic operation is an identity transformation**. Operators such as `diff(f, x)` or `solve(eq, x)` are rewritten by rule sets on a hash-consed DAG with e-graph equivalence classes until an extractor selects the best equal term for the target *phase* (a closed form or a number).

Non-identity computations (time stepping, discretized fields, particle systems) live in `sim`, alongside plain numeric `kernels` shared by both halves.

```
 Rust users ──────► api ──────────────────────────────────┐
 C / C++ / Python ► ffi (C ABI over api + sim) ───────────┤
 hosts, plugins ──► kernel (dynamic call ABI, manifest) ──┤
                                                          │ request term + target phase
               ┌──────────────────────────────────────────▼───┐
               │ graph: hash-consed DAG + e-classes           │  op registry, patterns,
               │ scheduler (tiers) · extraction · dimensions  │  facts, numeric balls
               └───────┬──────────────────────────┬───────────┘
          rewrites +   │                          │  numeric witnesses
          kernels      │                          │
               ┌───────▼───────┐   ┌──────────────▼─┐   ┌──────────────┐
               │ rules/<domain>│──▶│ kernels/<area> │   │ sim/ (data,  │
               │ (identities)  │   │ (plain numerics│   │ not terms)   │
               └───────────────┘   └────────────────┘   └──────────────┘
                      io/ (parse, LaTeX, Typst, pretty, plot) · backend/ (compile, JIT)
                      dim/ (units, constants, unit systems — facade over graph::dim)
```

## Module Map

| Module | Role |
|---|---|
| `graph` | Domain-free e-graph engine: `store` (hash-consed nodes, union-find classes, literal folding), `op` (open registry with `HEAVY`/`OPAQUE_ON_APPLY`/`COMMUTATIVE` flags, cost, binders), `pattern`/`subst`, `rule` (`RuleSet`, `Installer`), `schedule` (tiered Reduce → Normalize → Explore fixpoint with budgets/bans), `window` (local tree search), `extract` (`ClosedForm`, `SizeCost`), `facts` (assumptions, sign/integer analysis), `eval` (interval balls), `soundness`, and `dim*` (dimensional typing). |
| `rules` | One `RuleSet` per mathematical domain; `rules::standard()` installs all of them. |
| `kernels` | Plain numerical algorithms on `f64`, slices, and `ndarray`: quadrature, ODE/SDE/PDE (`pde.rs`, `pde/three.rs`), special functions, linear algebra (dense/CSR/BSR, Krylov, eigensolvers), FFT, optimization (`kernels::optimize`), finite fields, graphs, homology, Lie structures, `convex` (LP/QP/SOCP/SDP), `milp`, `glm`, `timeseries`, `survival`, `mcmc`, `control`, `meshgen` (2D/3D Delaunay), `cad` (STL/STEP, BSP CSG), `interpolate`, and `mp` (arbitrary precision). Checked APIs document validation; low-level constructors (e.g. `csr_from_triplets`) may panic on invalid shapes/indices ([docs/MATH_COVERAGE.md](docs/MATH_COVERAGE.md)). |
| `api` | `Session`, `Term<'s>` (`Copy` handle with overloaded operators), `Config`, `Answer`, `ComputeError`, `api::cache`, and `api::value` (`call_value`, `call_key`, `call_ho`, `fork`/`merge`, plugin registration). |
| `kernel` | Embedding surface: signature/type language, dynamic `Value` ABI, `CallCtx`/`OpError`, plugins/callbacks, `KernelManifest`, `bincode-next` codec, and `kernel::cdylib` (`plugin-cdylib`) ([docs/KERNEL.md](docs/KERNEL.md)). |
| `dim` | Public dimension/unit layer (`UnitDef`, `UnitRegistry`, CODATA constants, natural/Planck/atomic/geometrized systems), re-exporting `graph::dim` types. |
| `io` | Round-trippable text, LaTeX, Typst, Unicode pretty-printing, plotting, `.npy` (`npy`), `io::frame` (`DataFrame`, CSV, JSON), `io::parquet` (v1/v2 pages, nested lists/maps, Snappy/gzip/Zstd/LZ4), `io::hdf5` (superblocks 0/2/3, fractal heaps, v2 B-trees, chunked filters), and `io::netcdf` (CDF-1/2/5). Codecs treat file lengths/offsets as untrusted and return `Err` instead of panicking. |
| `backend` | `Backend` trait compiling closed forms to `Interpreter`, SSA tape IR, or Cranelift JIT (`jit`). `backend::device` provides typed kernel IR, CPU devices, wgpu (`gpu.rs`, `f32` and double-single `f64` in `ds.rs`), WGSL/GLSL/OpenCL/CUDA/Metal emitters, launch fusion (`plan.rs`), and array programs ([docs/BACKENDS.md](docs/BACKENDS.md)). |
| `sim` | Time stepping, fields, and particles: FDM/FEM/FVM/BEM, CFD, MD, SPH, spectral, multigrid, `fem_triangles`, `gas_dynamics`, `integrators`, Ising, and `sim::dimensional` (`[nx]`/`[nx, ny]`/`[nx, ny, nz]` kernels). `sim::scenario::run(name, json)` is the uniform entry point; `sim::schema::catalog()` describes all scenarios. |
| `ffi` | Small C ABI over `api`, `kernel::manifest`, and `sim::scenario`; cbindgen generates `rssn.h`/`rssn.hpp` (`bazel run //:refresh_headers`; [docs/BUILDING.md](docs/BUILDING.md)). |
| `prelude`, `constant` | `use rssn::prelude::*` convenience re-exports; build metadata and mathematical constants. |

### Rule Sets

Standard sessions share an installed rule template and immutable manifest for the current unit-registry revision while keeping private term graphs, native-object stores, and answer caches. `rules::standard()` installs:

- **Core algebra & analysis**: `arith`, `elementary`, `calculus`, `poly`, `solve`, `ode`, `linalg`, `sparse_symbolic` (modular determinants, Dixon lifting, Wiedemann), `holonomic` (D-finite equations/identities), `geometry`, `complex`, `number_theory`, `combinatorics`, `logic`, `special`, `stats`, `transforms`.
- **Numerical operators** (`numeric` sub-modules in `src/rules/numeric/`): `bridge.rs`, `decomp.rs`, `functions.rs`, `ode.rs`, `resample.rs`, `solvers.rs`, `special.rs`, `transform.rs`, `conic.rs`, `milp.rs`, `glm.rs`, `timeseries.rs`, `control.rs`, `survival.rs`, `mcmc.rs`, `mesh.rs`, `mp.rs` ([docs/MATH_COVERAGE.md](docs/MATH_COVERAGE.md)), plus generated dynamic-ABI families (`numerical_vector`, `numerical_stats`, `numerical_special`, `numerical_matrix`), `numerical_callback` ([docs/NUMERICAL_CALLBACKS.md](docs/NUMERICAL_CALLBACKS.md)), and `sde`.
- **Applied math & physics**: `variational`, `pde`, `functional`, `physics`, `physics_quantum`, `optimize`, `verify`.
- **Discrete & abstract algebra**: `discrete` (graphs, finite/permutation groups, homology, coding, crypto), `rewriting` (Knuth–Bendix), `lie`, `lie_structure`, `geometric_algebra`.
- **Units & utilities**: `units`, `radicals`, `cse`.

Each rule set defines **operators**, guarded **rewrites**, **definitions**, **kernels** (returning `Equal`, `Pinned`, `Approx`, or `Pass`), and operator **attributes** (dimensions, signatures, lowerings, integral tables).

Number-theoretic kernels use a two-tier design: a `u64` word-sized fast path (`kernels::modular` Montgomery/Barrett/NTT) and a generic exact path over the same `Domain` interface, differentially tested against each other. Heavy algorithms (Risch–Bronstein/Risch–Norman, Hermite reduction, Meijer G, Gosper, Zassenhaus/Kronecker factorization, Gröbner/F4/F5/FGLM/CAD, ODE/PDE solvers, quadratic sieve, Schreier–Sims/Todd–Coxeter, `kernels::permgroup`) are implemented as kernels that write results back as equalities.

## Phases

- **`Target::Symbolic`**: Extracts with `ClosedForm` (no heavy operators remaining). If none exists, returns the smallest term with `reduced = false`.
- **`Target::Numeric { tolerance }`**: Evaluates the best term, falling back to a kernel's numeric witness (interval ball). Approximate numbers are never mixed into the e-graph as exact.

## Sessions, Isolation, and Concurrency

Each `Session::compute` call clones the rule-configuration template graph, imports the request, runs the engine on the clone, and imports only the answer back. `Session` is single-threaded (`RefCell`); use `Session::fork` / `Session::merge` across threads, or `TermRef` (`Copy + Send + Sync`) in bindings.

- **Caching (`api::cache`, `CachePolicy`)**: Answers are memoized by request node, template, and result-relevant `Config` fingerprint (wall-clock timeouts are not cached). Opt-in **lemmas** cache closed forms of heavy sub-requests. Registering a plugin or unit clears both tables.
- **Batching (`Session::compute_all`)**: Snapshots hits/lemmas at entry, runs misses on worker clones (at most 2 private jobs per Rayon worker), and imports answers in input order. `execution` honors `RAYON_NUM_THREADS`; `kernels::reuse` pools bounded FFT plans and scratch buffers.
- **Incremental Matching & Indexing**: Classes track `Graph::touched` epochs so depth-`d` rewrites and window passes skip unchanged subgraphs (`RSSN_FULL_SEARCH=1` disables skipping). `Engine` builds a `TierIndex` once per program, lazily indexes operators for inspection, and memoizes `Graph::facts`, `looks_reduced`, and failed window cells by graph revision.
- **Parallel Tiers & Profiling**: `Engine::apply_tier` searches rewrite batches in parallel (from 512 candidate nodes) and commits serially in deterministic order. Opt-in two-phase kernels (`Kernel::two_phase`) prepare work, run `PreparedWork` in parallel, and commit in node order. Set `RSSN_PROFILE_RULES=1` for per-rule accounting (`graph::profile`, `examples/engine_timing.rs`).
- **Cancellation**: Runs terminate at a fixpoint, budget limit, or cooperative `rssn::cancellation` checkpoint.

## Embedding: `kernel` and `dim`

Hosts inspect `Session::manifest` (deterministic, hashable catalog of operators, signatures, units, constants, and capabilities) and invoke operators via `Session::call_value`, `call_key`, or `call_ho`. Every operator requires a typed signature or `DynamicReason` (`bazel run //:signature_audit`). Signatures are supplemented from test traces (`signatures.witnessed.json`) and reviewed overrides (`signatures.overlay.json`); `build.rs` requires both.

`kernel` and `dim` are modules rather than separate crates because they share the expression graph, operator table, and dimension arena without `pub(crate)` shortcuts. See [docs/KERNEL.md](docs/KERNEL.md), [docs/CALLBACKS.md](docs/CALLBACKS.md), and [docs/NUMERICAL_CALLBACKS.md](docs/NUMERICAL_CALLBACKS.md).

## Diagnostics

When an evaluator declines arguments, the engine leaves the call unreduced. `kernel::diag` surfaces the reason without affecting reduction:

- Evaluators call `diag::fail(kind, arg, ...)` (or `rules::numeric::bridge` helpers), recording `Diag { op, arg, kind, message }` in a **thread-local collector** only when active.
- `Session::apply_diag` and `Session::explain` re-run the unreduced call on the calling thread under `diag::capture`.
- Table-driven numeric operators (`rules::numeric::TABLES`) record exact argument/kernel findings; heavy operators derive `Arity`, `Shape`, `Symbolic`, or `Unsupported` via `diag::classify`. `From<Diag> for OpError` maps findings to `OpError::HostDiagnostic` (`rssn.diag.<kind>`).

## Arbitrary Precision

`graph::number::Number` holds only integers, rationals, or `f64` doubles. Arbitrary precision is a value-level facility in `kernels::mp`:

| Layer | Contents |
|---|---|
| `kernels::mp::real` | `Real` (`mant * 2^exp`) computed exactly in `BigInt` and rounded once per operation (`Round`). |
| `kernels::mp::funcs` | Transcendentals (`exp`, `ln`, trig/hyperbolic, `pow`, `cbrt`, `hypot`, etc.) returning rigorous `Ball` enclosures, plus cached integer-series constants (`pi`, `e`, `ln 2`, `sqrt 2`, golden ratio). |
| `kernels::mp::decimal` | Exact decimal I/O; renders a `Ball` only when both endpoints round identically. |

`rules::numeric::mp` (`mp_eval`, `mp_ball`, `mp_to_float`, etc.) carries values as **decimal string literals** in the term graph so high-precision values hash-cons and cache without altering `Value` or the engine. `mp_eval` raises precision via Ziv's strategy up to `MAX_BITS` (~9,800 decimal digits, exponent capped at `2^18`).

## Dimensions

`graph::dim` types every scalar term with exact rational or rational-function exponents (`x^n`, `x^(a+b)`, `x^(1/n)`), unifying unknown symbols and polymorphic plugin signatures (`dim_signature`). Every standard operator must attach a `DimRule`. `Advisory` mode records conflicts in `Answer::diagnostics`; `DimensionMode::Strict` raises `ComputeError::Dimension`; `Config::bypass_dimension_checks` disables checks ([docs/DIMENSIONS.md](docs/DIMENSIONS.md)).

## Compiled Evaluation and the JIT

`backend::current()` resolves `backend::with_backend` → `backend::set_default` → `Interpreter`.

- **Tape IR (`backend::tape`)**: Lowers terms to an SSA tape (constant-folded, DCE, structurally hashed) running over `f64`, dual numbers (forward/reverse autodiff and Hessian-vector products in `tape::autodiff`), and outward-rounded intervals. `Tape::eval_columns` processes up to 256-row blocks.
- **Cranelift JIT (`backend::jit`, feature `jit`)**: Compiles tapes to `fn(inputs, params, out)` cached by structural hash and manifest/target key; `TieredBackend` interprets first and compiles in the background. `unsafe` is confined to `abi` and the two call sites.
- **LIR (`backend::lir`, `backend::lir_jit`)**: Typed register machine with control flow for native row loops, `device::native`, and Aevia's `@jit` functions ([docs/BACKENDS.md](docs/BACKENDS.md)).
- **Callers**: Quadrature, infinite sums, `odeint`, `nminimize`, `chebyshev_approx`/`rational_approx`, numeric solvers, and `Term::compile`.

## C Interface

`src/ffi` (`term.rs`, `compute.rs`, `manifest.rs`) exposes a handle-and-JSON/`bincode-next` C ABI (`uint32_t` term handles per session, `rssn_term_apply`, `rssn_compute*`, `rssn_manifest_*`, `rssn_sim_run`). Guarded entry points catch unwinding panics under the `ffi` profile and return `RSSN_STATUS_PANIC`. See `examples/c/smoke.c` and [docs/BUILDING.md](docs/BUILDING.md).

## Extending rssn

- **New operator/domain**: Add a `RuleSet` via `Config::with` or `Session::with_rules`, including a `DimRule` and `Signature` (or `DynamicReason`).
- **External plugin**: Register via `Session::register` or `plugin-cdylib` ([docs/KERNEL.md](docs/KERNEL.md)).
- **Integral forms / Lowerings / Backends**: Extend `rules::calculus::IntegralTable`, attach `backend::Lowering`, or implement `backend::Backend`.

## Lint Policy

`src/lib.rs` denies `clippy::all`, `pedantic`, `nursery`, `unwrap_used`, `expect_used`, `indexing_slicing`, `arithmetic_side_effects`, `single_call_fn`, and `missing_docs`, with documented scientific-code allowances (`suboptimal_flops`, `cast_precision_loss`, `many_single_char_names`, `similar_names`, `tuple_array_conversions`, `needless_range_loop`, `unreadable_literal`, `excessive_precision`, `while_float`, and base cast lints). Local `#[allow]`s require a justification comment.

## Testing

- **Unit & soundness tests**: Per-module unit tests, `graph::soundness` (`no_standard_rewrite_is_refuted`), `graph::proptests`, and `graph::dim` coverage scanners.
- **Integration suites (`tests/`)**: `tests/kernels.rs`, `tests/sim.rs`, `tests/dimensions*.rs`, and `tests/coverage_*.rs` (`coverage_algebra`, `coverage_api_complex`, `coverage_numeric_a`, `coverage_numeric_b`). GPU tests in `src/backend/device/tests.rs` require `--features gpu-wgpu` and an adapter.
- **Generated bindings & benchmarks**: `bazel test //:generated_sources_check_test //:contracts_test` ([tools/rssn-api-gen/README.md](tools/rssn-api-gen/README.md)), `bazel run //:benchmark -- engine`. See [README.md](README.md#testing) for fast filtered runs.

## Where to Read Next

| Document | Subject |
|---|---|
| [README.md](README.md) | Features and quick start |
| [docs/KERNEL.md](docs/KERNEL.md) | Embedding kernel: manifest, dynamic ABI, plugins, units |
| [docs/DIMENSIONS.md](docs/DIMENSIONS.md) | Dimension algebra, inference, and embedding API |
| [docs/CALLBACKS.md](docs/CALLBACKS.md) · [docs/NUMERICAL_CALLBACKS.md](docs/NUMERICAL_CALLBACKS.md) | Host callbacks and numerical callback operators |
| [docs/BUILDING.md](docs/BUILDING.md) | Rust and C artifacts, headers, build metadata |
| [docs/MATH_COVERAGE.md](docs/MATH_COVERAGE.md) | Numerical edge-case conventions and coverage |
| [CONTRIBUTING.md](CONTRIBUTING.md) | Workflow and checks |
