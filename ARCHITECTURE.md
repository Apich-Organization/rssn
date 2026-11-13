# rssn architecture

rssn is a computer algebra and scientific computing library built around one
idea: **every symbolic operation is an identity transformation**.
Differentiating, integrating, solving or simplifying produces a term that is
*equal* to the request. Each of these is a heavy operator applied to
arguments, e.g. `diff(f, x)` or `solve(eq, x)`, and a rule rewrites that
operator into something cheaper. A single engine runs every rule set over a
hash-consed DAG of terms with e-graph equivalence classes. Then an extractor
picks the best equal term for the requested *phase*: a closed form, or a
number.

What is *not* an identity transformation (time stepping, discretised fields,
particle systems) lives in `sim`, next to the plain numeric `kernels` that both
halves share.

```
 Rust users ─────► api ─────────────────────────────────┐
 C / C++ / Python ► ffi (C ABI over api + sim) ─────────┤
 hosts, plugins ──► kernel (dynamic call ABI, manifest) ┤
                                                        │ request term + target phase
              ┌─────────────────────────────────────────▼───┐
              │ graph: hash-consed DAG + e-classes          │  op registry, patterns,
              │ scheduler (tiers) · extraction · dimensions │  facts, numeric balls
              └───────┬─────────────────────────┬───────────┘
         rewrites +   │                         │  numeric witnesses
         kernels      │                         │
              ┌───────▼───────┐   ┌─────────────▼──┐   ┌──────────────┐
              │ rules/<domain>│──▶│ kernels/<area> │   │ sim/ (data,  │
              │ (identities)  │   │ (plain numerics│   │ not terms)   │
              └───────────────┘   └────────────────┘   └──────────────┘
                     io/ (parse, LaTeX, Typst, pretty, plot) · backend/ (compile, JIT)
                     dim/ (units, constants, unit systems — facade over graph::dim)
```

## Module map

| Module | Role |
|---|---|
| `graph` | The engine, free of any mathematical domain. `store` holds hash-consed nodes, union-find classes and literal folding. `op` is an **open** operator registry with descriptors: arity, flags such as `HEAVY`, `OPAQUE_ON_APPLY` and `COMMUTATIVE`, cost, and binders. `pattern`/`subst` handle text patterns with guards; `rule` holds `RuleSet`, `Installer`, kernels and definitions. `schedule` runs the tiered fixpoint engine (Reduce → Normalize → Explore, with budgets and bans); `window` does destructive local search on tree windows. `extract` provides `ClosedForm` and `SizeCost`. `facts` stores assumptions and sign/integer analysis. `eval` evaluates interval-ball numerics. `soundness` checks every rewrite numerically. `dim`, `dim_expr`, `dim_infer` and `dim_signature` are the dimensional typing layer (below). |
| `rules` | One `RuleSet` per mathematical domain; `rules::standard()` installs all of them (list below). |
| `kernels` | Plain numerical algorithms on `f64`, slices and `ndarray`, with no terms involved: quadrature, ODE/SDE/PDE solvers, special functions, linear algebra (dense, CSR and block-sparse BSR with block-Jacobi/BILU(0) preconditioners, Krylov solvers, eigensolvers), FFT, optimisation, finite fields, graph, homology and Lie-structure routines and more. `kernels::pde` solves heat, wave, Poisson/Laplace and advection problems by finite differences in 1D, 2D and 3D on a unit box (`pde.rs` holds the 1D and 2D solvers and the `PdeProblem` dispatcher, `pde/three.rs` the 3D solvers, `Grid3D`/`Field3D` and the 1D Poisson/Laplace solvers); `kernels::interpolate` has tensor-product interpolation in any dimension. Rules call them to produce numeric witnesses, and users can call them directly. Kernels never panic on degenerate input (NaN, empty, huge, mismatched shapes); the conventions and the honest list of gaps are in [docs/MATH_COVERAGE.md](docs/MATH_COVERAGE.md). |
| `api` | `Session`, `Term<'s>` (a `Copy` handle with overloaded operators), `Config` (rule sets, target phase, bindings, assumptions, dimensions, budget), `Session::compute -> Answer`, `ComputeError`, and the answer/lemma cache (`api::cache`). `api::value` adds the dynamic-value entry points (`call_value`, `call_key`, `call_ho`, `fork`/`merge`, plugin registration). |
| `kernel` | The **embedding kernel**: everything a host needs to drive rssn from outside the crate. The signature/type language, the dynamic `Value` ABI, the calling convention (`CallCtx`, `OpError`), plugins and host callbacks, the self-describing `KernelManifest`, the deterministic `bincode-next` codec, and the optional shared-library plugin ABI (`cdylib`). See [docs/KERNEL.md](docs/KERNEL.md). |
| `dim` | The public dimension and unit layer: the one unit table (`UnitDef`, `UnitRegistry`), CODATA physical constants and the unit systems (natural, Planck, atomic, geometrized). Re-exports the dimension types that `graph::dim` implements. |
| `io` | Text parsing and printing (round-trippable), LaTeX, Typst, Unicode pretty-printing, markup helpers, plotting and rasterisation, `.npy` (feature `npy`). |
| `backend` | `Backend` trait: compiles a closed-form term into a callable function. The reference `Interpreter` ships with the crate, plus the tape IR and the optional Cranelift JIT. `backend::device` is the data-parallel axis: a typed kernel IR (elementwise core, index values, addressed loads, one fold), CPU devices, a wgpu device (`gpu.rs`; `f32`, and `f64` emulated by double-single pairs in `ds.rs`), source emitters for WGSL/GLSL/OpenCL/CUDA/Metal (`emit.rs`, `wgsl.rs`), a launch-fusion planner (`plan.rs`) and a library of array programs built on one launch primitive (`array.rs`) ([docs/BACKENDS.md](docs/BACKENDS.md)). |
| `sim` | Time stepping, fields and particles: FDM/FEM/FVM/BEM, CFD, MD (Verlet, BAOAB Langevin, Nosé–Hoover, cell lists), SPH, spectral, multigrid; P1/P2 finite elements on unstructured triangle meshes with adaptive refinement (`fem_triangles`); compressible gas dynamics with HLLC/MUSCL and an exact Riemann solver (`gas_dynamics`); symplectic and stiff integrators (`integrators`); Ising Metropolis and Wolff clusters; and the scenario models. `sim::dimensional` holds the dimension-generic kernels (flat `x`-fastest arrays with `dims = [nx]`, `[nx, ny]` or `[nx, ny, nz]`: Poisson iterations, transport, heat ADI, Schrödinger Crank–Nicolson, leapfrog wave) around which the dimension-specific names added to `physics_cfd`, `physics_fdm`, `physics_cnm`, `physics_mtm` and `classical` are thin wrappers; `physics_fvm` and `physics_fea` add their own entry points ([docs/MATH_COVERAGE.md](docs/MATH_COVERAGE.md), "Space-dimension coverage of simulation solvers"). `sim::scenario::run(name, json)` is the uniform entry point; `sim::schema::catalog()` describes every scenario (typed parameters with defaults, witnessed result fields, accepted method names), which is what hosts generate typed wrappers from. |
| `ffi` | A deliberately small C ABI over `api`, `kernel::manifest` and `sim::scenario`. cbindgen generates `rssn.h`/`rssn.hpp` from it (`DEV=1 cargo build --features dev-headers`; see [docs/BUILDING.md](docs/BUILDING.md)). |
| `prelude`, `constant` | `use rssn::prelude::*` re-exports the names most programs need; `constant` holds build metadata (version, commit, target) and mathematical constants. |

### Rule sets

`rules::standard()` installs all of these (a rule set declares the sets it
depends on):

- **Core algebra and analysis:** `arith`, `elementary`, `calculus`, `poly`, `solve`, `ode`, `linalg`, `sparse_symbolic` (sparse matrices with exact or symbolic entries: sparsity-preserving elimination, modular determinants, Dixon lifting, Wiedemann), `holonomic` (D-finite functions: differential equations, recurrences and a decision procedure for identities), `geometry`, `complex`, `number_theory`, `combinatorics`, `logic`, `special`, `stats`, `transforms`.
- **Numerical operators** (value-level operators backed by `kernels`): `numeric` (`src/rules/numeric.rs` registers the operators and the sub-modules in `src/rules/numeric/` hold them: `bridge.rs` is the table-driven bridge that decodes literal numbers, nested lists and method symbols into kernel arguments and encodes the answer (robust statistics, windows and spectra, one-level wavelets, tensors, sampled functions, finite-difference solvers); `decomp.rs` rank-revealing and multilinear decompositions; `functions.rs` operators taking expressions (`oscillatory_integral`, `anderson_fixed_point`, `newton_cg_minimize`); `ode.rs` initial-value problems with a named integrator; `resample.rs` Hyndman–Fan quantiles and bootstrap intervals; `solvers.rs` sparse Krylov and multigrid solves, nested dissection and homotopy; `special.rs` complex special functions; `transform.rs` multilevel and 2D wavelets and the fractional Fourier transform; the list is in [docs/MATH_COVERAGE.md](docs/MATH_COVERAGE.md), "Kernels exposed as operators"), and the dynamic-ABI families `numerical_vector`, `numerical_stats`, `numerical_special`, `numerical_matrix` (generated bindings, see `tools/rssn-api-gen`), plus `numerical_callback` (host callbacks, [docs/NUMERICAL_CALLBACKS.md](docs/NUMERICAL_CALLBACKS.md)) and `sde` (noise primitives and scheme metadata).
- **Applied mathematics and physics:** `variational`, `pde`, `functional`, `physics`, `physics_quantum`, `optimize`, `verify`.
- **Discrete and abstract algebra:** `discrete` (graphs, finite and permutation groups, character tables, point groups, topology and homology, coding, cryptography, fractals), `rewriting` (Knuth–Bendix), `lie` and `lie_structure` (Lie algebras, root systems, representations), `geometric_algebra`.
- **Units and utilities:** `units`, `radicals`, `cse` (common-subexpression extraction).

A rule set contains:

- **operators**, registered by name, so a plugin can add its own;
- **rewrites**: text patterns with guards such as `is_const(?a, ?x)` or `positive(?a)`;
- **definitions**: `name(a, b) := body`, expanded by a kernel;
- **kernels**: Rust functions over e-classes that return `Equal(term)`, `Pinned(term)`, `Approx(ball)` or `Pass`;
- **attributes** on operators: dimension rules, signatures, backend lowerings, integral tables.

Heavy algorithms are kernels. Examples are Risch–Norman integration and the substitution/reduction-formula integrator, Hermite reduction for algebraic integrands in `Q(x)(√D)`, Meijer G evaluation and symbolic reduction, Gosper summation and recurrence guessing, univariate (Zassenhaus) and multivariate (Kronecker) factorisation, Gröbner bases (Buchberger, F4, F5), FGLM, CAD and real quantifier elimination, the equation and inequality solver with its heuristic sub-solvers, the ODE solver (classification, variation of parameters of any order, Lie point symmetries, reductions, Fourier and Laplace methods), the PDE solver (characteristics in any dimension, separation of variables, Green's functions, conservation laws and symmetries), multivariate Taylor/Laurent series, the quadratic sieve, and Schreier–Sims / Todd–Coxeter for groups, with the structural permutation-group algorithms (backtrack search for centralisers, normalisers and set stabilisers, conjugacy classes, Sylow subgroups, block systems) in `kernels::permgroup`. Each kernel writes its result back as an equality.

## Phases

`Target::Symbolic` extracts with `ClosedForm`: no heavy operator may remain. When no closed form exists, the cheapest term by size is returned with `reduced = false`.

`Target::Numeric { tolerance }` evaluates the best term instead. If that fails, it uses a numeric witness (a ball) that a kernel produced, for example a quadrature result. Numbers are never mixed into the e-graph as if they were exact.

## Sessions, isolation and concurrency

A `Session` owns a plain term store. Each `compute` call works as follows:

1. Clone the template graph for the requested rule configuration.
2. Import the request into the clone.
3. Run the engine on the clone.
4. Import only the answer back into the store.

Requests therefore never interfere with each other. A session itself is **single-threaded** (its store is a `RefCell`, so `Session` is not `Sync`): use one per thread, move work between threads with `Session::fork`, and bring results back with `Session::merge`. A `Term<'s>` borrows its session; for language bindings `kernel::value::TermRef` (a session id plus a node id, `Copy + Send + Sync`) is refused by any session but the one that minted it.

Every distinct set of rule sets used by a `Config` installs a *template* (a copy of the graph plus an engine) the first time it is seen; the session keeps it for reuse.

A session caches what it computed (`api::cache`, `CachePolicy`). **Answers** are memoised by request node (the store is hash-consed, so equal terms are the same node), template and a fingerprint of every result-relevant `Config` field; a hit equals a fresh run, and runs stopped by the wall clock are not stored. **Lemmas** (opt-in) record the closed forms that heavy sub-requests reached, keyed by the configuration without its budget, and assert them in later jobs that contain the same sub-request. Registering a plugin or unit clears both tables; `Session::cache_stats` reports hits.

`Session::compute_all(&terms, &config)` runs several requests of one configuration at once. `compute` is split into `prepare` (clone the template and import the request, sequential because the store is borrowed), `run_job` (engine run and choice of the answer on the private graph; this is the parallel part, one rayon task per request) and `finish` (import the answer into the store, sequential and in input order). The result, including the store contents, equals a loop of `compute` calls and does not depend on the thread count; the engine, graphs and `Config` are `Send`/`Sync` (checked by compile-time assertions in `api`).

Matching is incremental. Every class carries the epoch of its last change (a new member, a union, new analysis data, a witness or a pin: `Graph::touched`), and the scheduler advances the epoch after each search phase. A rewrite whose left-hand side has depth `d` and whose guards read only the bound classes (`free_of`, `number`, `symbol`, `same`) skips a candidate node when no class within `d` child steps changed since the rewrite's last complete search and the node was already a candidate then; rewrites with fact-based or custom guards, and banned rewrites coming back, search in full. Window passes likewise re-run only on terms one of whose classes changed. The matcher itself allocates only for completed matches (`pattern::MatchScratch`). `RSSN_FULL_SEARCH=1` turns the skipping off for differential testing.

Inside one run, `Engine::apply_tier` searches each maximal run of consecutive rewrite rules of a tier against the frozen graph (in parallel from 512 candidate nodes, each rule keeping its own match limit, fuel and ban decision) and applies the builds and unions serially in rule order, then match order (egg's scheme); the batched semantics are used on both paths, so results do not depend on the cut-off or the thread count. Kernels keep their position and run sequentially, except kernels that opt in to **two-phase** execution (`Kernel::two_phase`, `prepare`, `commit`): the scheduler prepares every eligible node first, runs the returned `PreparedWork` items (plain data, no graph) in parallel and commits in node order. `groebner(..)` requests use it. Window optimisation stays sequential: it mutates the graph. Parallel sites are listed in `src/kernels/par.rs` and are kept only where they beat the sequential path (threshold about 1.1×).

There is no cancellation and no streaming: a run ends when it reaches a fixpoint or exhausts its `Budget` (`max_iterations`, `max_nodes`, optional `time_limit`).

## Embedding: `kernel` and `dim`

A host (a language binding, a code generator such as Aevia) never needs crate internals. It reads the **Kernel Manifest** (`Session::manifest`), a deterministic, hashable description of every operator, signature, overload, unit, constant, unit system and build capability, and calls operators by name or pinned `OpKey` with `Value`s (`Session::call_value`, `call_key`, `call_ho`). Every operator is typed or carries an explicit `DynamicReason`; `cargo run --bin signature-audit` fails otherwise. Beyond what a registration states, operators are typed from *evidence*: the engine logs the sorts of requests and answers while the test suite runs, and `signature-mine` turns the log into replayed signatures (`signatures.witnessed.json`, `SignatureSource::Witnessed`; [docs/KERNEL.md](docs/KERNEL.md) §3a). `signatures.overlay.json` is the version-controlled place where a human corrects an inferred signature. `build.rs` requires both files and fails with a message naming the missing one.

`kernel` and `dim` are modules rather than the `rssn-api` and `rssn-dim` crates a downstream project asked for, because both halves share the expression graph, the operator table and the dimension arena; a crate boundary would force a facade over half of `graph`. The interfaces have no `pub(crate)` shortcuts, so a later split is a packaging change. The full contract is in [docs/KERNEL.md](docs/KERNEL.md); host callbacks are in [docs/CALLBACKS.md](docs/CALLBACKS.md) and [docs/NUMERICAL_CALLBACKS.md](docs/NUMERICAL_CALLBACKS.md).

## Diagnostics: why a call did not reduce

An operator whose evaluator rejects its arguments (a shape mismatch, a tolerance out of range, a singular system, incompatible Neumann data) does not fail: the engine leaves the call as an unreduced term, which is the correct answer for a rewrite system and no help to a person. `kernel::diag` is the explanation channel, designed so that nothing about reduction changes. An evaluator keeps its `Option<Data>` signature and, on a failure it can name, calls `diag::fail(kind, Some(arg), format_args!("..."))` (or the helpers of `rules::numeric::bridge`: `bad`, `require(..)?`, `.reason(kind, arg)` on a kernel `Result`), which records a `Diag { op, arg, kind, message }` in a **thread-local collector and returns `None`**. With no collector installed the call does nothing (the message is not even formatted), so ordinary reduction pays nothing and the engine, the rewrite scheduler and its worker threads never see the channel. A host asks with `Session::apply_diag(name, args, ctx) -> Result<Value, Diag>` (the same value as `call_value` on success; on a refusal the finding) or `Session::explain(term) -> Option<Diag>`: both re-run the evaluator of the unreduced call on the *calling* thread under `diag::capture`, so the answer is the reason the evaluator itself recorded, independent of which worker thread the engine used. Only the table-driven numeric families (`rules::numeric::TABLES`) record reasons today; a heavy operator without any yields `DiagKind::Unspecified`. `From<Diag> for OpError` maps a finding to `OpError::HostDiagnostic` with the stable code `rssn.diag.<kind>`.

## Dimensions

Dimensions are a typing layer over the term graph (`graph::dim`): every scalar term has a dimension, every operator states how the dimensions of its arguments determine its value's dimension, and `Session::compute` reports the dimension of the answer. Base dimensions are interned and extensible, exponents are exact rationals or rational functions over atoms (so `x^n`, `x^(a+b)` and `x^(1/n)` are supported), and unknown symbols get dimension variables that are solved by unification. `dim_signature` gives plugin operators polymorphic dimension signatures (monomials over signature variables, solved per use by the generic `signature_rule`).

A `DimRule` is attached to an operator like `Partials` (`i.dim(op, DimRule::Same(0))`); a scanning test requires every standard operator to carry one. The default mode is `Advisory`: dimensions are inferred and conflicts reported in `Answer::diagnostics` without changing the answer. `DimensionMode::Strict` turns a conflict into `ComputeError::Dimension`; `Config::bypass_dimension_checks` is the explicit escape hatch. Unit expressions are resolved by the `units` rule set through a resolver installed on `quantity`, so `graph` does not depend on `rules`. [docs/DIMENSIONS.md](docs/DIMENSIONS.md) is the design document and lists the known limits.

## Compiled evaluation and the JIT

Every compiled evaluation goes through `backend::current()`: the reference `Interpreter` by default, replaced globally with `backend::set_default` or for one scope with `backend::with_backend`. New backends implement `backend::Backend` and are passed to `Term::compile_with`.

- **Tape IR** (`backend::tape`). Terms are lowered to an SSA tape that undoes the canonical form (subtraction, division, integer powers by square-and-multiply, `sqrt`/`cbrt`), folds constants, removes dead code and is hashed structurally. The same tape runs over `f64`, dual numbers (forward-mode gradients, plus reverse gradients and Hessian-vector products in `tape::autodiff`) and outward-rounded intervals (rigorous enclosures, which also decide signs of constant classes in `graph::facts`).
- **Cranelift JIT** (`backend::jit`, feature `jit`). The tape is compiled to machine code with the signature `fn(inputs, params, out)`. Compiled functions are cached by structural hash and manifest/target key; `TieredBackend` interprets first and switches to machine code compiled in the background. Operators can declare a `backend::Lowering` (an intrinsic, an `extern "C"` function or a gather); the special functions do, others go through a generic evaluation shim. The JIT's `unsafe` code is confined to its `abi` module and the two call sites that enter generated code.
- **LIR** (`backend::lir`, `backend::lir_jit`). A typed register machine with basic blocks, calls, checked integers, bounds-checked slices and a *bail* exit, compiled by Cranelift or interpreted with identical results. It is the common target of everything that needs control flow: batch evaluation of tapes in a native row loop, native CPU device kernels (`device::native`), and Aevia's `@jit` functions ([docs/BACKENDS.md](docs/BACKENDS.md), "LIR").
- **Column-wise interpretation.** `Tape::eval_columns` runs each instruction over a block of 256 rows; `TapeFunction` batches use it.
- **Users.** Numeric quadrature, infinite sums, `odeint` (with a compiled symbolic Jacobian for the stiff Radau IIA fallback), `nminimize` (compiled objective and gradient), `chebyshev_approx`/`rational_approx`, the numeric solver and `Term::compile` all go through the current backend, so turning on the JIT speeds up every one of them.

## C interface

The C interface lives in `src/ffi` (`term.rs`, `compute.rs`, `manifest.rs`). It is a JSON-and-handles ABI, not a typed one: terms are `uint32_t` handles within a session, answers cross as JSON text, the `RssnAnswer`/`RssnTensor` structs, or deterministic `bincode-next` buffers. There are no per-operator entry points; a host builds the request with `rssn_term_apply` and computes it.

| Area | Functions |
|---|---|
| Sessions | `rssn_session_new`, `rssn_session_free` |
| Building terms | `rssn_term_sym`, `rssn_term_int`, `rssn_term_rational`, `rssn_term_float`, `rssn_term_parse`, `rssn_term_apply(op_name, args, n)` |
| Configuration | `rssn_config_new`, `rssn_config_free`, `rssn_config_symbolic`, `rssn_config_numeric`, `rssn_config_bind`, `rssn_config_assume`, `rssn_config_budget` |
| Dimensions | `rssn_config_dimension_mode`, `rssn_config_dimension`, `rssn_config_define_base_dimension`, `rssn_config_bypass_dimensions`, `rssn_dimension`, `rssn_dimension_count`, `rssn_dimension_code` |
| Computing | `rssn_compute` (fills `RssnAnswer`), `rssn_compute_json`, `rssn_compute_bincode`, `rssn_compute_all` / `rssn_compute_all_bincode` (n requests in parallel; per-request answers and status codes), `rssn_simplify` |
| Reading results | `rssn_term_to_f64`, `rssn_term_eval`, `rssn_term_tensor_data` / `rssn_tensor_free` |
| Output | `rssn_term_to_string`, `rssn_term_to_latex`, `rssn_term_to_typst`, `rssn_term_to_pretty` |
| Manifest | `rssn_manifest_json`, `rssn_manifest_hash` (describe the *standard* build) |
| Simulations | `rssn_sim_run(name, json, &out)` |
| Memory, errors, version | `rssn_string_free`, `rssn_bytes_free`, `RssnStatus` codes, `rssn_last_error`, `rssn_version` |

No call unwinds: a panic becomes `RSSN_STATUS_PANIC`. The module documentation of `ffi` lists the ABI's limits in priority order. See `examples/c/smoke.c` for a complete client and [docs/BUILDING.md](docs/BUILDING.md) for building the static and shared libraries.

## Extending rssn

- **New operator or domain.** Write a `RuleSet` containing operators, rewrites, definitions and kernels. Add it with `Config::with`, or with `Session::with_rules` for a session. Every operator needs a dimension rule and a signature (or a `DynamicReason`); the soundness, dimension-coverage and signature-audit tests check them.
- **Plugin from outside the crate.** Register a `Plugin` through `Session::register`, or load a shared library with feature `plugin-cdylib` ([docs/KERNEL.md](docs/KERNEL.md)).
- **New integral forms.** Register extra antiderivative functions in the `rules::calculus::IntegralTable` attribute of the `integral` operator.
- **New backend.** Implement `backend::Backend` and pass it to `Term::compile_with`, or install it with `backend::set_default`.
- **Faster special functions.** Attach a `backend::Lowering` attribute to an operator so compiled code calls it directly.

## Lint policy

`src/lib.rs` applies the strict configuration: `clippy::all`, `pedantic`, `nursery`, `unwrap_used`, `expect_used`, `indexing_slicing`, `arithmetic_side_effects`, `single_call_fn` and `missing_docs` are denied (the later `allow(clippy::restriction)` group keeps the individual restriction lints at their group level); `dead_code`, `unsafe_code` and `warnings` warn. A few lints are deliberately allowed for scientific code, each with its reason next to it in `lib.rs`:

| allowed | why |
|---|---|
| `suboptimal_flops` | `mul_add` is a slow libm call without hardware FMA, and fused rounding makes results differ between targets |
| `cast_precision_loss` | counts and indices become `f64` throughout numerics and stay far below 2^53 |
| `many_single_char_names`, `similar_names` | formulas keep the notation of the literature |
| `tuple_array_conversions` | false positives on arrays built from separately computed bindings |
| `needless_range_loop` | stencils, tableaux and matrices index several arrays in step |
| `unreadable_literal`, `excessive_precision` | quadrature nodes and coefficients are copied digit for digit from their references |
| `while_float` | time-stepping loops advance a floating-point clock |
| `cast_possible_truncation`, `cast_sign_loss`, `cast_possible_wrap` | from the base configuration |

Local `#[allow]`s are used sparingly and carry a comment.

## Testing

- **Unit tests** sit next to each rule set and kernel. Each rule set's tests check identities symbolically and, where applicable, numerically.
- **`graph::soundness`** evaluates both sides of every rewrite of every standard rule set at sampled points (`rules::tests::no_standard_rewrite_is_refuted`).
- **`graph::proptests`** runs property tests of the kernel.
- **Dimension tests** (`graph::dim` unit tests, `tests/dimensions*.rs`) include the operator-coverage scanner and a per-rewrite dimension scan.
- **Integration tests** live in `tests/` (symbolic batteries, kernels, sim, io, callbacks and numerical ABI, JIT and tape contracts, plugin loading, signature overlay). `tests/kernels.rs` and `tests/sim.rs` are single test binaries that include the files under `tests/kernels/` and `tests/sim/` as modules (`dimension_coverage.rs` in each pins the 1D/2D/3D families against analytic solutions and lower-dimensional results). The `tests/coverage_*.rs` files exercise the public API only: `coverage_algebra` (Lie classification, exact homology), `coverage_api_complex` (complex analysis, the session cache, printers, a parser/evaluator differential test), `coverage_numeric_a` (Fresnel/Faddeeva, low-rank kernels, structure-preserving integrators) and `coverage_numeric_b` (BSR, exact sparse linear algebra). The device tests in `src/backend/device/tests.rs` run on a real GPU only with `--features gpu-wgpu` and an adapter. Some suites are expensive; how to run them cheaply is in [README.md](README.md#testing).
- **Generated bindings:** `cargo run -p rssn-api-gen -- --bindings-check` (and `--structured-check`, `--contracts-check`) keep the checked-in numerical ABI fragments in step with their policy tables ([tools/rssn-api-gen](tools/rssn-api-gen/README.md)).
- **Benchmarks:** `cargo bench --bench engine`.

## Where to read next

| Document | Subject |
|---|---|
| [README.md](README.md) | Features and a quick start |
| [docs/KERNEL.md](docs/KERNEL.md) | The embedding kernel: manifest, dynamic ABI, plugins, units, structured results |
| [docs/DIMENSIONS.md](docs/DIMENSIONS.md) | The dimension layer: design, API, coverage, limits |
| [docs/CALLBACKS.md](docs/CALLBACKS.md), [docs/NUMERICAL_CALLBACKS.md](docs/NUMERICAL_CALLBACKS.md) | Host callbacks and the numerical callback operators |
| [docs/BUILDING.md](docs/BUILDING.md) | Rust and C artifacts, headers, build metadata |
| [docs/MATH_COVERAGE.md](docs/MATH_COVERAGE.md) | Edge-case conventions of the numerical kernels, where each area is tested, known mathematical gaps |
| [CONTRIBUTING.md](CONTRIBUTING.md) | Workflow and checks |
