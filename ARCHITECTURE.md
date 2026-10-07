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

```
              ┌─────────────── api ───────────────┐      ┌──── ffi ────┐
 Rust users → │ Session · Term · Config · compute │  ←── │ C ABI, <1k  │ ← C / C++ / Python
              └────────────────┬──────────────────┘      └─────────────┘
                               │ request term + target phase
              ┌────────────────▼──────────────────┐
              │ graph: hash-consed DAG + e-classes │  op registry, patterns,
              │ scheduler (tiers) · extraction     │  facts, numeric balls
              └───────┬───────────────────┬───────┘
         rewrites +   │                   │  numeric witnesses
         kernels      │                   │
              ┌───────▼───────┐   ┌───────▼────────┐   ┌──────────────┐
              │ rules/<domain>│──▶│ kernels/<area> │   │ sim/ (data,  │
              │ (identities)  │   │ (plain numerics│   │ not terms)   │
              └───────────────┘   └────────────────┘   └──────────────┘
                     io/ (parse, LaTeX, Typst, pretty, plot) · backend/ (compile)
```

## Layers

| Module | Role |
|---|---|
| `graph` | The kernel. `store` holds hash-consed nodes, union-find classes and literal folding. `op` is an **open** operator registry with descriptors: arity, flags such as `HEAVY`, `OPAQUE_ON_APPLY` and `COMMUTATIVE`, cost, and binders. `pattern`/`subst` handle text patterns with guards; `rule` holds `RuleSet`, `Installer`, kernels and definitions. `schedule` runs the tiered fixpoint engine (Reduce → Normalize → Explore, with budgets and bans). `extract` provides `ClosedForm` and `SizeCost`. `facts` stores assumptions and sign/integer analysis. `eval` evaluates interval-ball numerics. `soundness` checks every rewrite numerically. |
| `rules` | One `RuleSet` per mathematical domain. See the list below. |
| `kernels` | Plain numerical algorithms on `f64`, slices and `ndarray`, with no terms involved: quadrature, ODE/PDE solvers, special functions, linear algebra, FFT, optimisation, finite fields and more. Rules call them to produce numeric witnesses, and users can call them directly. |
| `api` | `Session`, `Term<'s>` (a `Copy` handle with overloaded operators), `Config` (rule sets, target phase, bindings, assumptions, budget), `Session::compute -> Answer`, and `ComputeError`. |
| `io` | Text parsing and printing (round-trippable), LaTeX, Typst, Unicode pretty-printing, markup helpers, plotting, `.npy` (feature `npy`). |
| `backend` | `Backend` trait: it compiles a closed-form term into a callable function. The reference `Interpreter` ships with the crate. JIT backends live outside the core crate (see below). |
| `sim` | Time stepping, fields and particles: FDM/FEM/FVM/BEM, CFD, MD (Verlet, BAOAB Langevin, Nosé–Hoover, cell lists), SPH, spectral, multigrid; P1/P2 finite elements on unstructured triangle meshes with adaptive refinement (`fem_triangles`); compressible gas dynamics with HLLC/MUSCL and an exact Riemann solver (`gas_dynamics`); symplectic and stiff integrators (`integrators`); Ising Metropolis and Wolff clusters; and the scenario models. `sim::scenario::run(name, json)` is the uniform entry point. |
| `ffi` | A minimal C ABI over `api` and `sim::scenario`. cbindgen generates `rssn.h`/`rssn.hpp` from it (`DEV=1 cargo build`). |

### Rule sets

`rules::standard()` installs all of these:

- `arith`, `elementary`, `calculus`, `poly`, `solve`, `ode`, `linalg`, `geometry`, `complex`, `number_theory`, `combinatorics`, `logic`, `special`, `stats`, `transforms`, `variational`, `pde`, `functional`, `physics`, `physics_quantum`, `optimize`, `verify`, `discrete` (graphs, finite and permutation groups, character tables, point groups, topology and homology, coding, cryptography, fractals), `rewriting` (Knuth–Bendix), `lie` and `lie_structure` (Lie algebras, root systems, representations), `geometric_algebra`, `units`, `radicals`.

A rule set declares the sets it depends on. Its contents are:

- **operators**, registered by name, so a plugin can add its own;
- **rewrites**: text patterns with guards such as `is_const(?a, ?x)` or `positive(?a)`;
- **definitions**: `name(a, b) := body`, expanded by a kernel;
- **kernels**: Rust functions over e-classes that return `Equal(term)`, `Pinned(term)`, `Approx(ball)` or `Pass`.

Heavy algorithms are kernels. Examples are Risch–Norman integration and the substitution/reduction-formula integrator, Gosper summation and recurrence guessing, univariate (Zassenhaus) and multivariate (Kronecker) factorisation, Gröbner bases, CAD in any dimension, the equation and inequality solver with its heuristic sub-solvers, the ODE solver (classification, variation of parameters of any order, Lie point symmetries, reductions, Fourier and Laplace methods), the PDE solver (characteristics in any dimension, separation of variables, Green's functions, conservation laws and symmetries), and Schreier–Sims / Todd–Coxeter for groups. Each kernel writes its result back as an equality.

### Phases

`Target::Symbolic` extracts with `ClosedForm`: no heavy operator may remain. When no closed form exists, the cheapest term by size is returned with `reduced = false`.

`Target::Numeric { tolerance }` evaluates the best term instead. If that fails, it uses a numeric witness (a ball) that a kernel produced, for example a quadrature result. Numbers are never mixed into the e-graph as if they were exact.

### Sessions and isolation

A `Session` owns a plain term store. Each `compute` call works as follows:

1. Clone the template graph for the requested rule configuration.
2. Import the request into the clone.
3. Run the engine on the clone.
4. Import only the answer back into the store.

Requests therefore never interfere with each other. A session is single-threaded; use one per thread.

## C interface

The C interface lives in `src/ffi`, under 1000 lines. Terms are `uint32_t` handles within a session.

| Area | Functions |
|---|---|
| Sessions | `rssn_session_new`, `rssn_session_free` |
| Building terms | `rssn_term_sym`, `rssn_term_int`, `rssn_term_rational`, `rssn_term_float`, `rssn_term_parse`, `rssn_term_apply(op_name, args, n)` |
| Configuration | `rssn_config_new`, `rssn_config_free`, `rssn_config_symbolic`, `rssn_config_numeric`, `rssn_config_bind`, `rssn_config_assume`, `rssn_config_budget` |
| Computing | `rssn_compute` (fills `RssnAnswer`), `rssn_simplify` |
| Reading results | `rssn_term_to_f64`, `rssn_term_eval`, `rssn_term_tensor_data` / `rssn_tensor_free` |
| Output | `rssn_term_to_string`, `rssn_term_to_latex`, `rssn_term_to_typst`, `rssn_term_to_pretty`, `rssn_string_free` |
| Simulations | `rssn_sim_run(name, json, &out)` |
| Errors and version | `RssnStatus` codes, `rssn_last_error`, `rssn_version` |

No call unwinds: a panic becomes `RSSN_STATUS_PANIC`. See `examples/c/smoke.c` for a complete client.

## Extending rssn

- **New operator or domain.** Write a `RuleSet` containing operators, rewrites, definitions and kernels. Add it with `Config::with`, or with `Session::with_rules` for a session. The soundness test checks your rewrites numerically.
- **New integral forms.** Register extra antiderivative functions in the `rules::calculus::IntegralTable` attribute of the `integral` operator.
- **New backend.** Implement `backend::Backend` and pass it to `Term::compile_with`, or install it with `backend::set_default`.
- **Faster special functions.** Attach a `backend::Lowering` attribute to an operator so compiled code calls it directly.

## Compiled evaluation and the JIT

Every compiled evaluation goes through `backend::current()`: the reference
`Interpreter` by default, replaced globally with `backend::set_default` or
for one scope with `backend::with_backend`.

- **Tape IR** (`backend::tape`). Terms are lowered to an SSA tape that undoes
  the canonical form (subtraction, division, integer powers by
  square-and-multiply, `sqrt`/`cbrt`), folds constants, removes dead code and
  is hashed structurally. The same tape runs over `f64`, dual numbers
  (forward-mode gradients) and outward-rounded intervals (rigorous
  enclosures, which also decide signs of constant classes in `graph::facts`).
- **Cranelift JIT** (`backend::jit`, feature `jit`). The tape is compiled to
  machine code with the signature `fn(inputs, params, out)`. Compiled
  functions are cached by structural hash; `TieredBackend` interprets first
  and switches to machine code compiled in the background. Operators can
  declare a `backend::Lowering` (an intrinsic or an `extern "C"` function);
  the special functions do, others go through a generic evaluation shim.
  All `unsafe` code is confined to the `abi` module.
- **Users.** Numeric quadrature, infinite sums, `odeint` (with a compiled
  symbolic Jacobian for the stiff Radau IIA fallback), `nminimize` (compiled
  objective and gradient), `chebyshev_approx`/`rational_approx`, the numeric
  solver and `Term::compile` all go through the current backend, so turning
  on the JIT speeds up every one of them.

## Lint policy

`src/lib.rs` applies the strict configuration: `clippy::all`, `pedantic`,
`nursery`, `unwrap_used`, `expect_used`, `indexing_slicing`,
`arithmetic_side_effects`, `single_call_fn` and `missing_docs` are denied
(the later `allow(clippy::restriction)` group keeps the individual
restriction lints at their group level); `dead_code`, `unsafe_code` and
`warnings` warn. A few lints are deliberately allowed for scientific code,
each with its reason next to it in `lib.rs`:

| allowed | why |
|---|---|
| `suboptimal_flops` | `mul_add` is a slow libm call without hardware FMA, and fused rounding makes results differ between targets |
| `cast_precision_loss` | counts and indices become `f64` throughout numerics and stay far below 2^53 |
| `many_single_char_names`, `similar_names` | formulas keep the notation of the literature |
| `tuple_array_conversions` | false positives on arrays built from separately computed bindings |
| `cast_possible_truncation`, `cast_sign_loss`, `cast_possible_wrap` | from the base configuration |

Local `#[allow]`s are used sparingly and carry a comment.

## Testing

- **Unit tests** sit next to each rule set and kernel. Each rule set's tests check identities symbolically and, where applicable, numerically.
- **`graph::soundness`** evaluates both sides of every rewrite of every standard rule set at sampled points.
- **`graph::proptests`** runs property tests of the kernel.
- **Integration tests** live in `tests/` (io, kernels, sim).
- **Benchmarks:** `cargo bench --bench engine`.
