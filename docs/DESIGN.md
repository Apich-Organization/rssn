# Design notes for the open items

This document records the **architecture decisions** behind the remaining gaps in
`RSSN-INTEGRATION-REQUIREMENTS.md`. Each section states the problem, the options
that were considered, the decision, and what follows from it. Sections marked
*implemented* have code and tests; sections marked *designed* have a decision and
a scope but no code yet, and the reason is stated.

---

## 1. Overload dispatch: one decision, three entry points

### The problem

`OpTable` allows several registrations to share a name, and each registration may
declare several `Signature`s. Three different things were being called "the
overload", and they disagreed:

* `call_value` picked a **registration** by arity, then executed it without
  looking at the chosen overload's parameter types. The two overloads inside one
  registration (numeric `(Real) -> Real`, symbolic `(Term) -> Term`) were never
  selected — the operator's `eval` decided at run time whether to fold.
* `resolve_overload` could return an `OpKey` with an all-zero `signature_hash`
  for an operator that *has* signatures, and `call_key` rejects that key. A
  caller could therefore be handed a key it cannot use.
* `Kernel::call(&OpKey, …)` ignored the signature hash entirely and dispatched
  by name.

### Decision

**One resolution function, and every entry point goes through it.**

```
                    ┌──────────────────────────────┐
   name + args ────►│ resolve: arity, then type    │───► (OpId, overload index)
                    │ most specific; tie = error   │
                    └──────────────────────────────┘
                                   │
        ┌──────────────────────────┼──────────────────────────┐
        ▼                          ▼                          ▼
  call_value(name, …)      resolve_overload(name, …)     call_key(key, …)
  resolves, then calls     resolves, returns a key      *skips* resolution
```

* `resolve_overload` returns `OpKey` **only** when it has a definite answer. If
  the call is ambiguous it returns `Err(OpError::Ambiguous)`; it never returns a
  key that `call_key` would reject. The key it returns is always callable, which
  is the property a caching caller needs.
* `call_value` validates the arguments against the selected overload and returns
  `OpError::Type { expected, got, arg }` when one does not match. Executing a
  `(Real) -> Real` overload on a tensor is not a "conservative" behaviour, it is
  a wrong answer waiting to happen.
* `call_key` does *not* re-resolve: the whole point of the key is that the caller
  decided already. It looks the registration up by `(rule_set, name,
  signature_hash)` and calls exactly that overload.
* `Kernel::call` is implemented in terms of `call_key`, so the trait cannot
  silently downgrade to name dispatch.

An all-zero `signature_hash` keeps one meaning only: *the registration with no
signature at all* (a leaf). A signed operator's keys are never zero.

### Consequence

A host that resolves once and caches the key gets the same operator on every
later call. A host that calls by name gets the same operator the resolver would
have chosen, and an error instead of a wrong overload.

---

## 2. Signature precision: three layers, not one

### The problem

100% of operators are typed, but only 20 declared + 13 overlay are *semantically*
precise. The other 1349 are `Derived`, which distinguishes `(Real, …) -> Real`
from `(Term, …) -> Term` and nothing else: `det` is `(Term) -> Term` even though
it is a matrix determinant.

Hand-writing 1349 signatures is not the answer — it would be 1349 judgement calls
that no test can check, and it would rot the moment an operator changes.

### Decision

Three layers, each with a different trust story:

| Layer | Written by | Precision | Verification |
|---|---|---|---|
| `Declared` | the rule set, next to the operator | full | the operator's own tests |
| `Overlay` | a human, in `signatures.overlay.json` | full | the overlay's own test |
| `Derived` | the rules in §2.1 | shape-level | unit tests on the rules |

**The overlay is the mechanism for precision, and it is meant to grow.** A host
that needs `det: Matrix -> Real` gets it from the overlay; the crate ships the
thirteen that matter most today and the file is designed to be appended to.

§2.1 also gained one rule that raises precision for free: **a `Binder` makes its
argument a `Symbol`**. That alone fixes every `diff`/`solve`/`sum`/`product`/
`limit` signature, because the binder is where the operator declares which
argument is the variable.

### 2.1 What `Derived` may and may not claim

It may claim what the registration states:

* `eval` present → the numeric contract exists;
* `PREDICATE` → the result is a boolean;
* `Binder` at *k* → argument *k* is a symbol;
* `Variadic` → the tail repeats;
* `LEAF` → no call signature.

It may **not** claim a tensor rank, a dimension, or an element type, because the
registration does not state one. A host must therefore treat `Derived` as
*conservative*, which is exactly the label the manifest carries.

### Considered and rejected: deriving dimensions from `dim_rule`

`OpDescriptor.dim_rule` already says how an operator combines dimensions
(`Same(0)`, `Sum`, `Product`, `Quotient`, `Root(n)`, …), and 275 operators carry
one that is expressible as a `DimExpr`. It is tempting to turn that into a
`Derived` signature like `(Real[d]) -> Real[d]`.

**It should not be done, and here is why.** A signature is a promise the runtime
keeps. The overload matcher scores *arguments*, so a `ret` is never checked: a
call could pass dimensions the signature forbids and get a result with a
different one, and nothing would object. Emitting `(Real[d]) -> Real[d]` would
therefore be a claim the engine does not enforce — strictly worse than the
current silence, because a host would generate generic bindings on the strength
of it.

Enforcing it is a different feature — a post-call check of the result's dimension
against the instantiated `ret` — and that is a change to call semantics with its
own false-positive risk (an operator that legitimately returns a different
dimension because a rule fired).

The information a host actually needs is **already in the manifest**, in the
right place: `OpInfo.dim_rule` states how the operator combines dimensions, per
operator, without pretending the type system checks it. A host that wants
dimension-aware bindings reads `dim_rule`; it does not need it duplicated into a
signature.

---

## 3. Session and graph: what is shared, what is copied

### The problem

`Session` holds a template graph; `compute` clones it per request and installs
the request's terms into the clone. Two complaints follow: the clone costs
proportional to the graph, and `Session` is single-threaded.

### Decision

**The clone is the design, and the shared part is made as large as possible.**

* `Graph.ops` is `Arc<OpTable>` with copy-on-write (`ops_mut` clones only when a
  registration happens, which is never during a request). The operator table —
  the largest constant part — is therefore *shared*, not copied.
* `Graph.dims` is `Arc<RwLock<DimArena>>`, **shared until written**. The first
  mutation on a cloned graph detaches it (the manual equivalent of
  `Arc::make_mut`, which cannot be used because `RwLock<T>` is not `Clone`), so
  a base dimension defined during one request cannot reach the template — and a
  request that defines nothing shares the arena and copies nothing. This is
  what makes item 13 of the integrator's list true *and* cheap.
* The node arena and the interner *are* copied. That is the remaining cost.

**Threads**: `Session` is single-threaded on purpose. `fork` gives each thread its
own graph, `merge` folds the results back. The alternative — a lock around every
rewrite — would serialise exactly the work the engine parallelises internally.

### Designed, not implemented: a reusable scratch graph

The next step is a per-session scratch graph that is *reset* rather than
re-cloned: keep the node arena's capacity, clear the nodes, re-import the
template's terms. That removes the allocation churn without changing the
isolation story.

It is not implemented because `Term<'s>` borrows the session and the engine holds
`&Graph` across a request; reusing the storage means the scratch graph must
outlive the borrow, which forces `Term` to carry a generation counter and the
engine to take a guard. That is a lifetime redesign, not a refactor, and it must
not be done in the same round as the correctness fixes above. Scope: a
`Scratch` type in `api`, a generation counter on `Term`, and a benchmark that
shows the clone is actually the cost.

---

## 4. ODE numerics

### What exists

Four integrators with continuous extensions:

| Method | Dense output | Events |
|---|---|---|
| `dopri5` | `DenseKind::Poly` | yes |
| `rosenbrock23` | `DenseKind::Hermite` | no (no event loop) |
| `radau5` / `radau_with` | `DenseKind::Lagrange` | **plumbed but not exposed** |
| `bdf` / `bdf_with` | `DenseKind::Backward` | **plumbed but not exposed** |

Event location itself is method-independent: `process_events` evaluates the event
function at both ends of an accepted step, and Brent-iterates on the dense
extension when the sign changes. Radau and BDF already call it; their
`RadauOptions.events` / `BdfOptions.events` fields are simply not reachable from
the convenience entry points.

### Decision (implemented)

**Event location is a property of the dense output, so every integrator with a
dense output gets events, and the entry points say so.**

* `radau5`, `radau_with`, `bdf`, `bdf_with` gain `*_events` variants; the
  existing names delegate with an empty slice, so nothing breaks.
* `rosenbrock23` gains an event loop using its Hermite extension — it has the
  extension already, it was just not wired.
* `Rk4` is fixed-step and has no extension. Rather than refuse events, it
  **locates them by step bisection**: when the event function changes sign
  across a step, the step is re-integrated with halved `h` until the bracket is
  below a tolerance, then Brent-refined on a cubic Hermite interpolant built
  from the two endpoint derivatives. That is a real event locator, not an
  interpolation of convenience, and it is the standard treatment for
  fixed-step methods.

### Implemented: fixed-step event location

`solve_fixed_events` bisects the *step length* and re-integrates the partial
step from the same start. Because Euler, Heun and RK4 are one-step methods, the
state at any intermediate time is just another application of the same step, so
the located state is a point of the method's own trajectory rather than of an
interpolant — and no dense output is needed at all. The cost is ~50 extra
evaluations per event, which is the right trade for a method that is cheap by
construction. `Rk4`'s events are therefore as accurate as its trajectory, not
limited by a third-order Hermite fit.

### Designed, not implemented

* **BDF for index-2/3 DAEs.** The BDF corrector is already a Newton solve on a
  mass-matrix system; higher index needs the Hessenberg structure (`n2`, `n3`
  blocks, as in `tests/kernels/andrews_squeezer.rs`) and a different error
  estimator. This is a new method, not a completion: the order conditions, the
  error estimate and the step-size controller all change, and there is no
  existing benchmark in the tree to validate it against. It stays scoped here so
  it is not forgotten, and it is deliberately *not* attempted in the same round
  as the correctness fixes.
* **A variable mass matrix.** The BDF Jacobian assumes a constant `M`. Allowing
  `M(t)` needs `dM/dt` in the Newton matrix — a term the current code has no
  place for — and a different order condition, because the mass matrix enters
  the backward-difference identity. Same category: a new method, not a gap.

---

## 5. Complex numbers: one representation

### The problem

`Value::Complex` has `Exact { re, im }` and `Float { re, im }`, but `Number` has
no complex variant: a symbolic complex is an expression with `i` in it, and a
numeric one is `Complex64`. Two paths, and a host has to know which it is on.

### Decision

**`Number` stays real; `Value::Complex` is the boundary type, and the graph's
`i` is the symbolic form.**

The reasons are concrete:

* The e-graph's hash-consing, extraction cost and rewrite rules are all built on
  an ordered, exactly comparable `Number`. Adding `Complex` to it means a new
  total order on `Complex<BigRational>` (which has none that respects the field
  operations) and revisiting every rule that compares numbers.
* A symbolic complex number is already representable: `3 + 4*i` is a term.
  `Value::Complex` exists so that a *numeric* complex does not have to be
  smuggled through the term graph.
* The bridge is explicit: `Value::Complex` with a zero imaginary part converts to
  a real term; anything else stays at the boundary and reaches an operator
  through `EvalDyn`.

A host should therefore keep one `Complex` type on its side and convert at the
boundary, which is what `Value::Complex` is for.

---

## 6. Constants are values, not units

### Decision

**Constants live in their own manifest section and never enter `UnitRegistry`.**

A constant is a number with a dimension; a unit is a scale with a dimension. They
share the shape (`Scale` + `DimExpr`) but not the role: `c` is not a unit you can
write `quantity(3, c)` with, and `h` (the hour) is not a constant. Merging them
would make `quantity(2, c)` parse and mean something absurd.

The name clash the integrator flagged — `e` (elementary charge) versus `Math.e` —
is a consequence of putting both in one namespace, which is exactly what is
avoided here. A host that opens `Units` and `Constants` into one module must
qualify `e`; the manifest makes that visible because `constants` and `units` are
separate lists.

---

## 7. `manifest_hash` presentation

### The problem

`signature-audit` printed `manifest ce50f1e0fd13169f`, which is
`hash_u64()` — the first eight bytes of the 32-byte BLAKE3 digest read as a
**little-endian** `u64`. A reader naturally assumes it is the first eight bytes in
hex order, and the two differ.

### Decision

**The full 32-byte digest is the identity; `hash_u64` is a convenience and is
labelled as one.**

* `signature-audit` prints the 64-character lowercase hex of the full digest, and
  separately the `u64` with the word `u64` next to it.
* `KernelManifest::hash_hex()` is the canonical string form.
* JSON is never hashed; it is a human view. This is stated in the manifest's own
  documentation.

---

## 8. JIT policy

`Select`, `Extern1`, `Extern2` and `Gather` are lowered, and the cache key is
`mix(tape_hash, manifest_hash ^ target_key())`. The backend is behind the `jit`
feature because Cranelift is a large dependency and a build that never compiles a
kernel should not pay for it.

`Features.jit`, `jit_branch` and `jit_index` report the build, so a host that
wants `@jit` reads them and refuses with a precise capability error when they are
false. Making `jit` a default feature is a distribution decision, not an
architectural one; the manifest already carries the answer either way.

---

## 9. What the simulation schema needs

`SimResult` is a `Value::Record` with a `FlattenSchema`, so the *result* side is
structured. The *input* side is still a JSON scenario. The design is to give a
scenario the same treatment: a `ScenarioSchema` describing its parameters, with
the same `FlattenSchema` machinery, so that a host generates `Physics.*` from the
manifest instead of hand-writing it.

Not implemented: it needs a parameter-description type (name, kind, unit, range,
default) that the scenario loader and the manifest share, and the scenario loader
currently parses JSON directly into a struct. Scope: `sim::scenario::Schema`,
`KernelManifest.scenarios`, and a round-trip test.
