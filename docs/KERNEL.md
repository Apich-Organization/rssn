# The embedding kernel

This document describes `rssn::kernel` and `rssn::dim`: the two modules a host
program talks to when it drives rssn from outside the crate. It is written for
someone implementing that host, and for someone changing the engine without
breaking them.

## Why two modules and not two crates

A downstream project asked for `rssn-api` and `rssn-dim` as separate crates. We
did not split the crate, because both halves share the expression graph, the
operator table and the dimension arena: a crate boundary would mean either
re-exporting half of `graph` through a facade or duplicating it. The *interfaces*
live in these two modules instead, with no `pub(crate)` shortcuts, so the split
remains a packaging decision rather than an architectural one.

| Module | Would have been | What it is |
|---|---|---|
| `rssn::kernel` | `rssn-api` | manifest, dynamic value ABI, signatures, plugins, codec, structured results |
| `rssn::dim` | `rssn-dim` | dimension terms, the dimension arena, the unit registry |

## The four contracts

### 1. One wire format

`kernel::codec` is the only place that turns a value into bytes, and it always
uses `bincode-next`:

* `encode_deterministic` / `decode_deterministic` — deterministic bincode, for
  hashing and for storage. The same value always produces the same bytes.
* `encode_cbor` / `decode_cbor` — deterministic CBOR, for interchange.
* `digest` / `short_digest` — BLAKE3 over the deterministic encoding.

Nothing in rssn writes JSON to a peer. `KernelManifest::to_json` exists for
humans and for the C ABI's `rssn_manifest_json`, and is documented as such.

The one JSON file the crate *reads* is `signatures.overlay.json`. It is
version-controlled source (the `*.json` ignore rule has an explicit exception for
it); `build.rs` copies it into `$OUT_DIR` and the crate `include_str!`s that copy.
It is the one place a human corrects an inferred signature. A missing file fails
the build and an invalid signature fails manifest construction; see
[BUILDING.md](BUILDING.md).

### 2. One description

An operator is described once, by `OpDescriptor`, and that description becomes:

* a `Signature` (declared, from `signatures.overlay.json`, or inferred);
* an `OpKey` — `(rule_set, name, signature_hash)` — which is stable across
  builds and registration orders;
* an `OpInfo` in the manifest.

```rust
let manifest = session.manifest();
let sin = manifest.op("sin").expect("sin is a core operator");
assert!(!sin.overloads.is_empty());
let key = &sin.keys[0];                 // cache on (manifest_hash, key)
```

### 3. Every operator is typed, and the type is derived from the registration

A host generating typed bindings cannot work from names. Every non-leaf
operator therefore has at least one signature, and it is *derived from what the
registration states*, never guessed:

| The operator declares | Therefore |
|---|---|
| a scalar evaluator (`eval`) | `n` reals in, a real out |
| `OpFlags::PREDICATE` | the result is `Bool` |
| a `Binder` at index `k` | argument `k` is the bound `Symbol` |
| `Arity::Variadic` | the last parameter repeats |
| `OpFlags::LEAF` | no call signature; a `dynamic_reason` instead |

and every non-leaf operator also has the **symbolic contract** — `n` terms in,
one term out. An operator with a scalar evaluator therefore has two overloads,
`(Real, …) -> Real` and `(Term, …) -> Term`, which is exactly what the engine
does: `sin(0.5)` folds, `sin(x)` stays a term. `SignatureSource` says where a
signature came from: `Declared` (the rule set wrote it), `Overlay` (a human
corrected it in `signatures.overlay.json`), `Derived` (the rules above),
`Witnessed` (derived, plus overloads the engine was seen to satisfy; next
section), or `None`.

### 3a. Witnessed signatures: the engine types its own operators

The derivation above is sound and nearly empty: for `det`, `factorint` or
`groebner` the registration only says "a heavy operator of arity `n`", so the
derived signature is the symbolic contract `(Term, ..) -> Term`. A host cannot
generate a useful typed binding from that, and writing a thousand signatures by
hand is the drift the manifest exists to prevent.

`kernel::sort` closes the gap with evidence the engine already produces. A
finished run leaves a request such as `factorint(360)` in the same equivalence
class as its answer `list(list(2, 3), list(3, 2), list(5, 1))`; that class is a
typing judgement. The module

1. defines a small lattice of **sorts** (`Int`, `Rat`, `Float`, `Sym`, `Expr`,
   `Bool`, `Prop`, `Obj<name>`, lists with a rectangularity flag, tuples) with a
   join, and classifies any e-class into one using only what the graph states
   (literals, list constructors, the `PREDICATE` flag, a non-scalar dimension
   rule);
2. with `RSSN_SORT_TRACE=<file>` set, logs one *fact* per operator application
   at the end of every top-level run: operator, argument sorts, result sort,
   and the call that showed it;
3. `signature-mine` treats the log as a source of *calls*, not of truths (a
   test may have loaded only some rule sets, and an application is never
   evidence about itself: `gcd(12, 18)` alone in its class says nothing about
   `gcd`). It runs every distinct call again on the full engine, keeps the
   ones that are answered, clusters them by the kinds of their arguments and
   joins each cluster (anti-unification).
4. Each candidate is then **probed**: the numbers in its witnesses' reduced
   arguments are shifted, shapes kept, and the call is run again. This is what
   stops `sin(0) = 0` from typing `sin` as `(Int) -> Int`, or one matrix of
   perfect squares from typing `matsqrt` as integer-valued; with fewer than
   three facts and a symbolic argument, exact number leaves of the result are
   widened to expressions as well.
5. A signature is kept when it says something about the *answer* and its
   parameter list is comparable with, or disjoint from, every other overload
   of the operator. One whose parameters equal an overload the registration
   already states is kept as a **refinement** of that overload's result
   (`holonomic_zero: (Term, Term) -> Bool` next to the symbolic contract),
   listed separately in the manifest and taking no part in resolution.

```text
RSSN_SORT_TRACE=sorts.tsv cargo test --lib       # the test suite is the corpus
cargo run --profile tools --bin signature-mine -- sorts.tsv
cargo run --profile tools --bin signature-mine -- --check   # replay the shipped table
```

The table is versioned like the overlay (the manifest of a clean clone must not
depend on having run the tests) and regenerated, never edited. In the manifest
such an operator has `signature_source: witnessed`; `overload_sources` says per
overload whether it is derived or witnessed, `refinements` lists the sharper
results of existing overloads, and `witnesses` carries the calls. At
3.0.0-rc.1 the table types the answers of about 700 operators that were
`(Term, ..) -> Term` before.
Every operator that has a test is thereby typed, and a new operator arrives
typed as soon as it is tested.

What a witnessed signature promises is what a gradual type promises at a
dynamic boundary. It is evidence, not a totality proof: the parameter side may
be narrower than what the operator accepts (the symbolic contract stays as the
most general overload, so nothing becomes uncallable), and a host must validate
the result it decodes. Overload sets stay unambiguous by construction: a
witnessed signature is kept only when its parameter list is comparable with, or
disjoint from, every other overload (`sort::coherent`), and the overload
matcher breaks a tie between comparable overloads by pointwise dominance.

### 4. No silent unknowns

Every operator is either typed or explained. The explanation is a
`DynamicReason`:

| Reason | Meaning |
|---|---|
| `request` | reduced by rules; the result type is whatever the rule that fires produces |
| `variadic` | the arity is not fixed |
| `fully_polymorphic` | it accepts any type |
| `shape_dependent` | a container whose shape depends on the arguments |
| `opaque_payload` | it carries host data the kernel does not interpret |
| `not_yet_specified` | a declaration is missing; the audit treats this as a defect |

`cargo run --bin signature-audit` reports the split between declared, overlay
and derived signatures and exits non-zero if any operator is neither typed nor
explained. At 3.0.0-rc.1 it reports 1476 operators in 37 rule sets, every one
typed (124 declared, 25 overlay and 1352 derived signatures), 82 units, 30
constants (18 exact) and 0 unexplained operators. The counts change whenever a
rule set changes; the audit (and `KernelManifest`) is the authority, not this
paragraph.

### 5. Determinism

A manifest built from the same operator table and the same rule sets is
byte-identical, in any process, in any registration order. The mechanism:

* every list in the manifest is sorted and deduplicated;
* no ids, addresses, timestamps or hash-map iteration orders appear;
* `manifest_hash` is BLAKE3 over the deterministic encoding of the manifest
  with the hash field zeroed.

A host keys caches on `(manifest_hash, OpKey, target)`.

## Calling an operator dynamically

```rust
use rssn::kernel::{CallCtx, OpError, Value};

let session = rssn::api::Session::new();

// By name, with values.
let value = session.call_value("sin", &[Value::Float(0.5)], &CallCtx::new())?;

// An answer that does not reduce comes back as a handle to a term.
let x = Value::Symbol("x".to_owned());
let derivative = session.call_value("diff", &[x.clone(), x], &CallCtx::new())?;
```

`CallCtx` carries the target (symbolic or numeric), the budget, the dimension
mode, assumptions, symbol bindings, substitutions, host callbacks, and an
optional trace. It is read-only for the operator being called; only the trace
is appended to.

Errors are structured, not strings:

```rust
match session.call_value("sin", &[], &CallCtx::new()) {
    Err(OpError::Arity { expected, got }) => { /* ... */ },
    Err(OpError::UnknownOp(key)) => { /* ... */ },
    _ => unreachable!(),
}
```

### Choosing between overloads

`call_value` resolves by arity and then by type, most specific first; a tie is
`OpError::Ambiguous` rather than a silent choice. A generated facade should
resolve once and then call through the key it got, so the static choice and the
executed operator cannot disagree:

```rust
let args = [Value::Float(0.5)];
let key = session.resolve_overload("sin", &args)?;   // OpKey
session.call_key(&key, &args, &CallCtx::new())?;     // exact dispatch
```

`Session::overloads(name)` lists every signature with its key. An operator with
no signature has no overloads and is reachable by name only; its key carries an
all-zero `signature_hash`.

The match lattice is `Exact` > `Widened` (an `Int` into a `Real`, a `Real` into
a `Complex`) > `Shaped` (same constructor, looser argument) > `Generic` (a type
variable or `Dyn`), and a call is scored by its *worst* parameter.

## Plugins

A plugin is a bundle of extensions registered through `Registry`:

```rust
use rssn::kernel::{Registry, SimplePlugin};

let plugin = SimplePlugin::new("acme/twice", "0.1.0", |reg: &mut Registry<'_>| {
    reg.register_op(
        rssn::graph::OpDescriptor::new("acme_twice", rssn::graph::Arity::Fixed(1))
            .signature(rssn::kernel::Signature::scalars(1, rssn::kernel::TypeExpr::Real)),
    )?;
    reg.register_rewrites(
        rssn::graph::Tier::Normalize,
        &["acme/twice-zero: acme_twice(0) => 0"],
    )?;
    Ok(())
});
session.register(Box::new(plugin))?;
```

Registration refuses an operator that declares neither a signature nor a
`DynamicReason`: a plugin cannot make the manifest incomplete. It also refuses
a name that clashes with an existing operator of a different shape; adding an
overload of the same name requires `.overload()` on the descriptor.

## Physical constants

`dim::constants` carries the CODATA 2022 constants with an exact rational value,
a dimension and an SI unit, and distinguishes the ones that are **exact by
definition** (the seven SI defining constants and everything derived from them)
from the **measured** ones, which carry a standard uncertainty. `c` really is
the integer `299792458`; `hbar` is `h/(2π)` exactly, carried as a rational times
`π^-1`.

`KernelManifest.constants` lists them, so a host can generate
`Units.Constants` without a second table. Note that `G` is the gravitational
constant here and the gauss in the unit table — different namespaces, on
purpose.

## Loading a plugin from a shared library

`kernel::cdylib` defines a byte-level C ABI: one `bincode-next` buffer of
`Vec<Value>` in, one of `Result<Value, String>` out. A `#[repr(C)]` mirror of
`Value` would have to track every variant across two independently compiled
artifacts and would break silently the day one is added; the deterministic
encoding is already the stable description of a value.

It is behind the `plugin-cdylib` feature because a shared library runs
arbitrary code in the host process. `Features::plugin_kinds` reports whether the
build can load one, so a host that meets a `kind = "cdylib"` plugin in a build
without the capability refuses it with a precise diagnostic.
`tests/cdylib_plugin.rs` builds a plugin in C with `cc` and calls it end to end.

## Build capabilities

`KernelManifest.features` is what a host gates on before generating bindings:

| Field | Meaning |
|---|---|
| `dynamic_call`, `parallel_symbolic`, `higher_order` | the dynamic value ABI, the parallel symbolic phase, host callbacks |
| `jit`, `jit_branch`, `jit_index` | the optional JIT and what it lowers |
| `exact_units`, `affine_units`, `user_base_dimensions` | the exact unit layer |
| `structured_results` | `Solutions` / `OdeSolution` records cross the ABI |
| `sde` | the `sde` rule set ships, so the noise operators below are callable |
| `plugin_abi`, `plugin_kinds` | the plugin ABI, and how a plugin may reach this build |

`sde` is specifically about *noise*: `sde_wiener(seed, m, steps, h)` returns the
`steps × m` matrix of seeded Brownian increments (dimension `sqrt(T)` for a step
of dimension `T`) and `sde_levy_areas(seed, h, dw, terms)` the iterated Ito
integrals conditioned on given increments. `sde_scheme_info(name)` answers what a
scheme promises — `list(strong_order, ito, commutative_only, adaptive,
levy_areas)` for the eight solvers — without running it; the name crosses as a
string value. The solvers themselves are not in the manifest and never will be:
they take host closures for `a(t, x)` and `g(t, x)`, which a `Value` cannot
carry. A host that wants an ensemble keeps its own loop and uses the noise
operators so that every scheme in the loop sees one path.

## Which values can enter the graph

`KernelManifest.value_kinds` lists every `Value` variant with
`in_term_graph: bool`. A variant with `in_term_graph = true` becomes a node of
the expression graph and can be passed to any operator; one with `false` reaches
only an operator with a dynamic evaluator, or comes back as a result.

`Int`, `Real`, `Float`, `Bool`, `String`, `Symbol`, `Unit`, `Quantity`,
`Tensor`, `Term` and `Dyn` enter the graph. `Complex` enters only when its
imaginary part is zero. `Sparse`, `Record`, `Variant`, `Function` and `Opaque`
do not, and the manifest says why in each `note`.

## Structured results

A result that is more than one number is a `Value::Record` with a declared
`FlattenSchema`. That is deliberate: it crosses the dynamic ABI unchanged, it
serialises with `bincode-next` for free, and the layout travels with the value.

| Type | Fields |
|---|---|
| `Solutions` | `variable`, `roots`, `multiplicities`, `conditions`, `complete` |
| `OdeSolutionView` | `t`, `y`, `flatten_order`, `events`, `nfev`, `rejected` |
| `SimResult` | `scenario`, `t`, `fields`, `schema` |

`Session::solve_structured` lifts the `solve` operator's `list(...)` into a
`Solutions`, and `OdeSolution::view` names an integrator's anonymous columns.

An empty `multiplicities` means the solver did not distinguish them; it does not
mean "all one" as a claim. `complete` is true only when the solver answered in
its exhaustive `list(...)` form; a single root found on the way is reported with
`complete = false`, because it is not a claim that the others do not exist.

`OdeSolutionView::dimensions` travels inside the record, so a host reading the
`Value` does not have to go back to Rust for it.

## Dimensions and units

`rssn::dim` owns the only unit table in the crate.

* `UnitDef` — name, aliases, `DimExpr` dimension, exact `Scale` (`rat * pi^pi`),
  optional affine offset, prefixability, display, category.
* `UnitRegistry` — `list`, `lookup`, `resolve`, `register`, `format`, `snapshot`.
  `StandardRegistry` is the SI table; `register` adds to a process-wide overlay.
* Registration is checked: a name or alias that already belongs to a different
  unit is refused, so a plugin cannot silently redefine the metre.
* Everything is exact. The scale of a foot is the rational `381/1250`, the scale
  of a degree is `pi/180`, the offset of `degC` is `27315/100`. A float would
  lose the distinction between `1 in` and `0.0254 m` after a few conversions.
* The table covers SI, the non-SI units accepted with it, the IAU astronomy
  units (`au`, `ly`, `pc`, `Rsun`, `Msun`, `Mearth`, `jup`), the CGS set
  (`dyn`, `erg`, `Gal`, `P`, `St`, `G`, `Oe`, `Mx`, `Ba`), and the common
  non-SI units of the Aevia design standard (`Da`/`u`, `kn`, `mph`, `acre`,
  and the dimensionless information units `bit`, `B`, `KiB`) — 82 units at
  3.0.0-rc.1.
* `KernelManifest.units` carries `aliases`, `prefixable`, `category` and
  `since` as well as the name, dimension and scale, so a host can generate its
  unit table from the manifest alone.

A typed front end should declare dimensions with `Config::dimension_term`,
which takes a `DimExpr` and never parses text:

```rust
let config = Config::new().dimension_term("x", DimExpr::base("Currency"));
```

Unlike `Config::dimension(symbol, "Currency")` it cannot confuse a unit symbol
with a base dimension, and a base the arena has never seen is defined on the
way in.

`DimExpr` is the *type-level* mirror of a dimension: names instead of arena-local
ids, so it can go in a manifest and be read back by another process.
`DimExpr::to_term` interns unknown names on demand, which is how a host defines
a base dimension (`Currency`) the engine has never seen.

## Threads and handles

A `Session` is single-threaded (its store sits behind a `RefCell`, so it is not
`Sync`): use one session per thread. `Term<'s>` borrows a session and cannot
outlive it; `TermRef` is `Copy + Send + Sync` and owns nothing — it is a
`(session_id, node)` pair, refused by any session but the one that minted it.
That turns a use-after-free of a term handle into an error instead of a wrong
answer.

`fork` gives a session its own graph, which can be moved to another thread;
`merge` imports every node of another and translates dimension assumptions
between the two arenas, failing on a disagreement rather than picking a winner.
`Session::compute_all` is the supported way to compute several requests in
parallel inside one session.

## The C ABI

```c
const char *json = rssn_manifest_json();   /* caller frees */
uint64_t    hash = rssn_manifest_hash();
```

Both describe the *standard* build. A host that installed plugins reads the
richer per-session manifest through the Rust API. The rest of the C interface is
summarised in [ARCHITECTURE.md](../ARCHITECTURE.md#c-interface).

## What the JIT can do

* Arithmetic, `sqrt`, `abs`, `floor`, `ceil`, comparisons, `select` and integer
  powers are native instructions.
* `Inst::Gather` is the load/index instruction: `inputs[base + offset]`, with
  the index clamped before the address is computed, so an out-of-range index
  reads NaN and never memory outside the buffer.
* An operator can declare how it compiles with `OpDescriptor::jit(Lowering)`:
  `Intrinsic`, `Extern1`, `Extern2` or `Gather`. Without a declaration, the
  backend decides: a scalar evaluator is called through a shim, everything else
  runs in the interpreter.
* The cache key is `mix(tape_structural_hash, manifest_hash ^ target_key())`,
  so two builds for different micro-architectures, or with different operator
  semantics, never share machine code.

## Adding an operator

1. Register it with `Installer::op` or `Registry::register_op`.
2. Declare a `Signature` (`.signature(..)`), or a `DynamicReason`
   (`.dynamic_reason(..)`).
3. Add `.doc(..)`, `.category(..)` and `.since(..)` — the manifest is the
   documentation a host reads.
4. Run `cargo run --bin signature-audit` and `cargo test`.
