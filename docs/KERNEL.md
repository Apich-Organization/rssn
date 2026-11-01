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

### 3. No silent unknowns

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

`cargo run --bin signature-audit` reports the split and exits non-zero if any
operator is neither typed nor explained. Today: 1369 operators, 799 typed
(58.4%), 570 explained, 0 unexplained.

### 4. Determinism

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
let value = session.call_value("sin", &[Value::Real(0.5.into())], &CallCtx::new())?;

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
mean "all one" as a claim.

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

`DimExpr` is the *type-level* mirror of a dimension: names instead of arena-local
ids, so it can go in a manifest and be read back by another process.
`DimExpr::to_term` interns unknown names on demand, which is how a host defines
a base dimension (`Currency`) the engine has never seen.

## Threads and handles

`Session` is `Send + Sync` and is the unit of concurrency. `Term<'s>` borrows a
session and is not `Send`; `TermRef` is `Copy + Send + Sync` and owns nothing —
it is a `(session_id, node)` pair, refused by any session but the one that minted
it. That turns a use-after-free of a term handle into an error instead of a wrong
answer.

`fork` gives a session its own graph; `merge` imports every node of another and
translates dimension assumptions between the two arenas, failing on a
disagreement rather than picking a winner.

## The C ABI

```c
const char *json = rssn_manifest_json();   /* caller frees */
uint64_t    hash = rssn_manifest_hash();
```

Both describe the *standard* build. A host that installed plugins reads the
richer per-session manifest through the Rust API.

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
