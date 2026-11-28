# Checked host callbacks

[中文](CALLBACKS.zh.md)

RSSN's dynamic ABI supports serial, synchronous host callbacks. Register a
`HostCallback` in a `CallCtx`, then pass its returned handle to `Session::call_ho`.
The handles are appended after the ordinary arguments in their supplied order;
they are ordinary function parameters in the operator's declared signature.
For callback slots before or between ordinary arguments, put
`Value::Function(handle)` at each position declared in the manifest and use
`Session::call_value` or `Session::call_key` with the same context.

```rust
use std::sync::Arc;
use rssn::api::Session;
use rssn::kernel::{CallCtx, OpError, Value};
use rssn::kernel::plugin::NativeCallback;
use rssn::kernel::signature::{Signature, TypeExpr};

let session = Session::new();
let mut context = CallCtx::new();
let square = context.add_callback(Arc::new(NativeCallback::new(
    |args: &[Value]| {
        let x = args.first().and_then(Value::as_f64)
            .ok_or_else(|| OpError::Host("expected a floating-point input".into()))?;
        Ok(Value::Float(x * x))
    },
    Signature::new(vec![TypeExpr::Float], TypeExpr::Float),
)))?;
let result = session.call_ho(
    "num_integrate",
    &[
        Value::Float(0.0), Value::Float(1.0),
        Value::Float(1e-10), Value::Float(1e-10), Value::int(128),
    ],
    &[square], &context,
)?;
# Ok::<(), OpError>(())
```

The result is a record with `value`, `error`, `evaluations` and `converged`
fields. `error` is the adaptive Gauss–Kronrod error estimate, not a rigorous
enclosure. Exhausting the panel budget produces `converged = false`; equal
endpoints return a zero integral without calling the host, and reversed bounds
reverse the integral's sign. Callback errors propagate with their original
error code and message. The evaluator stops invoking the host after the first
callback failure.

The callback bank snapshots and validates each signature once at registration.
Signatures must be concrete, fixed-arity and unconstrained; structured tensor,
record and variant values are recursively checked at invocation. A callback
handle's type comes from this bank rather than a fabricated function type in
`Value::type_expr()`. Inputs are checked before host code runs, and outputs are
checked before returning them to the numerical kernel. The function type check
uses argument contravariance, return covariance and exact argument count.

A callback may make a nested call into the same `Session`: RSSN releases the
session graph borrow before dispatching a dynamic evaluator. No bank lock or
graph borrow is held while host signature or invocation code runs. A context
must retain the registered callback object for as long as its handle is used.
`NativeCallback` requires `Send + Sync`; RSSN does not make callbacks parallel
implicitly, and operators needing callbacks advertise session-bound execution
and host effects in the manifest. They have no symbolic rewrite or JIT/FFI/GPU
implementation unless the descriptor explicitly provides one.

## Physical dimensions

`num_integrate` declares the polymorphic contract

```
(a: Float[D], b: Float[D], abs_tol: Float[D*E], rel_tol: Float,
 max_panels: Int, f: (Float[D]) -> Float[E])
  -> {value: Float[D*E], error: Float[D*E], evaluations: Int, converged: Bool}
```

Supply concrete callback dimensions. The adapter converts quantities to
coherent magnitudes before invoking the host and validates the callback's
returned dimension. Bounds must have the same dimension, relative tolerance is
dimensionless, and a nonzero absolute tolerance has the integral's dimension.
An exactly zero absolute tolerance is dimension-polymorphic. Tolerances are
nonnegative and cannot both be zero. Unknown dimensions never establish a
concrete callback contract.

This boundary does not turn floating-point magnitudes into exact real values.
A language bridge must preserve the floating representation, closure lifetime,
thread restrictions and concrete type/dimension evidence of its own values.
A language bridge can create a fresh scoped callback handle for each kernel
call; a persistent source-language callback object must remain distinct from
RSSN's context-local raw handle.

## From C, and from the language bindings

A C host has no closure to hand over, so the C ABI registers a plain function
pointer instead (`src/ffi/callback.rs`):

```c
typedef RssnStatus (*RssnCallbackFn)(void *user,
                                     const struct RssnSession *session,
                                     const uint32_t *args, size_t count,
                                     uint32_t *out);
```

`args` are `count` live cell ids holding the arguments; `out` must be set to a
cell holding the result. Both are read and written with the ordinary
`rssn_value_*` readers and constructors, so a callback needs no second API.

| Entry point | What it does |
|---|---|
| `rssn_callback_new(session, signature_json, fn, user, out)` | registers `fn` with an explicit `{"params": [...], "ret": ...}` signature |
| `rssn_callback_new_for(session, op, index, fn, user, out)` | copies the signature `op` declares at argument `index` |
| `rssn_callback_signatures(session, op, index, out_json)` | lists what `op` accepts there, for the few that declare more than one function type |
| `rssn_callback_free(session, cell)` | retires a registration into a reusable slot |

The registration lives in the session's callback bank, and `rssn_value_call`
runs with a clone of that bank, so the cell holding the `Value::Function` is an
ordinary argument. A declared signature that mentions a dimension variable is
instantiated as dimensionless: a C callback sees plain numbers, which is what
the engine passes it anyway. The three bindings build on exactly this surface;
`bindings/examples/` has one program per language.

## Manifest compatibility

The current manifest API version is 2. Operator metadata includes optional
`rust_path` provenance, and callbacks appear as `HostCallback` evaluators.
Binary consumers must negotiate the API version before decoding: adding a
field changes the binary struct layout. JSON readers can default the missing
provenance field when reading older documents. Source provenance describes the
adapter's origin; it does not promise that every backend or Rust signature is
executable through the value ABI.
