# Checked host callbacks

RSSN's dynamic ABI supports serial, synchronous host callbacks. Register a
`HostCallback` in a `CallCtx`, then pass its returned handle to `Session::call_ho`.
The handles are appended after the ordinary arguments in their supplied order;
they are ordinary function parameters in the operator's declared signature.

```rust
use std::sync::Arc;
use rssn::api::Session;
use rssn::kernel::{CallCtx, Value};
use rssn::kernel::plugin::NativeCallback;
use rssn::kernel::signature::{Signature, TypeExpr};

let session = Session::new();
let mut context = CallCtx::new();
let square = context.add_callback(Arc::new(NativeCallback::new(
    |args: &[Value]| {
        let x = args[0].as_f64().unwrap();
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
# Ok::<(), rssn::kernel::OpError>(())
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

## Manifest compatibility

The current manifest API version is 2. Operator metadata includes optional
`rust_path` provenance, and callbacks appear as `HostCallback` evaluators.
Binary consumers must negotiate the API version before decoding: adding a
field changes the binary struct layout. JSON readers can default the missing
provenance field when reading older documents. Source provenance describes the
adapter's origin; it does not promise that every backend or Rust signature is
executable through the value ABI.
