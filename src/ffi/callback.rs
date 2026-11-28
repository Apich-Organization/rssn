//! Host callbacks: the C form of a [`HostCallback`].
//!
//! The value ABI in [`super::value`] can carry every [`Value`] the kernel has
//! except one: a [`Value::Function`] names a Rust closure through a
//! [`CallbackRef`], and a C host has no closure to hand over. That single gap
//! is what kept the 189 higher-order operators (`num_bfgs`,
//! `num_analysis_integrate_adaptive_scaled`, `num_bootstrap`, ...) out of
//! reach of every binding, and it is what this module closes.
//!
//! # The protocol
//!
//! A host callback is a C function with this shape:
//!
//! ```c
//! RssnStatus f(void *user, const RssnSession *session,
//!             const uint32_t *args, size_t count, uint32_t *out);
//! ```
//!
//! * `user` is whatever pointer the host registered; rssn never reads it.
//! * `session` is the session the callback was registered on, so the callback
//!   can read its arguments and build its result with the ordinary
//!   `rssn_value_*` readers and constructors.
//! * `args` are `count` live cell ids holding the arguments the operator
//!   passed. They are owned by rssn and are released when the callback
//!   returns; a host that wants to keep one must copy it with
//!   `rssn_value_new_*` or `rssn_value_to_*`.
//! * `out` must be set to a cell id holding the result. The cell may be one
//!   the host created (it stays the host's, and outlives the call) or one of
//!   the argument cells (rssn copies the value out before releasing them).
//! * the return value is an [`RssnStatus`]. Anything but `Ok` aborts the
//!   enclosing operator with that status; call [`rssn_callback_fail`] first to
//!   attach a message, which becomes the [`OpError::Host`] text.
//!
//! Registration stores the callback as a `Value::Function` cell, exactly the
//! value the operator's signature expects, so a host passes that cell as an
//! ordinary argument to [`rssn_value_call`](super::rssn_value_call).
//!
//! ```c
//! uint32_t f;
//! rssn_callback_new_for(session, "num_analysis_integrate_adaptive_scaled", 5,
//!                       square, NULL, &f);
//! uint32_t args[6] = { xs, rules, lo, hi, depth, f };
//! uint32_t out;
//! rssn_value_call(session, "num_analysis_integrate_adaptive_scaled",
//!                 args, 6, &out);
//! ```
//!
//! # Signatures
//!
//! A callback is type-checked against the signature it was registered with,
//! and that signature also participates in overload resolution for the
//! operator that receives it. Two ways to say what it is:
//!
//! * [`rssn_callback_new_for`] copies the declared type of the operator's
//!   function parameter at a given argument position. This is what a binding
//!   uses, because it cannot guess the host's types and the engine already
//!   states them.
//! * [`rssn_callback_new`] takes the signature explicitly, as JSON:
//!   `{"params": ["real"], "ret": "real"}`. The `TypeExpr` spelling is the
//!   manifest's own (`"real"`, `"float"`, `"int"`, `"complex"`, `"bool"`,
//!   `"string"`, `"symbol"`, `{"tensor": {"rank": {"lit": 1}, "shape": null,
//!   "elem": "float"}}`, `{"record": [["value", "float"]]}`, ...).
//!
//! A declared signature that mentions a dimension *variable* (`Real[D]`) is
//! instantiated with the dimensionless dimension before it is handed to the
//! host: a C callback sees plain numbers and tensors of them, which is what
//! the engine passes it anyway (the `_scaled` family normalises first).
//!
//! # Lifetime and threads
//!
//! A callback lives in the session's bank until [`rssn_callback_free`] retires
//! it or the session is freed, whichever comes first. Freeing the cell alone
//! keeps the registration alive; retire the callback when a host builds a
//! fresh closure per iteration. Retired slots are reused, so a long-lived
//! process does not grow the bank.
//!
//! Callbacks run on the thread that entered the computation, and the engine
//! does not call them in parallel (see
//! [`CallCtx::invoke_callback`](crate::kernel::call::CallCtx::invoke_callback)).
//! They may re-enter rssn on the same session — a callback that calls
//! `rssn_value_call` is legal — but they must not outlive the session.

use std::ffi::c_void;
use std::os::raw::c_char;
use std::sync::Arc;

use super::RssnSession;
use super::RssnStatus;
use super::fail;
use super::guard;
use super::last_error_text;
use super::set_error;
use super::to_c;
use crate::kernel::call::OpError;
use crate::kernel::plugin::HostCallback;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::NatExpr;
use crate::kernel::signature::Signature;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::CallbackRef;
use crate::kernel::value::Value;

/// A host function rssn may call.
///
/// `user` is the registered context pointer. `session` owns the `count`
/// argument cells in `args`, which remain valid until the callback returns.
/// Write the result cell ID to `out` and return [`RssnStatus::Ok`] on success.
/// The result can be one of the argument cells or a cell the host created.
pub type RssnCallbackFn = extern "C" fn(
    user: *mut c_void,
    session: *const RssnSession,
    args: *const u32,
    count: usize,
    out: *mut u32,
) -> RssnStatus;

/// One registered host function.
struct HostFn {
    function: RssnCallbackFn,
    user: *mut c_void,
    session: *const RssnSession,
    signature: Signature,
}

// SAFETY: a `HostFn` is reachable from `CallbackBank`, which is `Clone` and
// stored in a `CallCtx`, so the trait requires `Send + Sync`. The raw session
// pointer is only dereferenced while the session is alive and only from the
// thread that entered the computation (the engine never invokes a callback in
// parallel), and the host's `user` pointer is opaque to rssn. The host
// therefore owns the obligation that its callback is safe to call on the
// computing thread and that it does not outlive the session.
unsafe impl Send for HostFn {}
unsafe impl Sync for HostFn {}

impl HostCallback for HostFn {
    fn call(
        &self,
        args: &[Value],
    ) -> Result<Value, OpError> {
        // SAFETY: the session outlives every callback registered on it, and
        // the bank is owned by that session.
        let session = unsafe { self.session.as_ref() }.ok_or_else(|| {
            OpError::Host("the session of this callback has been freed".to_owned())
        })?;
        // The arguments become ordinary cells so the host can read them with
        // the `rssn_value_*` readers it already knows. No borrow is held while
        // host code runs, so a callback may call back into rssn.
        let cells: Vec<u32> = {
            let mut bank = session.values.borrow_mut();
            args.iter()
                .cloned()
                .map(|value| bank.insert(value))
                .collect()
        };
        let mut out = 0u32;
        let status = (self.function)(
            self.user,
            self.session,
            cells.as_ptr(),
            cells.len(),
            &raw mut out,
        );
        // Read the result before releasing the argument cells: a host is
        // allowed to answer with one of its own arguments.
        let result = if status == RssnStatus::Ok {
            let bank = session.values.borrow();
            bank.get(out).cloned()
        } else {
            None
        };
        {
            let mut bank = session.values.borrow_mut();
            for cell in &cells {
                bank.remove(*cell);
            }
        }
        if status != RssnStatus::Ok {
            let message = last_error_text();
            return Err(OpError::Host(if message.is_empty() {
                format!("host callback returned status {}", status as i32)
            } else {
                message
            }));
        }
        result.ok_or_else(|| {
            OpError::Host(if out == 0 {
                "host callback returned no result cell".to_owned()
            } else {
                format!("host callback returned cell {out}, which is not live")
            })
        })
    }

    fn signature(&self) -> Signature {
        self.signature.clone()
    }
}

/// The JSON form [`rssn_callback_new`] reads.
#[derive(serde::Deserialize)]
struct CallbackSignature {
    params: Vec<TypeExpr>,
    ret: TypeExpr,
}

/// Records the message of a failed host callback.
///
/// Call this from inside a callback before returning a non-`Ok` status: the
/// text becomes the `OpError::Host` message, which the enclosing call reports
/// through `rssn_last_error`. The string is copied; free nothing.
///
/// # Safety
/// `message` must be null or a NUL-terminated UTF-8 string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_callback_fail(message: *const c_char) {
    if message.is_null() {
        return;
    }
    // SAFETY: caller contract.
    if let Ok(text) = unsafe { super::text(message) } {
        set_error(text);
    }
}

/// Registers `function` with an explicit signature and writes a cell holding
/// the callback to `*out`.
///
/// `signature` is JSON: `{"params": ["real"], "ret": "real"}`. The
/// [`TypeExpr`](crate::kernel::signature::TypeExpr) spelling follows the manifest.
/// The callback must be concrete (no type parameters, no variadic tail, no
/// dimension variables once instantiated) or registration fails with
/// [`RssnStatus::InvalidArgument`].
///
/// # Safety
/// `session` must be a live session, `signature` a NUL-terminated string and
/// `out` a valid pointer. `function` must be a function pointer that stays
/// valid for as long as the session, and `user` must be safe to pass to it on
/// the computing thread.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_callback_new(
    session: *const RssnSession,
    signature: *const c_char,
    function: RssnCallbackFn,
    user: *mut c_void,
    out: *mut u32,
) -> RssnStatus {
    guard(|| {
        // SAFETY: caller contract.
        let Some(session) = (unsafe { session.as_ref() }) else {
            return fail(RssnStatus::NullArgument, "null session");
        };
        if out.is_null() {
            return fail(RssnStatus::NullArgument, "null output pointer");
        }
        // SAFETY: caller contract.
        let text = match unsafe { super::text(signature) } {
            | Ok(text) => text,
            | Err(status) => return status,
        };
        let parsed: CallbackSignature = match serde_json::from_str(text) {
            | Ok(parsed) => parsed,
            | Err(error) => {
                return fail(
                    RssnStatus::InvalidArgument,
                    &format!("callback signature is not JSON of the form {{params, ret}}: {error}"),
                );
            },
        };
        register(
            session,
            &Signature::new(parsed.params, parsed.ret),
            function,
            user,
            out,
        )
    })
}

/// Registers `function` as the callback for argument `index` of operator `op`,
/// taking the signature from the operator's own declaration.
///
/// This is the call a binding should make: the engine already states that the
/// argument is a `([real]) -> real` (or a tensor-to-tensor, record-to-bool,
/// ...) function, so the host does not have to repeat it. Dimension variables
/// in the declared type are instantiated as dimensionless.
///
/// When the operator has several overloads with different function types at
/// `index` (an objective and a gradient, say), the first one wins;
/// [`rssn_callback_signatures`] lists them all so a host can register the
/// others explicitly with [`rssn_callback_new`].
///
/// # Safety
/// As [`rssn_callback_new`], plus `op` must be a NUL-terminated string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_callback_new_for(
    session: *const RssnSession,
    op: *const c_char,
    index: usize,
    function: RssnCallbackFn,
    user: *mut c_void,
    out: *mut u32,
) -> RssnStatus {
    guard(|| {
        // SAFETY: caller contract.
        let Some(session) = (unsafe { session.as_ref() }) else {
            return fail(RssnStatus::NullArgument, "null session");
        };
        if out.is_null() {
            return fail(RssnStatus::NullArgument, "null output pointer");
        }
        // SAFETY: caller contract.
        let name = match unsafe { super::text(op) } {
            | Ok(name) => name,
            | Err(status) => return status,
        };
        let Some(signature) = candidate_signatures(&session.session, name, index)
            .into_iter()
            .next()
        else {
            return fail(
                RssnStatus::InvalidArgument,
                &format!("operator `{name}` declares no function parameter at index {index}"),
            );
        };
        register(session, &signature, function, user, out)
    })
}

/// Writes the JSON array of signatures a callback for argument `index` of `op`
/// may be registered with; free it with `rssn_string_free`.
///
/// Each element has the shape [`rssn_callback_new`] accepts, so a binding can
/// iterate: register the first, retry with the next when the call reports a
/// callback type mismatch. Dimension variables are already instantiated.
///
/// # Safety
/// `session` must be a live session, `op` a NUL-terminated string and `out` a
/// valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_callback_signatures(
    session: *const RssnSession,
    op: *const c_char,
    index: usize,
    out: *mut *mut c_char,
) -> RssnStatus {
    guard(|| {
        // SAFETY: caller contract.
        let Some(session) = (unsafe { session.as_ref() }) else {
            return fail(RssnStatus::NullArgument, "null session");
        };
        if out.is_null() {
            return fail(RssnStatus::NullArgument, "null output pointer");
        }
        // SAFETY: caller contract.
        let name = match unsafe { super::text(op) } {
            | Ok(name) => name,
            | Err(status) => return status,
        };
        let signatures: Vec<serde_json::Value> = candidate_signatures(
            &session.session,
            name,
            index,
        )
        .into_iter()
        .map(|signature| serde_json::json!({ "params": signature.params, "ret": signature.ret }))
        .collect();
        let json = serde_json::Value::Array(signatures).to_string();
        let text = to_c(json);
        if text.is_null() {
            return fail(
                RssnStatus::Panic,
                "the signature list is not representable as C text",
            );
        }
        // SAFETY: `out` is non-null by the check above.
        unsafe { *out = text };
        RssnStatus::Ok
    })
}

/// Releases the callback named by the cell `cell_id` and the cell itself.
///
/// The registration and the cell are separate: freeing the cell keeps the
/// callback alive (an operator may still hold it in a term it kept), so this
/// call is the one to use when a host builds a fresh closure per iteration.
/// Retiring an already-retired callback is not an error.
///
/// # Safety
/// `session` must be a live session.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_callback_free(
    session: *const RssnSession,
    cell_id: u32,
) {
    let _ = guard(|| {
        // SAFETY: caller contract.
        let Some(session) = (unsafe { session.as_ref() }) else {
            return fail(RssnStatus::NullArgument, "null session");
        };
        let retired = session.values.borrow_mut().take(cell_id).and_then(|value| {
            match value {
                | Value::Function(CallbackRef(index)) => Some(index),
                | _ => None,
            }
        });
        if let Some(index) = retired {
            session.callbacks.borrow_mut().retire(index as usize);
        }
        RssnStatus::Ok
    });
}

/// The number of live callbacks in `session`; a diagnostic for leak hunting.
///
/// # Safety
/// `session` must be a live session.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_callback_count(session: *const RssnSession) -> usize {
    // SAFETY: caller contract.
    unsafe { session.as_ref() }.map_or(0, |session| session.callbacks.borrow().len())
}

/// Stores `Value::Function(handle)` in a new cell.
fn store(
    session: &RssnSession,
    handle: CallbackRef,
    out: *mut u32,
) -> RssnStatus {
    let cell = session.values.borrow_mut().insert(Value::Function(handle));
    // SAFETY: `out` was checked for null by the caller.
    unsafe { *out = cell };
    RssnStatus::Ok
}

/// Adds `signature` to the session's bank and stores the resulting value.
fn register(
    session: &RssnSession,
    signature: &Signature,
    function: RssnCallbackFn,
    user: *mut c_void,
    out: *mut u32,
) -> RssnStatus {
    let signature = Signature::new(
        signature.params.iter().map(erase_type).collect(),
        erase_type(&signature.ret),
    );
    let callback = Arc::new(HostFn {
        function,
        user,
        session: std::ptr::from_ref::<RssnSession>(session),
        signature,
    });
    match session.callbacks.borrow_mut().add(callback) {
        | Ok(handle) => store(session, handle, out),
        | Err(error) => {
            fail(
                RssnStatus::InvalidArgument,
                &format!("callback signature is not registrable: {error}"),
            )
        },
    }
}

/// Every distinct function signature `op` declares at argument `index`, in
/// overload order, with dimension variables instantiated.
fn candidate_signatures(
    session: &crate::api::Session,
    op: &str,
    index: usize,
) -> Vec<Signature> {
    let mut found: Vec<Signature> = Vec::new();
    for (_, signature) in session.overloads(op) {
        let Some(TypeExpr::Function(params, ret)) = signature.params.get(index) else {
            continue;
        };
        let candidate = Signature::new(params.iter().map(erase_type).collect(), erase_type(ret));
        if !found.contains(&candidate) {
            found.push(candidate);
        }
    }
    found
}

/// Replaces every dimension variable with the dimensionless dimension.
///
/// A host callback is a concrete function: it receives plain numbers, not a
/// type parameter. The engine passes it values whose dimensions it has already
/// normalised (the `_scaled` family) or whose dimension is a bound variable, so
/// erasing the variable to `Dimless` is both the most permissive and the
/// truthful description of what the host sees.
fn erase_type(ty: &TypeExpr) -> TypeExpr {
    match ty {
        | TypeExpr::Dim(inner, dim) => {
            match erase_dim(dim) {
                // `Real[dimensionless]` is `Real` to the matcher, so say so:
                // a host that reads the listed signature should see the plain
                // scalar it is actually handed.
                | DimExpr::Dimless => erase_type(inner),
                | dim => TypeExpr::Dim(Box::new(erase_type(inner)), dim),
            }
        },
        | TypeExpr::Unit(dim) => TypeExpr::Unit(erase_dim(dim)),
        | TypeExpr::Point(dim) => TypeExpr::Point(erase_dim(dim)),
        | TypeExpr::Tensor { rank, shape, elem } => {
            TypeExpr::Tensor {
                rank: erase_nat(rank),
                shape: shape
                    .as_ref()
                    .map(|shape| shape.iter().map(erase_nat).collect()),
                elem: Box::new(erase_type(elem)),
            }
        },
        | TypeExpr::Sparse { format, rank, elem } => {
            TypeExpr::Sparse {
                format: format.clone(),
                rank: erase_nat(rank),
                elem: Box::new(erase_type(elem)),
            }
        },
        | TypeExpr::Record(fields) => {
            TypeExpr::Record(
                fields
                    .iter()
                    .map(|(name, field)| (name.clone(), erase_type(field)))
                    .collect(),
            )
        },
        | TypeExpr::Variant(variants) => {
            TypeExpr::Variant(
                variants
                    .iter()
                    .map(|(tag, fields)| (tag.clone(), fields.iter().map(erase_type).collect()))
                    .collect(),
            )
        },
        | TypeExpr::Function(params, ret) => {
            TypeExpr::Function(
                params.iter().map(erase_type).collect(),
                Box::new(erase_type(ret)),
            )
        },
        | other => other.clone(),
    }
}

/// A dimension with its variables replaced by `Dimless`.
fn erase_dim(dim: &DimExpr) -> DimExpr {
    match dim {
        | DimExpr::Dimless | DimExpr::Base(_) => dim.clone(),
        | DimExpr::Mul(parts) => DimExpr::Mul(parts.iter().map(erase_dim).collect()),
        | DimExpr::Inverse(inner) => DimExpr::Inverse(Box::new(erase_dim(inner))),
        | DimExpr::Pow(inner, exponent) => {
            match exponent {
                | crate::kernel::signature::Exp::Rat(_) => {
                    DimExpr::Pow(Box::new(erase_dim(inner)), exponent.clone())
                },
                // A symbolic exponent has no concrete value to erase to.
                | crate::kernel::signature::Exp::Symbol(_) => DimExpr::Dimless,
            }
        },
        | DimExpr::Any | DimExpr::NonMetric | DimExpr::Container | DimExpr::Var(_) => {
            DimExpr::Dimless
        },
    }
}

/// A natural-number expression with its variables replaced by `1`.
fn erase_nat(nat: &NatExpr) -> NatExpr {
    match nat {
        | NatExpr::Lit(_) => nat.clone(),
        | NatExpr::Var(_) => NatExpr::Lit(1),
        | NatExpr::Sum(a, b) => NatExpr::Sum(Box::new(erase_nat(a)), Box::new(erase_nat(b))),
        | NatExpr::Product(a, b) => {
            NatExpr::Product(Box::new(erase_nat(a)), Box::new(erase_nat(b)))
        },
    }
}
