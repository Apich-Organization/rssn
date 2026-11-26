//! Dynamic values: the C form of [`kernel::Value`](crate::kernel::Value).
//!
//! The term ABI in [`super::term`] can only carry values that are nodes of the
//! expression graph. An operator whose result is a native object — an AD tape,
//! a compiled backend, a solver handle — returns a [`Value::Opaque`] instead,
//! and an operator that takes such an object as an argument can only be called
//! with one. Those operators are the ones the manifest marks
//! `backend_support.interpreter` but whose arguments have no term form, and
//! this module is how they become reachable from C. A function argument is the
//! same story: it travels as a `Value::Function`, which [`super::callback`]
//! registers from a C function pointer.
//!
//! # Cells
//!
//! A [`Value`] does not cross the boundary by value: it is stored in a
//! session-owned *cell* and addressed by a `uint32_t` id. `0` is never a valid
//! cell. A cell lives until [`rssn_value_free`] or until the session is freed,
//! and only the session that minted it can read it.
//!
//! The usual shape of a native-object session is
//!
//! ```c
//! uint32_t backend;
//! rssn_value_call(session, "num_backend_current", NULL, 0, &backend);
//! uint32_t name_cell;
//! rssn_value_call(session, "num_backend_name", &backend, 1, &name_cell);
//! char *name;
//! rssn_value_to_text(session, name_cell, &name);      /* rssn_string_free(name) */
//! ```
//!
//! # Inspecting a value
//!
//! [`rssn_value_kind`] gives the variant, [`rssn_value_to_f64`],
//! [`rssn_value_to_text`], [`rssn_value_to_tensor`], [`rssn_value_len`] and
//! [`rssn_value_at`] read it, [`rssn_value_to_term`] turns it into a term when
//! the value has a graph form, and [`rssn_value_json`] serialises the whole
//! value for a host that would rather parse than branch.
//!
//! # Errors
//!
//! Every entry point reports an [`RssnStatus`]; the message is in
//! `rssn_last_error`. A call that fails to reduce reports the code of the
//! [`OpError`] that stopped it: [`RssnStatus::UnknownOperator`] for an unknown
//! name or a wrong argument count, [`RssnStatus::InvalidArgument`] for a type
//! mismatch or an out-of-domain value, [`RssnStatus::Dimension`] for a strict
//! dimensional conflict.

use std::collections::BTreeMap;
use std::os::raw::c_char;
use std::os::raw::c_int;

use super::RssnSession;
use super::RssnStatus;
use super::RssnTensor;
use super::fail;
use super::guard;
use super::to_c;
use crate::kernel::call::CallCtx;
use crate::kernel::call::OpError;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use crate::kernel::value::Value;
use crate::kernel::value::ValueKind;

/// A session-owned table of dynamic values.
#[derive(Debug)]
pub(crate) struct ValueBank {
    next: u32,
    cells: BTreeMap<u32, Value>,
}

impl ValueBank {
    pub(crate) fn new() -> Self {
        Self {
            next: 1,
            cells: BTreeMap::new(),
        }
    }

    /// Stores `value` and returns its cell id. Ids are never reused, so a
    /// stale id can never name a different value.
    pub(crate) fn insert(
        &mut self,
        value: Value,
    ) -> u32 {
        let id = self.next;
        self.next = self.next.saturating_add(1);
        self.cells.insert(id, value);
        id
    }

    pub(crate) fn get(
        &self,
        id: u32,
    ) -> Option<&Value> {
        self.cells.get(&id)
    }

    /// Removes and returns the value in cell `id`.
    pub(crate) fn take(
        &mut self,
        id: u32,
    ) -> Option<Value> {
        self.cells.remove(&id)
    }

    pub(crate) fn remove(
        &mut self,
        id: u32,
    ) -> bool {
        self.cells.remove(&id).is_some()
    }

    pub(crate) fn len(&self) -> usize {
        self.cells.len()
    }
}

/// The stable code of a [`ValueKind`], as [`rssn_value_kind`] reports it.
const fn kind_code(kind: ValueKind) -> u8 {
    match kind {
        | ValueKind::Int => 0,
        | ValueKind::Real => 1,
        | ValueKind::Float => 2,
        | ValueKind::Complex => 3,
        | ValueKind::Bool => 4,
        | ValueKind::String => 5,
        | ValueKind::Symbol => 6,
        | ValueKind::Unit => 7,
        | ValueKind::Quantity => 8,
        | ValueKind::Tensor => 9,
        | ValueKind::Sparse => 10,
        | ValueKind::Record => 11,
        | ValueKind::Variant => 12,
        | ValueKind::Term => 13,
        | ValueKind::Function => 14,
        | ValueKind::Dyn => 15,
        | ValueKind::Opaque => 16,
    }
}

/// The status of a failed dynamic call.
const fn status_of(error: &OpError) -> RssnStatus {
    match error {
        | OpError::UnknownOp(_) | OpError::Arity { .. } | OpError::Ambiguous { .. } => {
            RssnStatus::UnknownOperator
        },
        | OpError::Dimension(_) => RssnStatus::Dimension,
        | OpError::NotNumeric(_) => RssnStatus::NotNumeric,
        | OpError::ForeignTerm { .. } => RssnStatus::NoSuchTerm,
        | OpError::BudgetExceeded => RssnStatus::Cancelled,
        | OpError::Internal(_) => RssnStatus::Panic,
        | OpError::NotReduced(_)
        | OpError::Domain(_)
        | OpError::Type { .. }
        | OpError::Unsupported { .. }
        | OpError::Plugin(_)
        | OpError::Codec(_)
        | OpError::CallbackType { .. }
        | OpError::HostDiagnostic { .. }
        | OpError::Host(_) => RssnStatus::InvalidArgument,
    }
}

/// Records `error` and returns its status.
fn op_error(error: &OpError) -> RssnStatus {
    fail(status_of(error), &error.to_string())
}

/// The context a value call runs with: the symbolic target, the default
/// budget, advisory dimensions, and the host callbacks registered on the
/// session.
///
/// The bank is cloned (a vector of `Arc`s) rather than borrowed, so no borrow
/// is held while the operator runs: a callback that calls back into rssn on
/// the same session is legal.
fn context(session: &RssnSession) -> CallCtx<'static> {
    let mut ctx = CallCtx::new();
    ctx.callbacks = session.callbacks.borrow().clone();
    ctx
}

/// Borrows a session together with its value bank.
///
/// # Safety
/// `session` must be null or a live session pointer.
unsafe fn borrow<'a>(session: *const RssnSession) -> Result<&'a RssnSession, RssnStatus> {
    // SAFETY: caller contract.
    unsafe { session.as_ref() }.ok_or_else(|| fail(RssnStatus::NullArgument, "null session"))
}

/// Runs `body` with the session behind `raw`, mapping its result to a status.
fn with_session(
    raw: *const RssnSession,
    body: impl FnOnce(&RssnSession) -> Result<(), RssnStatus>,
) -> RssnStatus {
    guard(|| {
        // SAFETY: caller contract.
        match unsafe { borrow(raw) }.and_then(|session| body(session)) {
            | Ok(()) => RssnStatus::Ok,
            | Err(status) => status,
        }
    })
}

/// Runs `body` on the value in cell `id`. The borrow ends with `body`.
fn with_cell<T>(
    session: &RssnSession,
    id: u32,
    body: impl FnOnce(&Value) -> Result<T, RssnStatus>,
) -> Result<T, RssnStatus> {
    let bank = session.values.borrow();
    let value = bank.get(id).ok_or_else(|| {
        fail(
            RssnStatus::NoValue,
            "value cell does not belong to this session",
        )
    })?;
    body(value)
}

/// Copies the values in the `count` cells `ids`.
fn collect_values(
    session: &RssnSession,
    ids: &[u32],
) -> Result<Vec<Value>, RssnStatus> {
    let bank = session.values.borrow();
    let mut values = Vec::with_capacity(ids.len());
    for &id in ids {
        let value = bank
            .get(id)
            .ok_or_else(|| fail(RssnStatus::NoValue, "argument is not a live cell"))?;
        values.push(value.clone());
    }
    Ok(values)
}

/// Stores `value` in the session's bank and writes its cell id to `out`.
fn store(
    session: &RssnSession,
    value: Value,
    out: *mut u32,
) -> Result<(), RssnStatus> {
    if out.is_null() {
        return Err(fail(RssnStatus::NullArgument, "null output pointer"));
    }
    let id = session.values.borrow_mut().insert(value);
    // SAFETY: `out` is non-null by the check above.
    unsafe { *out = id };
    Ok(())
}

/// A rank-one tensor of `values`, typed with their common type when they have
/// one (the term type otherwise, which every operator accepts).
fn list_value(values: Vec<Value>) -> Value {
    let elem = values
        .first()
        .map(Value::type_expr)
        .filter(|first| values.iter().all(|value| value.type_expr() == *first))
        .unwrap_or_else(crate::kernel::infer::term);
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values,
        elem,
    })
}

/// Fills a `RssnTensor` with owned buffers, as `rssn_tensor_free` expects.
///
/// # Safety
/// `out` must be a valid pointer.
unsafe fn tensor_into(
    out: *mut RssnTensor,
    shape: &[usize],
    data: &[f64],
) {
    let rank = shape.len();
    let len = data.len();
    let shape = Box::into_raw(shape.to_vec().into_boxed_slice()).cast::<usize>();
    let data = Box::into_raw(data.to_vec().into_boxed_slice()).cast::<f64>();
    // SAFETY: non-null and valid by contract.
    unsafe {
        *out = RssnTensor {
            rank,
            shape,
            len,
            data,
        }
    };
}

/// The text of a string or symbol value; any other value's display form.
fn text_of(value: &Value) -> String {
    match value {
        | Value::String(text) | Value::Symbol(text) => text.clone(),
        | other => other.to_string(),
    }
}

/// The number a value denotes, looking through the wrappers that carry one.
///
/// A `Dyn` is a number plus a dimension and a `Quantity` is a number plus a
/// unit; both are what an operator hands a host callback for a dimensioned
/// parameter, and both have a plain numeric part. Reporting that part is what
/// makes a callback able to read its own argument without a dimension reader
/// the ABI does not have.
fn number_of(value: &Value) -> Option<f64> {
    match value {
        | Value::Dyn(inner) => number_of(&inner.value),
        | Value::Quantity(quantity) => number_of(&quantity.value),
        | other => other.as_f64(),
    }
}

/// Looks through a `Dyn` or `Quantity` wrapper around a container.
fn unwrap_dimensioned(value: &Value) -> &Value {
    match value {
        | Value::Dyn(inner) => unwrap_dimensioned(&inner.value),
        | Value::Quantity(quantity) => unwrap_dimensioned(&quantity.value),
        | other => other,
    }
}

// ---------------------------------------------------------------------------
// construction
// ---------------------------------------------------------------------------

/// Stores the exact integer `value` and writes its cell id to `*out`.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_new_int(
    session: *const RssnSession,
    value: i64,
    out: *mut u32,
) -> RssnStatus {
    with_session(session, |session| {
        store(session, Value::Int(value.into()), out)
    })
}

/// Stores the float `value` and writes its cell id to `*out`.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_new_float(
    session: *const RssnSession,
    value: f64,
    out: *mut u32,
) -> RssnStatus {
    with_session(session, |session| store(session, Value::Float(value), out))
}

/// Stores the text `value` and writes its cell id to `*out`.
///
/// # Safety
/// `session` must be a live session, `value` a NUL-terminated string, and
/// `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_new_text(
    session: *const RssnSession,
    value: *const c_char,
    out: *mut u32,
) -> RssnStatus {
    with_session(session, |session| {
        // SAFETY: caller contract.
        let text = unsafe { super::text(value) }?.to_owned();
        store(session, Value::String(text), out)
    })
}

/// Stores the free symbol `name` and writes its cell id to `*out`.
///
/// # Safety
/// `session` must be a live session, `name` a NUL-terminated string, and
/// `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_new_symbol(
    session: *const RssnSession,
    name: *const c_char,
    out: *mut u32,
) -> RssnStatus {
    with_session(session, |session| {
        // SAFETY: caller contract.
        let text = unsafe { super::text(name) }?.to_owned();
        store(session, Value::Symbol(text), out)
    })
}

/// Stores the term `handle` as a value and writes its cell id to `*out`.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_new_term(
    session: *const RssnSession,
    handle: u32,
    out: *mut u32,
) -> RssnStatus {
    with_session(session, |session| {
        let term = super::resolve(&session.session, handle)?;
        store(session, Value::Term(term.handle()), out)
    })
}

/// Stores a list of the `count` cells `cells` and writes its cell id to `*out`.
///
/// The elements keep their individual kinds; a list of a float and a symbol is
/// a rank-one tensor with a dynamic element type, which is what an operator
/// with a `list(...)` parameter reads.
///
/// # Safety
/// `cells` must hold `count` cell ids (null is allowed only when `count` is 0),
/// `session` must be a live session, and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_new_list(
    session: *const RssnSession,
    cells: *const u32,
    count: usize,
    out: *mut u32,
) -> RssnStatus {
    with_session(session, |session| {
        if count > 0 && cells.is_null() {
            return Err(fail(RssnStatus::NullArgument, "null element array"));
        }
        // SAFETY: `cells` holds `count` ids by the caller's contract; a null
        // pointer is allowed only when `count` is 0, and `from_raw_parts`
        // requires a non-null pointer even then.
        let ids: &[u32] = if count == 0 {
            &[]
        } else {
            unsafe { std::slice::from_raw_parts(cells, count) }
        };
        let values = collect_values(session, ids)?;
        store(session, list_value(values), out)
    })
}

/// Stores a dense `f64` tensor and writes its cell id to `*out`.
///
/// # Safety
/// `data` must hold `len` values, `shape` must hold `rank` extents whose
/// product is `len`, `session` must be a live session, and `out` a valid
/// pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_new_tensor(
    session: *const RssnSession,
    data: *const f64,
    len: usize,
    shape: *const usize,
    rank: usize,
    out: *mut u32,
) -> RssnStatus {
    with_session(session, |session| {
        if len > 0 && data.is_null() {
            return Err(fail(RssnStatus::NullArgument, "null tensor data"));
        }
        if rank > 0 && shape.is_null() {
            return Err(fail(RssnStatus::NullArgument, "null tensor shape"));
        }
        // SAFETY: caller contract; a null pointer is allowed only for an empty
        // extent list, and `from_raw_parts` requires a non-null pointer.
        let values: &[f64] = if len == 0 {
            &[]
        } else {
            unsafe { std::slice::from_raw_parts(data, len) }
        };
        let extents: &[usize] = if rank == 0 {
            &[]
        } else {
            unsafe { std::slice::from_raw_parts(shape, rank) }
        };
        if rank > 0 && extents.iter().product::<usize>() != len {
            return Err(fail(
                RssnStatus::InvalidArgument,
                "tensor shape does not match the data length",
            ));
        }
        let tensor = TensorValue {
            shape: extents.to_vec(),
            data: values.iter().copied().map(Value::Float).collect(),
            // The data is `f64`, so the element type is `Float`, not the exact
            // `Real`: an operator that declares `Tensor<Float>` (every numeric
            // kernel does) accepts it, and one that declares `Tensor<Real>`
            // still does through the widening arm of `matches_type`.
            elem: TypeExpr::Float,
        };
        store(session, Value::Tensor(tensor), out)
    })
}

/// Stores a truth value (`0` false, anything else true) and writes its cell id
/// to `*out`.
///
/// A `bool` has no numeric form an operator would accept in its place, so a
/// host callback that answers a `Bool` parameter needs this constructor.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_new_bool(
    session: *const RssnSession,
    value: c_int,
    out: *mut u32,
) -> RssnStatus {
    with_session(session, |session| {
        store(session, Value::Bool(value != 0), out)
    })
}

/// Stores the complex number `re + im*i` and writes its cell id to `*out`.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_new_complex(
    session: *const RssnSession,
    re: f64,
    im: f64,
    out: *mut u32,
) -> RssnStatus {
    with_session(session, |session| {
        store(
            session,
            Value::Complex(crate::kernel::value::ComplexValue::Float { re, im }),
            out,
        )
    })
}

/// Stores a record of the `count` fields `names[i]` / `cells[i]` and writes its
/// cell id to `*out`.
///
/// Field order is preserved; the kernel matches records by field name, so an
/// operator's `record{energy: real, ...}` parameter is satisfied by the same
/// names in any order. A name may not repeat.
///
/// # Safety
/// `names` must hold `count` NUL-terminated strings, `cells` must hold `count`
/// live cell ids, `session` must be a live session, and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_new_record(
    session: *const RssnSession,
    names: *const *const c_char,
    cells: *const u32,
    count: usize,
    out: *mut u32,
) -> RssnStatus {
    with_session(session, |session| {
        if count > 0 && (names.is_null() || cells.is_null()) {
            return Err(fail(RssnStatus::NullArgument, "null record field array"));
        }
        // SAFETY: the caller promises `count` names and `count` cells.
        let name_slice: &[*const c_char] = if count == 0 {
            &[]
        } else {
            unsafe { std::slice::from_raw_parts(names, count) }
        };
        let ids: &[u32] = if count == 0 {
            &[]
        } else {
            unsafe { std::slice::from_raw_parts(cells, count) }
        };
        let values = collect_values(session, ids)?;
        let mut fields = Vec::with_capacity(count);
        for (name, value) in name_slice.iter().zip(values) {
            // SAFETY: caller contract for each name.
            let name = unsafe { super::text(*name) }?.to_owned();
            if fields.iter().any(|(known, _)| known == &name) {
                return Err(fail(
                    RssnStatus::InvalidArgument,
                    &format!("duplicate record field `{name}`"),
                ));
            }
            fields.push((name, value));
        }
        store(
            session,
            Value::Record(crate::kernel::value::RecordValue { fields }),
            out,
        )
    })
}

// ---------------------------------------------------------------------------
// calling
// ---------------------------------------------------------------------------
/// Calls the operator `op` with the `count` argument cells and stores the
/// result as a new cell.
///
/// This is [`Session::call_value`](crate::api::Session::call_value): the
/// operator is resolved by name and argument count, an operator with a dynamic
/// evaluator is run by it, and every other operator is built into a term and
/// reduced. It is the only way to call an operator whose arguments have no
/// term form — a native object handle, a record, a host callback — because
/// those travel as cells and stay owned by the session.
///
/// The call is symbolic, with the default budget and advisory dimensions.
/// Convert the result with [`rssn_value_to_term`] and compute it with a
/// configuration when a numeric target is wanted.
///
/// # Safety
/// `op` must be a NUL-terminated string, `args` must hold `count` cell ids
/// (null is allowed only when `count` is 0), `session` must be a live session,
/// and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_call(
    session: *const RssnSession,
    op: *const c_char,
    args: *const u32,
    count: usize,
    out: *mut u32,
) -> RssnStatus {
    with_session(session, |session| {
        // SAFETY: caller contract.
        let name = unsafe { super::text(op) }?;
        if count > 0 && args.is_null() {
            return Err(fail(RssnStatus::NullArgument, "null argument array"));
        }
        // SAFETY: `args` holds `count` ids by the caller's contract; a null
        // pointer is allowed only when `count` is 0.
        let ids: &[u32] = if count == 0 {
            &[]
        } else {
            unsafe { std::slice::from_raw_parts(args, count) }
        };
        let values = collect_values(session, ids)?;
        let ctx = context(session);
        let result = session
            .session
            .call_value(name, &values, &ctx)
            .map_err(|error| op_error(&error))?;
        store(session, result, out)
    })
}

/// Calls the operator `op` and writes the result as JSON text.
///
/// The JSON is the `serde` form of a value: a native object appears as
/// `{"opaque": {...}}`, a tensor as `{"tensor": {"shape": [...], "data": [...]}}`.
/// Free the string with `rssn_string_free`.
///
/// # Safety
/// As [`rssn_value_call`], plus `out_json` must be a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_call_json(
    session: *const RssnSession,
    op: *const c_char,
    args: *const u32,
    count: usize,
    out_json: *mut *mut c_char,
) -> RssnStatus {
    with_session(session, |session| {
        if out_json.is_null() {
            return Err(fail(RssnStatus::NullArgument, "null output pointer"));
        }
        // SAFETY: caller contract.
        let name = unsafe { super::text(op) }?;
        if count > 0 && args.is_null() {
            return Err(fail(RssnStatus::NullArgument, "null argument array"));
        }
        // SAFETY: `args` holds `count` ids by the caller's contract; a null
        // pointer is allowed only when `count` is 0.
        let ids: &[u32] = if count == 0 {
            &[]
        } else {
            unsafe { std::slice::from_raw_parts(args, count) }
        };
        let values = collect_values(session, ids)?;
        let ctx = context(session);
        let result = session
            .session
            .call_value(name, &values, &ctx)
            .map_err(|error| op_error(&error))?;
        let json = serde_json::to_string(&result)
            .map_err(|error| fail(RssnStatus::InvalidArgument, &error.to_string()))?;
        let text = to_c(json);
        if text.is_null() {
            return Err(fail(
                RssnStatus::Panic,
                "the JSON is not representable as C text",
            ));
        }
        // SAFETY: `out_json` is non-null by the check above.
        unsafe { *out_json = text };
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// inspection
// ---------------------------------------------------------------------------

/// Writes the kind code of the value in `cell_id` to `*out`.
///
/// The codes are the declaration order of `ValueKind`: `0` int, `1` real,
/// `2` float, `3` complex, `4` bool, `5` string, `6` symbol, `7` unit,
/// `8` quantity, `9` tensor, `10` sparse, `11` record, `12` variant,
/// `13` term, `14` function, `15` dyn, `16` opaque.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_kind(
    session: *const RssnSession,
    cell_id: u32,
    out: *mut u8,
) -> RssnStatus {
    with_session(session, |session| {
        if out.is_null() {
            return Err(fail(RssnStatus::NullArgument, "null output pointer"));
        }
        let code = with_cell(session, cell_id, |value| {
            Ok(kind_code(unwrap_dimensioned(value).kind()))
        })?;
        // SAFETY: `out` is non-null by the check above.
        unsafe { *out = code };
        Ok(())
    })
}

/// Writes the numeric value of `cell_id` to `*out` when it has one.
///
/// A dimensioned scalar (`dyn`) and a quantity are read as the number they
/// carry; `rssn_value_kind` still says which wrapper it was.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_to_f64(
    session: *const RssnSession,
    cell_id: u32,
    out: *mut f64,
) -> RssnStatus {
    with_session(session, |session| {
        if out.is_null() {
            return Err(fail(RssnStatus::NullArgument, "null output pointer"));
        }
        let number = with_cell(session, cell_id, |value| {
            number_of(value)
                .ok_or_else(|| fail(RssnStatus::NotNumeric, "the value has no numeric form"))
        })?;
        // SAFETY: `out` is non-null by the check above.
        unsafe { *out = number };
        Ok(())
    })
}

/// Writes the text of `cell_id` — a string, a symbol, or any value's display
/// form — as a caller-owned string; free it with `rssn_string_free`.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_to_text(
    session: *const RssnSession,
    cell_id: u32,
    out: *mut *mut c_char,
) -> RssnStatus {
    with_session(session, |session| {
        if out.is_null() {
            return Err(fail(RssnStatus::NullArgument, "null output pointer"));
        }
        let text = with_cell(session, cell_id, |value| Ok(text_of(value)))?;
        let text = to_c(text);
        if text.is_null() {
            return Err(fail(
                RssnStatus::Panic,
                "the text is not representable as C text",
            ));
        }
        // SAFETY: `out` is non-null by the check above.
        unsafe { *out = text };
        Ok(())
    })
}

/// Stores the graph form of `cell_id` as a term and writes its handle to `*out`.
///
/// Fails with [`RssnStatus::NoValue`] when the value has no term form — a
/// native object, a record, a host callback. The manifest's `value_kinds` lists
/// which kinds do, under `in_term_graph`.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_to_term(
    session: *const RssnSession,
    cell_id: u32,
    out: *mut u32,
) -> RssnStatus {
    with_session(session, |session| {
        if out.is_null() {
            return Err(fail(RssnStatus::NullArgument, "null output pointer"));
        }
        let value = with_cell(session, cell_id, |value| Ok(value.clone()))?;
        let term = session.session.term_of(&value).map_err(|_| {
            fail(
                RssnStatus::NoValue,
                "the value has no term form; read it with rssn_value_to_text or rssn_value_json",
            )
        })?;
        // SAFETY: `out` is non-null by the check above.
        unsafe { *out = term.id() };
        Ok(())
    })
}

/// Extracts `cell_id` as a dense `f64` array; release it with `rssn_tensor_free`.
///
/// A tensor value is flattened row-major; a scalar value gives rank 0.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_to_tensor(
    session: *const RssnSession,
    cell_id: u32,
    out: *mut RssnTensor,
) -> RssnStatus {
    with_session(session, |session| {
        if out.is_null() {
            return Err(fail(RssnStatus::NullArgument, "null output pointer"));
        }
        let (shape, data) = with_cell(session, cell_id, |value| {
            match unwrap_dimensioned(value) {
                | Value::Tensor(tensor) => {
                    let mut data = Vec::with_capacity(tensor.data.len());
                    for element in &tensor.data {
                        data.push(number_of(element).ok_or_else(|| {
                            fail(
                                RssnStatus::NotNumeric,
                                "the tensor has a non-numeric element",
                            )
                        })?);
                    }
                    Ok((tensor.shape.clone(), data))
                },
                | Value::Term(handle) => {
                    session
                        .session
                        .term_by_id(handle.node)
                        .and_then(|term| term.as_tensor())
                        .ok_or_else(|| {
                            fail(
                                RssnStatus::NoValue,
                                "the term is not a rectangular list of reals",
                            )
                        })
                },
                | other => {
                    let number = number_of(other)
                        .ok_or_else(|| fail(RssnStatus::NoValue, "the value is not a tensor"))?;
                    Ok((Vec::new(), vec![number]))
                },
            }
        })?;
        // SAFETY: `out` is valid by the check above.
        unsafe { tensor_into(out, &shape, &data) };
        Ok(())
    })
}

/// Writes the number of elements of `cell_id` (a tensor's length, a record's
/// field count) to `*out`.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_len(
    session: *const RssnSession,
    cell_id: u32,
    out: *mut usize,
) -> RssnStatus {
    with_session(session, |session| {
        if out.is_null() {
            return Err(fail(RssnStatus::NullArgument, "null output pointer"));
        }
        let length = with_cell(session, cell_id, |value| {
            match unwrap_dimensioned(value) {
                | Value::Tensor(tensor) => Ok(tensor.data.len()),
                | Value::Record(record) => Ok(record.fields.len()),
                | Value::Variant(variant) => Ok(variant.payload.len()),
                | _ => Err(fail(RssnStatus::NoValue, "the value has no elements")),
            }
        })?;
        // SAFETY: `out` is non-null by the check above.
        unsafe { *out = length };
        Ok(())
    })
}

/// Stores element `index` of `cell_id` as a new cell and writes its id to `*out`.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_at(
    session: *const RssnSession,
    cell_id: u32,
    index: usize,
    out: *mut u32,
) -> RssnStatus {
    with_session(session, |session| {
        let element = with_cell(session, cell_id, |value| {
            Ok(match unwrap_dimensioned(value) {
                | Value::Tensor(tensor) => tensor.data.get(index).cloned(),
                | Value::Record(record) => record.fields.get(index).map(|(_, v)| v.clone()),
                | Value::Variant(variant) => variant.payload.get(index).cloned(),
                | _ => None,
            })
        })?;
        let element = element.ok_or_else(|| fail(RssnStatus::NoValue, "index out of range"))?;
        store(session, element, out)
    })
}

/// Writes the value in `cell_id` as JSON text; free it with `rssn_string_free`.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_json(
    session: *const RssnSession,
    cell_id: u32,
    out: *mut *mut c_char,
) -> RssnStatus {
    with_session(session, |session| {
        if out.is_null() {
            return Err(fail(RssnStatus::NullArgument, "null output pointer"));
        }
        let value = with_cell(session, cell_id, |value| Ok(value.clone()))?;
        let json = serde_json::to_string(&value)
            .map_err(|error| fail(RssnStatus::InvalidArgument, &error.to_string()))?;
        let text = to_c(json);
        if text.is_null() {
            return Err(fail(
                RssnStatus::Panic,
                "the JSON is not representable as C text",
            ));
        }
        // SAFETY: `out` is non-null by the check above.
        unsafe { *out = text };
        Ok(())
    })
}

/// Writes the name of field `index` of a record to `*out` as a caller-owned
/// string; free it with `rssn_string_free`.
///
/// Together with [`rssn_value_len`] and [`rssn_value_at`] this reads a record
/// without going through JSON: the length is the field count, the name is the
/// label and `_at` is the value, all in declaration order.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_field_name(
    session: *const RssnSession,
    cell_id: u32,
    index: usize,
    out: *mut *mut c_char,
) -> RssnStatus {
    with_session(session, |session| {
        if out.is_null() {
            return Err(fail(RssnStatus::NullArgument, "null output pointer"));
        }
        let name = with_cell(session, cell_id, |value| {
            match unwrap_dimensioned(value) {
                | Value::Record(record) => {
                    record
                        .fields
                        .get(index)
                        .map(|(name, _)| name.clone())
                        .ok_or_else(|| fail(RssnStatus::NoValue, "field index out of range"))
                },
                | _ => Err(fail(RssnStatus::NoValue, "the value is not a record")),
            }
        })?;
        let text = to_c(name);
        if text.is_null() {
            return Err(fail(
                RssnStatus::Panic,
                "the field name is not representable as C text",
            ));
        }
        // SAFETY: `out` is non-null by the check above.
        unsafe { *out = text };
        Ok(())
    })
}

/// Writes the tag of a variant value to `*out` as a caller-owned string; free
/// it with `rssn_string_free`.
///
/// The payload is the value's elements: [`rssn_value_len`] and
/// [`rssn_value_at`] read them.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_tag(
    session: *const RssnSession,
    cell_id: u32,
    out: *mut *mut c_char,
) -> RssnStatus {
    with_session(session, |session| {
        if out.is_null() {
            return Err(fail(RssnStatus::NullArgument, "null output pointer"));
        }
        let tag = with_cell(session, cell_id, |value| {
            match unwrap_dimensioned(value) {
                | Value::Variant(variant) => Ok(variant.tag.clone()),
                | _ => Err(fail(RssnStatus::NoValue, "the value is not a variant")),
            }
        })?;
        let text = to_c(tag);
        if text.is_null() {
            return Err(fail(
                RssnStatus::Panic,
                "the variant tag is not representable as C text",
            ));
        }
        // SAFETY: `out` is non-null by the check above.
        unsafe { *out = text };
        Ok(())
    })
}

/// Releases a cell. A live session releases the rest with itself.
///
/// # Safety
/// `session` must be a live session.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_free(
    session: *const RssnSession,
    cell_id: u32,
) {
    let _ = guard(|| {
        // SAFETY: caller contract.
        match unsafe { borrow(session) } {
            | Ok(session) => {
                session.values.borrow_mut().remove(cell_id);
                RssnStatus::Ok
            },
            | Err(status) => status,
        }
    });
}

/// The number of live cells in `session`; a diagnostic for leak hunting.
///
/// # Safety
/// `session` must be a live session.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_value_count(session: *const RssnSession) -> usize {
    // SAFETY: caller contract.
    unsafe { borrow(session) }
        .map(|session| session.values.borrow().len())
        .unwrap_or(0)
}
