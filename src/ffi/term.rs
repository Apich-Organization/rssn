//! Sessions and terms.

use std::os::raw::c_char;

use super::RSSN_NO_TERM;
use super::RssnStatus;
use super::compute_error;
use super::fail;
use super::guard;
use super::text;
use super::to_c;
use crate::api::Session;
use crate::api::Term;

/// An rssn session: owns every term made through it. Opaque.
pub struct RssnSession(pub(crate) Session);

/// A dense numeric array extracted from a term, row-major.
#[repr(C)]
#[derive(Debug)]
pub struct RssnTensor {
    /// Number of dimensions; 0 for a scalar.
    pub rank: usize,
    /// `rank` extents.
    pub shape: *mut usize,
    /// Number of entries (the product of the extents).
    pub len: usize,
    /// `len` values.
    pub data: *mut f64,
}

/// Creates a session with every standard rule set. Free it with
/// `rssn_session_free`. Returns null if the built-in rules fail to install.
#[unsafe(no_mangle)]
pub extern "C" fn rssn_session_new() -> *mut RssnSession {
    std::panic::catch_unwind(|| Box::into_raw(Box::new(RssnSession(Session::new()))))
        .unwrap_or(std::ptr::null_mut())
}

/// Releases a session and all its terms. Null is ignored.
///
/// # Safety
/// `session` must be null or come from `rssn_session_new` and not be used
/// afterwards.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_session_free(session: *mut RssnSession) {
    if !session.is_null() {
        // SAFETY: allocated by `rssn_session_new`.
        drop(unsafe { Box::from_raw(session) });
    }
}

/// Borrows the session behind a pointer.
///
/// # Safety
/// `session` must be null or a live session pointer.
pub(crate) unsafe fn borrow_session<'a>(session: *const RssnSession) -> Result<&'a Session, RssnStatus> {
    // SAFETY: caller contract.
    unsafe { session.as_ref() }
        .map(|s| &s.0)
        .ok_or_else(|| fail(RssnStatus::NullArgument, "null session"))
}

/// Resolves a term handle.
pub(crate) fn resolve(
    session: &Session,
    id: u32,
) -> Result<Term<'_>, RssnStatus> {
    session
        .term_by_id(id)
        .ok_or_else(|| fail(RssnStatus::NoSuchTerm, "term handle does not belong to this session"))
}

/// Builds a term, mapping any failure to [`RSSN_NO_TERM`].
fn make(
    raw: *const RssnSession,
    build: impl FnOnce(&Session) -> Result<u32, RssnStatus>,
) -> u32 {
    let mut id = RSSN_NO_TERM;
    let _ = guard(|| {
        // SAFETY: caller contract of every constructor.
        match unsafe { borrow_session(raw) }.and_then(build) {
            | Ok(made) => {
                id = made;
                RssnStatus::Ok
            },
            | Err(status) => status,
        }
    });
    id
}

/// The symbol `name`, or `RSSN_NO_TERM`.
///
/// # Safety
/// `session` must be a live session; `name` a NUL-terminated string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_term_sym(
    session: *const RssnSession,
    name: *const c_char,
) -> u32 {
    // SAFETY: caller contract.
    make(session, |s| Ok(s.sym(unsafe { text(name) }?).id()))
}

/// An exact integer.
///
/// # Safety
/// `session` must be a live session.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_term_int(
    session: *const RssnSession,
    value: i64,
) -> u32 {
    make(session, |s| Ok(s.int(value).id()))
}

/// The exact fraction `numerator / denominator`, or `RSSN_NO_TERM` for a
/// zero denominator.
///
/// # Safety
/// `session` must be a live session.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_term_rational(
    session: *const RssnSession,
    numerator: i64,
    denominator: i64,
) -> u32 {
    make(session, |s| {
        s.rational(numerator, denominator)
            .map(Term::id)
            .ok_or_else(|| fail(RssnStatus::NoValue, "zero denominator"))
    })
}

/// A floating point number.
///
/// # Safety
/// `session` must be a live session.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_term_float(
    session: *const RssnSession,
    value: f64,
) -> u32 {
    make(session, |s| Ok(s.float(value).id()))
}

/// Parses infix text such as `"diff(sin(x)^2, x)"`.
///
/// # Safety
/// `session` must be a live session; `source` a NUL-terminated string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_term_parse(
    session: *const RssnSession,
    source: *const c_char,
) -> u32 {
    make(session, |s| {
        // SAFETY: caller contract.
        let source = unsafe { text(source) }?;
        s.parse(source).map(Term::id).map_err(|e| compute_error(&e))
    })
}

/// Applies the operator named `op` to `count` argument handles.
///
/// # Safety
/// `session` must be a live session, `op` a NUL-terminated string and
/// `args` point to `count` handles (it may be null when `count` is 0).
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_term_apply(
    session: *const RssnSession,
    op: *const c_char,
    args: *const u32,
    count: usize,
) -> u32 {
    make(session, |s| {
        // SAFETY: caller contract.
        let op = unsafe { text(op) }?;
        let ids: &[u32] = match (args.is_null(), count) {
            | (_, 0) => &[],
            | (true, _) => return Err(fail(RssnStatus::NullArgument, "null argument array")),
            // SAFETY: `args` points to `count` handles by contract.
            | (false, n) => unsafe { std::slice::from_raw_parts(args, n) },
        };
        let terms = ids.iter().map(|&id| resolve(s, id)).collect::<Result<Vec<_>, _>>()?;
        s.call(op, &terms).map(Term::id).map_err(|e| compute_error(&e))
    })
}

/// Stores the value of a literal number term in `*out`.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_term_to_f64(
    session: *const RssnSession,
    term: u32,
    out: *mut f64,
) -> RssnStatus {
    guard(|| {
        // SAFETY: caller contract.
        let s = match unsafe { borrow_session(session) } {
            | Ok(s) => s,
            | Err(e) => return e,
        };
        let value = match resolve(s, term) {
            | Ok(t) => t.as_f64(),
            | Err(e) => return e,
        };
        match (value, out.is_null()) {
            | (_, true) => fail(RssnStatus::NullArgument, "null output pointer"),
            | (None, _) => fail(RssnStatus::NoValue, "term is not a literal number"),
            | (Some(v), false) => {
                // SAFETY: non-null and valid by contract.
                unsafe { *out = v };
                RssnStatus::Ok
            },
        }
    })
}

/// Evaluates a term numerically with `count` symbol bindings
/// `names[i] = values[i]`, without running any rules.
///
/// # Safety
/// `session` must be a live session; `names` and `values` must hold
/// `count` entries (null when `count` is 0); `out` must be valid.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_term_eval(
    session: *const RssnSession,
    term: u32,
    names: *const *const c_char,
    values: *const f64,
    count: usize,
    out: *mut f64,
) -> RssnStatus {
    guard(|| {
        if out.is_null() || (count > 0 && (names.is_null() || values.is_null())) {
            return fail(RssnStatus::NullArgument, "null pointer argument");
        }
        // SAFETY: caller contract.
        let s = match unsafe { borrow_session(session) } {
            | Ok(s) => s,
            | Err(e) => return e,
        };
        let t = match resolve(s, term) {
            | Ok(t) => t,
            | Err(e) => return e,
        };
        let mut bindings = Vec::with_capacity(count);
        for i in 0..count {
            // SAFETY: both arrays hold `count` entries by contract.
            let (name, value) = unsafe { (*names.add(i), *values.add(i)) };
            // SAFETY: each name is a NUL-terminated string by contract.
            match unsafe { text(name) } {
                | Ok(name) => bindings.push((name, value)),
                | Err(e) => return e,
            }
        }
        match t.eval(&bindings) {
            | Some(v) => {
                // SAFETY: non-null and valid by contract.
                unsafe { *out = v };
                RssnStatus::Ok
            },
            | None => fail(RssnStatus::NotNumeric, "term does not evaluate to a number"),
        }
    })
}

/// Renders a term with `render`; null on failure.
fn render(
    raw: *const RssnSession,
    term: u32,
    render: fn(Term<'_>) -> String,
) -> *mut c_char {
    let mut out = std::ptr::null_mut();
    let _ = guard(|| {
        // SAFETY: caller contract of every renderer.
        match unsafe { borrow_session(raw) }.and_then(|s| resolve(s, term)) {
            | Ok(t) => {
                out = to_c(render(t));
                RssnStatus::Ok
            },
            | Err(e) => e,
        }
    });
    out
}

/// The term as infix text (parseable by `rssn_term_parse`); free with
/// `rssn_string_free`. Null on failure.
///
/// # Safety
/// `session` must be a live session.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_term_to_string(
    session: *const RssnSession,
    term: u32,
) -> *mut c_char {
    render(session, term, |t| t.to_string())
}

/// The term as LaTeX math; free with `rssn_string_free`.
///
/// # Safety
/// `session` must be a live session.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_term_to_latex(
    session: *const RssnSession,
    term: u32,
) -> *mut c_char {
    render(session, term, |t| t.to_latex())
}

/// The term as Typst math; free with `rssn_string_free`.
///
/// # Safety
/// `session` must be a live session.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_term_to_typst(
    session: *const RssnSession,
    term: u32,
) -> *mut c_char {
    render(session, term, |t| t.to_typst())
}

/// The term as a multi-line Unicode drawing; free with `rssn_string_free`.
///
/// # Safety
/// `session` must be a live session.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_term_to_pretty(
    session: *const RssnSession,
    term: u32,
) -> *mut c_char {
    render(session, term, |t| t.to_pretty())
}

/// Extracts a (possibly nested) list of numbers as a dense array into
/// `*out`; release it with `rssn_tensor_free`.
///
/// # Safety
/// `session` must be a live session and `out` a valid pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_term_tensor_data(
    session: *const RssnSession,
    term: u32,
    out: *mut RssnTensor,
) -> RssnStatus {
    guard(|| {
        if out.is_null() {
            return fail(RssnStatus::NullArgument, "null output pointer");
        }
        // SAFETY: caller contract.
        let tensor = match unsafe { borrow_session(session) }.and_then(|s| resolve(s, term)) {
            | Ok(t) => t.as_tensor(),
            | Err(e) => return e,
        };
        let Some((shape, data)) = tensor else {
            return fail(RssnStatus::NoValue, "term is not a dense array of numbers");
        };
        let (rank, len) = (shape.len(), data.len());
        let shape = Box::into_raw(shape.into_boxed_slice()).cast::<usize>();
        let data = Box::into_raw(data.into_boxed_slice()).cast::<f64>();
        // SAFETY: non-null and valid by contract.
        unsafe { *out = RssnTensor { rank, shape, len, data } };
        RssnStatus::Ok
    })
}

/// Releases the buffers of a tensor filled by `rssn_term_tensor_data` and
/// zeroes it. Null is ignored.
///
/// # Safety
/// `tensor` must be null or filled by `rssn_term_tensor_data` and not yet
/// freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rssn_tensor_free(tensor: *mut RssnTensor) {
    // SAFETY: caller contract.
    let Some(t) = (unsafe { tensor.as_mut() }) else {
        return;
    };
    if !t.shape.is_null() {
        // SAFETY: allocated as a boxed slice of `rank` extents.
        drop(unsafe { Box::from_raw(std::ptr::slice_from_raw_parts_mut(t.shape, t.rank)) });
    }
    if !t.data.is_null() {
        // SAFETY: allocated as a boxed slice of `len` values.
        drop(unsafe { Box::from_raw(std::ptr::slice_from_raw_parts_mut(t.data, t.len)) });
    }
    *t = RssnTensor {
        rank: 0,
        shape: std::ptr::null_mut(),
        len: 0,
        data: std::ptr::null_mut(),
    };
}
