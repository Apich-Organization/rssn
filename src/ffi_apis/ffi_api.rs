//! FFI API for the rssn library.
//!
//! This module provides a C-compatible foreign function interface (FFI) for interacting
//! with the core data structures and functions of the `rssn` library.
//!
//! The primary design pattern used here is the "handle-based" interface. Instead of
//! exposing complex Rust structs directly, which is unsafe and unstable, we expose
#![allow(unsafe_code)]
#![allow(clippy::indexing_slicing)]
#![allow(
    clippy::no_mangle_with_rust_abi
)]
#![allow(deprecated)]

use std::cell::RefCell;
use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;
use std::ptr;
use std::sync::Mutex;

use crate::plugins::manager::PluginManager;

thread_local! {
    static LAST_ERROR : RefCell < Option < CString >> = const { RefCell::new(None) };
}

/// A private helper to update the last error message for the current thread.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub(crate) unsafe fn update_last_error(
    err: String
) {

    let c_string = CString::new(err)
        .unwrap_or_else(|_| {

            CString::new(
                "Failed to create \
                 error message",
            )
            .expect(
                "Fallback error \
                 string should never \
                 fail to be created",
            )
        });

    LAST_ERROR.with(|prev| {

        *prev.borrow_mut() =
            Some(c_string);
    });
}

/// Retrieves the last error message set by an FFI function on the current thread.
///
/// The returned pointer is valid until the next call to an FFI function on the same thread.
/// The caller should not free this pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_get_last_error()
-> *const c_char {

    LAST_ERROR.with(|prev| {

        match *prev.borrow() {
            | Some(ref err) => {
                err.as_ptr()
            },
            | None => ptr::null(),
        }
    })
}

use std::sync::Arc;

use serde::Deserialize;
use serde::Serialize;

use crate::symbolic::core::Expr;

/// A macro to generate the boilerplate for a handle-based FFI API.
///
/// This macro creates three functions for a given type `$T`:
/// - `_from_json`: Deserializes a JSON string into a `$T` object and returns it as a raw pointer (handle).
/// - `_to_json`: Serializes a `$T` object (given by its handle) into a JSON string.
/// - `_free`: Frees the memory of a `$T` object given by its handle.

macro_rules! impl_handle_api {
    ($T:ty, $from_json:ident, $to_json:ident, $free:ident) => {
        /// Deserializes a JSON string into a `$T` object and returns a raw pointer (handle) to it.
        ///
        /// The caller is responsible for freeing the returned handle using `$free`.
        #[unsafe(no_mangle)]
        #[deprecated(
            since = "0.1.19",
            note = "This module is deprecated \
                    and please use special \
                    ffi api modules instead."
        )]

        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        pub unsafe extern "C" fn $from_json(json_ptr : *const c_char) -> *mut $T {

            if json_ptr.is_null() {

                return ptr::null_mut();
            }

            let json_str = match unsafe {

                CStr::from_ptr(json_ptr).to_str()
            } {
                | Ok(s) => s,
                | Err(_) => return ptr::null_mut(),
            };

            let result : Result<$T, _> = serde_json::from_str(json_str);

            match result {
                | Ok(obj) => Arc::into_raw(Arc::new(obj)).cast_mut(),
                | Err(_) => ptr::null_mut(),
            }
        }

        /// Serializes a `$T` object (given by its handle) into a JSON string.
        ///
        /// The caller is responsible for freeing the returned string using `free_string`.
        #[unsafe(no_mangle)]
        #[deprecated(
            since = "0.1.19",
            note = "This module is deprecated \
                    and please use special \
                    ffi api modules instead."
        )]

        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        pub unsafe extern "C" fn $to_json(handle : *mut $T) -> *mut c_char {

            if handle.is_null() {

                return ptr::null_mut();
            }

            let obj = unsafe {

                &*handle
            };

            match serde_json::to_string(obj) {
                | Ok(json_str) => {
                    match CString::new(json_str) {
                        | Ok(c_str) => c_str.into_raw(),
                        | Err(_) => ptr::null_mut(),
                    }
                },
                | Err(_) => ptr::null_mut(),
            }
        }

        /// Frees the memory associated with a `$T` handle.
        #[deprecated(
            since = "0.1.6",
            note = "Please use the handle-based `rssn_expr_free` instead."
        )]
        #[unsafe(no_mangle)]

        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        pub unsafe extern "C" fn $free(handle : *mut $T) {

            if !handle.is_null() {

                unsafe {

                    let _ = Arc::from_raw(handle);
                }
            }
        }
    };
}

impl_handle_api!(
    Expr,
    expr_from_json,
    expr_to_json,
    expr_free
);

/// Returns the string representation of an `Expr` handle.
///
/// The caller is responsible for freeing the returned string using `free_string`.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn expr_to_string(
    handle: *mut Expr
) -> *mut c_char {

    if handle.is_null() {

        return ptr::null_mut();
    }

    let expr = unsafe {

        &*handle
    };

    match CString::new(expr.to_string())
    {
        | Ok(s) => s.into_raw(),
        | Err(_) => ptr::null_mut(),
    }
}

use crate::symbolic::handles::HANDLE_MANAGER;
use crate::symbolic::simplify_dag::simplify;
use crate::symbolic::unit_unification::unify_expression;

/// Creates an expression from a JSON string and returns a thread-safe handle.
///
/// Returns 0 if the JSON is invalid.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_expr_create(
    json_ptr: *const c_char
) -> usize {

    unsafe {

        if json_ptr.is_null() {

            update_last_error(
                "Null pointer passed \
                 to rssn_expr_create"
                    .to_string(),
            );

            return 0;
        }

        let json_str =
            match CStr::from_ptr(
                json_ptr,
            )
            .to_str()
            {
                | Ok(s) => s,
                | Err(e) => {

                    update_last_error(
                        format!(
                "Invalid UTF-8 in \
                 input string: {e}"
            ),
                    );

                    return 0;
                },
            };

        match serde_json::from_str(
            json_str,
        ) {
            | Ok(expr) => {
                HANDLE_MANAGER
                    .insert(expr)
            },
            | Err(e) => {

                update_last_error(
                    format!(
                "JSON deserialization \
                 error: {e}"
            ),
                );

                0
            },
        }
    }
}

/// Frees the memory associated with an expression handle.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_expr_free(
    handle: usize
) {

    if handle != 0 {

        HANDLE_MANAGER.free(handle);
    }
}

/// Simplifies an expression handle and returns a handle to the new, simplified expression.
///
/// Returns 0 on error (e.g., invalid handle).
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_expr_simplify(
    handle: &usize
) -> usize {

    unsafe {

        match HANDLE_MANAGER
            .get(*handle)
        {
            | Some(expr) => {

                let simplified_expr =
                    simplify(
                        &(*expr)
                            .clone(),
                    );

                HANDLE_MANAGER.insert(
                    simplified_expr,
                )
            },
            | None => {

                update_last_error(format!(
                "Invalid handle passed to rssn_expr_simplify: {handle}"
            ));

                0
            },
        }
    }
}

#[derive(Serialize, Deserialize)]
/// A generic FFI-compatible result type, used for returning either a successful value or an error.

pub struct FfiResult<T, E> {
    /// The successful result value, if any.
    pub ok: Option<T>,
    /// The error message or error object, if any.
    pub err: Option<E>,
}

impl<T, E> FfiResult<T, E> {
    /// Creates a new `FfiResult` with a successful value.

    pub const fn ok(val: T) -> Self {

        Self {
            ok: Some(val),
            err: None,
        }
    }

    /// Creates a new `FfiResult` with an error.

    pub const fn err(err: E) -> Self {

        Self {
            ok: None,
            err: Some(err),
        }
    }
}

/// Simplifies an `Expr` and returns a handle to the new, simplified expression.
///
/// The caller is responsible for freeing the returned handle using `expr_free`.
#[deprecated(
    since = "0.1.6",
    note = "Please use the \
            handle-based \
            `rssn_expr_simplify` \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn expr_simplify(
    handle: *mut Expr
) -> *mut Expr {

    unsafe {

        // 1. Check the pointer itself (now *mut Expr)
        if handle.is_null() {

            return ptr::null_mut();
        }

        // 2. Dereference the raw pointer to get the actual Expr value.
        // The raw pointer is treated as a shared reference (&*handle)
        // to access the data without consuming it.
        let expr_ref = &*handle;

        // 3. Clone the *value* of the expression, and take a reference to the clone.
        // This gives `&Expr` required by `simplify`.
        let simplified_expr =
            simplify(&expr_ref.clone());

        // 4. Return the new expression as a raw pointer.
        Arc::into_raw(Arc::new(
            simplified_expr,
        ))
        .cast_mut()
    }
}

/// Attempts to unify the units within an expression.
///
/// Returns a JSON string representing a `FfiResult` which contains either the
/// new `Expr` object in the `ok` field or an error message in the `err` field.
/// The caller can then use `expr_from_json` to get a handle to the new expression.
/// The caller is responsible for freeing the returned string using `free_string`.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn expr_unify_expression(
    handle: *mut Expr
) -> *mut c_char {

    if handle.is_null() {

        let result = FfiResult {
            ok : None::<Expr>,
            err : Some("Null pointer passed to expr_unify_expression".to_string()),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let expr = unsafe {

        &*handle
    };

    let unification_result =
        unify_expression(expr);

    let ffi_result =
        match unification_result {
            | Ok(unified_expr) => {
                FfiResult {
                    ok: Some(
                        unified_expr,
                    ),
                    err: None,
                }
            },
            | Err(e) => {
                FfiResult {
                    ok: None,
                    err: Some(e),
                }
            },
        };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

/// Allocates and returns a test string ("pong") to the caller.
///
/// This function serves as a more advanced health check for the FFI interface.
/// It allows the client to verify two things:
/// 1. That the FFI function can be called successfully.
/// 2. That memory allocated in Rust can be safely passed to and then freed by the client
///    by calling `free_string` on the returned pointer.
///
/// Returns a pointer to a null-terminated C string. The caller is responsible for freeing this string.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_test_string_passing()
-> *mut c_char {

    match CString::new("pong") {
        | Ok(s) => s.into_raw(),
        | Err(_) => ptr::null_mut(),
    }
}

/// Frees a C string that was allocated by this library.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn free_string(
    s: *mut c_char
) {

    if !s.is_null() {

        unsafe {

            let _ =
                CString::from_raw(s);
        }
    }
}

use crate::output::latex::to_latex;
use crate::output::pretty_print::pretty_print;

/// Converts an expression to a LaTeX string.
///
/// The caller is responsible for freeing the returned string using `free_string`.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn expr_to_latex(
    handle: *mut Expr
) -> *mut c_char {

    if handle.is_null() {

        return ptr::null_mut();
    }

    let expr = unsafe {

        &*handle
    };

    let latex_str = to_latex(expr);

    match CString::new(latex_str) {
        | Ok(s) => s.into_raw(),
        | Err(_) => ptr::null_mut(),
    }
}

/// Converts an expression to a formatted, pretty-printed string.
///
/// The caller is responsible for freeing the returned string using `free_string`.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn expr_to_pretty_string(
    handle: *mut Expr
) -> *mut c_char {

    if handle.is_null() {

        return ptr::null_mut();
    }

    let expr = unsafe {

        &*handle
    };

    let pretty_str = pretty_print(expr);

    match CString::new(pretty_str) {
        | Ok(s) => s.into_raw(),
        | Err(_) => ptr::null_mut(),
    }
}

#[derive(Deserialize)]

struct LagrangeInput {
    points: Vec<(f64, f64)>,
}

#[derive(Serialize)]

struct FfiPolynomial {
    coeffs: Vec<f64>,
}

#[derive(Deserialize)]

struct BezierInput {
    control_points: Vec<Vec<f64>>,
    t: f64,
}

/// Computes a Lagrange interpolating polynomial and returns its coefficients as a JSON string.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn interpolate_lagrange(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        let result = FfiResult {
            ok : None::<FfiPolynomial>,
            err : Some("Null pointer passed to interpolate_lagrange".to_string()),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let input: Result<
        LagrangeInput,
        _,
    > = serde_json::from_str(json_str);

    let ffi_result = match input {
        | Ok(input_data) => {
            match interp_module::lagrange_interpolation(&input_data.points) {
                | Ok(poly) => {

                    let ffi_poly = FfiPolynomial {
                        coeffs : poly.coeffs,
                    };

                    FfiResult {
                        ok : Some(ffi_poly),
                        err : None::<String>,
                    }
                },
                | Err(e) => {
                    FfiResult {
                        ok : None,
                        err : Some(e),
                    }
                },
            }
        },
        | Err(e) => {
            FfiResult {
                ok : None,
                err : Some(format!(
                    "JSON deserialization error: {e}"
                )),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

/// Evaluates a point on a Bézier curve and returns the coordinates as a JSON string.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn interpolate_bezier_curve(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        let result = FfiResult {
            ok : None::<Vec<f64>>,
            err : Some("Null pointer passed to interpolate_bezier_curve".to_string()),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let input: Result<BezierInput, _> =
        serde_json::from_str(json_str);

    let ffi_result = match input {
        | Ok(input_data) => {

            let point = interp_module::bezier_curve(
                &input_data.control_points,
                input_data.t,
            );

            FfiResult {
                ok: Some(point),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(format!(
                "JSON deserialization \
                 error: {e}"
            )),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

use crate::numerical::combinatorics;
use crate::numerical::vector;

#[derive(Deserialize)]

struct VecInput {
    v: Vec<f64>,
}

#[derive(Deserialize)]

struct TwoVecInput {
    v1: Vec<f64>,
    v2: Vec<f64>,
}

#[derive(Deserialize)]

struct VecScalarInput {
    v: Vec<f64>,
    s: f64,
}

#[derive(Deserialize)]

struct U64Input {
    n: u64,
}

#[derive(Deserialize)]

struct TwoU64Input {
    n: u64,
    k: u64,
}

macro_rules! impl_ffi_1_vec_in_f64_out {
    ($fn_name:ident, $wrapped_fn:ident, $note:expr_2021) => {
        #[deprecated(since = "0.1.6", note = $note)]
        #[unsafe(no_mangle)]
        /// Implements a FFI function that takes a JSON string representing a vector (Vec<f64>),
        /// calls a wrapped function that operates on it, and returns the f64 result as a JSON string.
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        pub unsafe extern "C" fn $fn_name(json_ptr : *const c_char) -> *mut c_char {

            if json_ptr.is_null() {

                return match CString::new("{\"err\":\"Null pointer passed to function\"}") {
                    | Ok(s) => s.into_raw(),
                    | Err(_) => ptr::null_mut(),
                };
            }

            let json_str = match unsafe {

                CStr::from_ptr(json_ptr).to_str()
            } {
                | Ok(s) => s,
                | Err(_) => return ptr::null_mut(),
            };

            let input : Result<VecInput, _> = serde_json::from_str(json_str);

            let ffi_result = match input {
                | Ok(input_data) => {
                    FfiResult {
                        ok : Some(vector::$wrapped_fn(
                            &input_data.v,
                        )),
                        err : None::<String>,
                    }
                },
                | Err(e) => {
                    FfiResult {
                        ok : None,
                        err : Some(format!(
                            "JSON error: {}",
                            e
                        )),
                    }
                },
            };

            match serde_json::to_string(&ffi_result) {
                | Ok(json_str) => {
                    match CString::new(json_str) {
                        | Ok(c_str) => c_str.into_raw(),
                        | Err(_) => ptr::null_mut(),
                    }
                },
                | Err(_) => ptr::null_mut(),
            }
        }
    };
}

macro_rules! impl_ffi_2_vec_in_f64_out {
    ($fn_name:ident, $wrapped_fn:ident, $note:expr_2021) => {
        #[deprecated(since = "0.1.6", note = $note)]
        #[unsafe(no_mangle)]
        /// Implements a FFI function that takes a JSON string representing two vectors (Vec<f64>),
        /// calls a wrapped function that operates on them, and returns the f64 result as a JSON string.
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        pub unsafe extern "C" fn $fn_name(json_ptr : *const c_char) -> *mut c_char {

            if json_ptr.is_null() {

                return match CString::new("{\"err\":\"Null pointer passed to function\"}") {
                    | Ok(s) => s.into_raw(),
                    | Err(_) => ptr::null_mut(),
                };
            }

            let json_str = match unsafe {

                CStr::from_ptr(json_ptr).to_str()
            } {
                | Ok(s) => s,
                | Err(_) => return ptr::null_mut(),
            };

            let input : Result<TwoVecInput, _> = serde_json::from_str(json_str);

            let ffi_result = match input {
                | Ok(input_data) => {

                    match vector::$wrapped_fn(
                        &input_data.v1,
                        &input_data.v2,
                    ) {
                        | Ok(val) => {
                            FfiResult {
                                ok : Some(val),
                                err : None,
                            }
                        },
                        | Err(e) => {
                            FfiResult {
                                ok : None,
                                err : Some(e),
                            }
                        },
                    }
                },
                | Err(e) => {
                    FfiResult {
                        ok : None,
                        err : Some(format!(
                            "JSON error: {}",
                            e
                        )),
                    }
                },
            };

            match serde_json::to_string(&ffi_result) {
                | Ok(json_str) => {
                    match CString::new(json_str) {
                        | Ok(c_str) => c_str.into_raw(),
                        | Err(_) => ptr::null_mut(),
                    }
                },
                | Err(_) => ptr::null_mut(),
            }
        }
    };
}

macro_rules! impl_ffi_2_vec_in_vec_out {
    ($fn_name:ident, $wrapped_fn:ident, $note:expr_2021) => {
        #[deprecated(since = "0.1.6", note = $note)]
        #[unsafe(no_mangle)]
        /// Implements a FFI function that takes a JSON string representing two vectors (Vec<f64>),
        /// calls a wrapped function that operates on them, and returns the Vec<f64> result as a JSON string.
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        pub unsafe extern "C" fn $fn_name(json_ptr : *const c_char) -> *mut c_char {

            if json_ptr.is_null() {

                return match CString::new("{\"err\":\"Null pointer passed to function\"}") {
                    | Ok(s) => s.into_raw(),
                    | Err(_) => ptr::null_mut(),
                };
            }

            let json_str = match unsafe {

                CStr::from_ptr(json_ptr).to_str()
            } {
                | Ok(s) => s,
                | Err(_) => return ptr::null_mut(),
            };

            let input : Result<TwoVecInput, _> = serde_json::from_str(json_str);

            let ffi_result = match input {
                | Ok(input_data) => {

                    match vector::$wrapped_fn(
                        &input_data.v1,
                        &input_data.v2,
                    ) {
                        | Ok(val) => {
                            FfiResult {
                                ok : Some(val),
                                err : None,
                            }
                        },
                        | Err(e) => {
                            FfiResult {
                                ok : None,
                                err : Some(e),
                            }
                        },
                    }
                },
                | Err(e) => {
                    FfiResult {
                        ok : None,
                        err : Some(format!(
                            "JSON error: {}",
                            e
                        )),
                    }
                },
            };

            match serde_json::to_string(&ffi_result) {
                | Ok(json_str) => {
                    match CString::new(json_str) {
                        | Ok(c_str) => c_str.into_raw(),
                        | Err(_) => ptr::null_mut(),
                    }
                },
                | Err(_) => ptr::null_mut(),
            }
        }
    };
}

macro_rules! impl_ffi_1_u64_in_f64_out {
    ($fn_name:ident, $wrapped_fn:ident, $note:expr_2021) => {
        #[deprecated(since = "0.1.6", note = $note)]
        #[unsafe(no_mangle)]
        /// Implements a FFI function that takes a JSON string representing a u64,
        /// calls a wrapped function that operates on it, and returns the f64 result as a JSON string.
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        pub unsafe extern "C" fn $fn_name(json_ptr : *const c_char) -> *mut c_char {

            if json_ptr.is_null() {

                return match CString::new("{\"err\":\"Null pointer passed to function\"}") {
                    | Ok(s) => s.into_raw(),
                    | Err(_) => ptr::null_mut(),
                };
            }

            let json_str = match unsafe {

                CStr::from_ptr(json_ptr).to_str()
            } {
                | Ok(s) => s,
                | Err(_) => return ptr::null_mut(),
            };

            let input : Result<U64Input, _> = serde_json::from_str(json_str);

            let ffi_result = match input {
                | Ok(input_data) => {
                    FfiResult {
                        ok : Some(combinatorics::$wrapped_fn(input_data.n)),
                        err : None::<String>,
                    }
                },
                | Err(e) => {
                    FfiResult {
                        ok : None,
                        err : Some(format!(
                            "JSON error: {}",
                            e
                        )),
                    }
                },
            };

            match serde_json::to_string(&ffi_result) {
                | Ok(json_str) => {
                    match CString::new(json_str) {
                        | Ok(c_str) => c_str.into_raw(),
                        | Err(_) => ptr::null_mut(),
                    }
                },
                | Err(_) => ptr::null_mut(),
            }
        }
    };
}

macro_rules! impl_ffi_2_u64_in_f64_out {
    ($fn_name:ident, $wrapped_fn:ident, $note:expr_2021) => {
        #[deprecated(since = "0.1.6", note = $note)]
        #[unsafe(no_mangle)]
        /// Implements a FFI function that takes a JSON string representing two u64 values,
        /// calls a wrapped function that operates on them, and returns the f64 result as a JSON string.
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        pub unsafe extern "C" fn $fn_name(json_ptr : *const c_char) -> *mut c_char {

            if json_ptr.is_null() {

                return match CString::new("{\"err\":\"Null pointer passed to function\"}") {
                    | Ok(s) => s.into_raw(),
                    | Err(_) => ptr::null_mut(),
                };
            }

            let json_str = match unsafe {

                CStr::from_ptr(json_ptr).to_str()
            } {
                | Ok(s) => s,
                | Err(_) => return ptr::null_mut(),
            };

            let input : Result<TwoU64Input, _> = serde_json::from_str(json_str);

            let ffi_result = match input {
                | Ok(input_data) => {
                    FfiResult {
                        ok : Some(
                            combinatorics::$wrapped_fn(
                                input_data.n,
                                input_data.k,
                            ),
                        ),
                        err : None::<String>,
                    }
                },
                | Err(e) => {
                    FfiResult {
                        ok : None,
                        err : Some(format!(
                            "JSON error: {}",
                            e
                        )),
                    }
                },
            };

            match serde_json::to_string(&ffi_result) {
                | Ok(json_str) => {
                    match CString::new(json_str) {
                        | Ok(c_str) => c_str.into_raw(),
                        | Err(_) => ptr::null_mut(),
                    }
                },
                | Err(_) => ptr::null_mut(),
            }
        }
    };
}

impl_ffi_1_vec_in_f64_out!(
    vector_norm,
    norm,
    "Please use rssn_vec_norm instead."
);

impl_ffi_2_vec_in_f64_out!(
    vector_dot_product,
    dot_product,
    "Please use rssn_vec_dot_product \
     instead."
);

impl_ffi_2_vec_in_f64_out!(
    vector_distance,
    distance,
    "Please use rssn_vec_distance \
     instead."
);

impl_ffi_2_vec_in_f64_out!(
    vector_angle,
    angle,
    "Please use rssn_vec_angle \
     instead."
);

impl_ffi_2_vec_in_vec_out!(
    vector_add,
    vec_add,
    "Please use rssn_vec_add instead."
);

impl_ffi_2_vec_in_vec_out!(
    vector_sub,
    vec_sub,
    "Please use rssn_vec_sub instead."
);

impl_ffi_2_vec_in_vec_out!(
    vector_cross_product,
    cross_product,
    "Please use \
     rssn_vec_cross_product instead."
);

#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_vec_scalar_mul \
            instead."
)]
#[unsafe(no_mangle)]

/// Multiplies a vector by a scalar.
///
/// This function is deprecated.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn vector_scalar_mul(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        return match CString::new(
            "{\"err\":\"Null pointer \
             passed to function\"}",
        ) {
            | Ok(s) => s.into_raw(),
            | Err(_) => ptr::null_mut(),
        };
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let input: Result<
        VecScalarInput,
        _,
    > = serde_json::from_str(json_str);

    let ffi_result = match input {
        | Ok(input_data) => {
            FfiResult {
                ok: Some(
                    vector::scalar_mul(
                        &input_data.v,
                        input_data.s,
                    ),
                ),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(format!(
                    "JSON error: {e}"
                )),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

impl_ffi_1_u64_in_f64_out!(
    combinatorics_factorial,
    factorial,
    "Please use rssn_comb_factorial \
     instead."
);

impl_ffi_2_u64_in_f64_out!(
    combinatorics_permutations,
    permutations,
    "Please use \
     rssn_comb_permutations instead."
);

impl_ffi_2_u64_in_f64_out!(
    combinatorics_combinations,
    combinations,
    "Please use \
     rssn_comb_combinations instead."
);

/// Computes the L2 norm of a vector.
///
/// # Arguments
/// * `data` - A pointer to the first element of the vector.
/// * `len` - The length of the vector.
/// * `result` - A pointer to store the computed norm.
///
/// # Returns
/// 0 on success, -1 on error (e.g., null pointer).
///
/// # Safety
/// The `data` pointer must point to a valid array of `len` f64 elements.
/// The `result` pointer must point to a valid f64 location.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_vec_norm(
    data: *const f64,
    len: usize,
    result: *mut f64,
) -> i32 {

    unsafe {

        if data.is_null()
            || result.is_null()
        {

            update_last_error(
                "Null pointer passed \
                 to rssn_vec_norm"
                    .to_string(),
            );

            return -1;
        }

        let vec_slice = {

            std::slice::from_raw_parts(
                data, len,
            )
        };

        {

            *result =
                vector::norm(vec_slice);
        };

        0
    }
}

/// Computes the dot product of two vectors.
///
/// # Arguments
/// * `d1` - A pointer to the first element of the first vector.
/// * `l1` - The length of the first vector.
/// * `d2` - A pointer to the first element of the second vector.
/// * `l2` - The length of the second vector.
/// * `result` - A pointer to store the computed dot product.
///
/// # Returns
/// 0 on success, -1 on error (e.g., null pointer or mismatched lengths).
///
/// # Safety
/// The `d1` and `d2` pointers must point to valid arrays of `l1` and `l2` f64 elements respectively.
/// The `result` pointer must point to a valid f64 location.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_vec_dot_product(
    d1: *const f64,
    l1: usize,
    d2: *const f64,
    l2: usize,
    result: *mut f64,
) -> i32 {

    unsafe {

        if d1.is_null()
            || d2.is_null()
            || result.is_null()
        {

            update_last_error(
            "Null pointer passed to \
             rssn_vec_dot_product"
                .to_string(),
        );

            return -1;
        }

        let v1 = {

            std::slice::from_raw_parts(
                d1, l1,
            )
        };

        let v2 = {

            std::slice::from_raw_parts(
                d2, l2,
            )
        };

        match vector::dot_product(
            v1, v2,
        ) {
            | Ok(val) => {

                {

                    *result = val;
                };

                0
            },
            | Err(e) => {

                update_last_error(e);

                -1
            },
        }
    }
}

/// Computes the factorial of a number `n`.
///
/// # Arguments
/// * `n` - The number to compute the factorial of.
/// * `result` - A pointer to store the computed factorial.
///
/// # Returns
/// 0 on success, -1 on error (e.g., null pointer).
///
/// # Safety
/// The `result` pointer must point to a valid f64 location.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_comb_factorial(
    n: u64,
    result: *mut f64,
) -> i32 {

    unsafe {

        if result.is_null() {

            update_last_error(
            "Null pointer passed to \
             rssn_comb_factorial"
                .to_string(),
        );

            return -1;
        }

        {

            *result =
            combinatorics::factorial(n);
        };

        0
    }
}

/// Computes the number of permutations (nPk).
///
/// # Arguments
/// * `n` - The total number of items.
/// * `k` - The number of items to choose.
/// * `result` - A pointer to store the computed number of permutations.
///
/// # Returns
/// 0 on success, -1 on error (e.g., null pointer).
///
/// # Safety
/// The `result` pointer must point to a valid f64 location.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_comb_permutations(
    n: u64,
    k: u64,
    result: *mut f64,
) -> i32 {

    unsafe {

        if result.is_null() {

            update_last_error(
            "Null pointer passed to \
             rssn_comb_permutations"
                .to_string(),
        );

            return -1;
        }

        {

            *result =
            combinatorics::permutations(
                n, k,
            );
        };

        0
    }
}

/// Computes the number of combinations (nCk).
///
/// # Arguments
/// * `n` - The total number of items.
/// * `k` - The number of items to choose.
/// * `result` - A pointer to store the computed number of combinations.
///
/// # Returns
/// 0 on success, -1 on error (e.g., null pointer).
///
/// # Safety
/// The `result` pointer must point to a valid f64 location.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_comb_combinations(
    n: u64,
    k: u64,
    result: *mut f64,
) -> i32 {

    unsafe {

        if result.is_null() {

            update_last_error(
            "Null pointer passed to \
             rssn_comb_combinations"
                .to_string(),
        );

            return -1;
        }

        {

            *result =
            combinatorics::combinations(
                n, k,
            );
        };

        0
    }
}

use crate::numerical::number_theory as nt;

#[derive(Deserialize)]

struct TwoU64NtInput {
    a: u64,
    b: u64,
}

#[derive(Deserialize)]

struct ModPowInput {
    base: u64,
    exp: u64,
    modulus: u64,
}

#[derive(Deserialize)]

struct TwoI64NtInput {
    a: i64,
    b: i64,
}

#[derive(Deserialize)]

struct U64NtInput {
    n: u64,
}

macro_rules! impl_ffi_2_u64_in_u64_out {
    ($fn_name:ident, $wrapped_fn:ident, $note:expr_2021) => {
        /// Computes `$wrapped_fn` from two `u64` values passed as a JSON string.
        ///
        /// This function is deprecated.
        #[deprecated(since = "0.1.6", note = $note)]
        #[unsafe(no_mangle)]
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        pub unsafe extern "C" fn $fn_name(json_ptr: *const c_char) -> *mut c_char {

            if json_ptr.is_null() {

                return match CString::new("{\"err\":\"Null pointer passed to function\"}") {
                    | Ok(s) => s.into_raw(),
                    | Err(_) => ptr::null_mut(),
                };
            }

            let json_str = match unsafe {

                CStr::from_ptr(json_ptr).to_str()
            } {
                | Ok(s) => s,
                | Err(_) => return ptr::null_mut(),
            };

            let input : Result<TwoU64NtInput, _> = serde_json::from_str(json_str);

            let ffi_result = match input {
                | Ok(input_data) => {
                    FfiResult {
                        ok : Some(nt::$wrapped_fn(
                            input_data.a,
                            input_data.b,
                        )),
                        err : None::<String>,
                    }
                },
                | Err(e) => {
                    FfiResult {
                        ok : None,
                        err : Some(format!(
                            "JSON error: {}",
                            e
                        )),
                    }
                },
            };

            match serde_json::to_string(&ffi_result) {
                | Ok(json_str) => {
                    match CString::new(json_str) {
                        | Ok(c_str) => c_str.into_raw(),
                        | Err(_) => ptr::null_mut(),
                    }
                },
                | Err(_) => ptr::null_mut(),
            }
        }
    };
}

macro_rules! impl_ffi_1_u64_in_bool_out {
    ($fn_name:ident, $wrapped_fn:ident, $note:expr_2021) => {
        /// Computes `$wrapped_fn` from a single `u64` value passed as a JSON string and returns a boolean.
        ///
        /// This function is deprecated.
        #[deprecated(since = "0.1.6", note = $note)]
        #[unsafe(no_mangle)]
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        pub unsafe extern "C" fn $fn_name(json_ptr: *const c_char) -> *mut c_char {

            if json_ptr.is_null() {

                return match CString::new("{\"err\":\"Null pointer passed to function\"}") {
                    | Ok(s) => s.into_raw(),
                    | Err(_) => ptr::null_mut(),
                };
            }

            let json_str = match unsafe {

                CStr::from_ptr(json_ptr).to_str()
            } {
                | Ok(s) => s,
                | Err(_) => return ptr::null_mut(),
            };

            let input : Result<U64NtInput, _> = serde_json::from_str(json_str);

            let ffi_result = match input {
                | Ok(input_data) => {
                    FfiResult {
                        ok : Some(nt::$wrapped_fn(
                            input_data.n,
                        )),
                        err : None::<String>,
                    }
                },
                | Err(e) => {
                    FfiResult {
                        ok : None,
                        err : Some(format!(
                            "JSON error: {}",
                            e
                        )),
                    }
                },
            };

            match serde_json::to_string(&ffi_result) {
                | Ok(json_str) => {
                    match CString::new(json_str) {
                        | Ok(c_str) => c_str.into_raw(),
                        | Err(_) => ptr::null_mut(),
                    }
                },
                | Err(_) => ptr::null_mut(),
            }
        }
    };
}

impl_ffi_2_u64_in_u64_out!(
    nt_gcd,
    gcd,
    "Please use rssn_nt_gcd instead."
);

impl_ffi_1_u64_in_bool_out!(
    nt_is_prime_miller_rabin,
    is_prime_miller_rabin,
    "Please use rssn_nt_is_prime \
     instead."
);

/// Computes the modular exponentiation (base^exp % modulus) using input from a JSON string.
///
/// This function is deprecated.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_nt_mod_pow instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn nt_mod_pow(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        return match CString::new(
            "{\"err\":\"Null pointer \
             passed to function\"}",
        ) {
            | Ok(s) => s.into_raw(),
            | Err(_) => ptr::null_mut(),
        };
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let input: Result<ModPowInput, _> =
        serde_json::from_str(json_str);

    let ffi_result = match input {
        | Ok(d) => {
            FfiResult {
                ok: Some(nt::mod_pow(
                    u128::from(d.base),
                    d.exp,
                    d.modulus,
                )),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(format!(
                    "JSON error: {e}"
                )),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

/// Computes the modular inverse (a^-1 mod b) using input from a JSON string.
///
/// This function is deprecated.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_nt_mod_inverse \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn nt_mod_inverse(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        return match CString::new(
            "{\"err\":\"Null pointer \
             passed to function\"}",
        ) {
            | Ok(s) => s.into_raw(),
            | Err(_) => ptr::null_mut(),
        };
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let input: Result<
        TwoI64NtInput,
        _,
    > = serde_json::from_str(json_str);

    let ffi_result = match input {
        | Ok(d) => {
            FfiResult {
                ok: nt::mod_inverse(
                    d.a, d.b,
                ),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(format!(
                    "JSON error: {e}"
                )),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

use crate::numerical::special::{
    self as special_module,
};

#[derive(Deserialize)]

struct SpecialFunc1Input {
    x: f64,
}

#[derive(Deserialize)]

struct SpecialFunc2Input {
    a: f64,
    b: f64,
}

macro_rules! impl_special_fn_one_arg {
    ($fn_name:ident, $wrapped_fn:ident, $note:expr_2021) => {
        #[deprecated(since = "0.1.6", note = $note)]
        #[unsafe(no_mangle)]
        /// Implements a FFI function that takes a JSON string representing a single f64 argument,
        /// calls a wrapped special function that operates on it, and returns the f64 result as a JSON string.
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        pub unsafe extern "C" fn $fn_name(json_ptr : *const c_char) -> *mut c_char {

            if json_ptr.is_null() {

                let result = FfiResult {
                    ok : None::<f64>,
                    err : Some("Null pointer passed to function".to_string()),
                };

                return match serde_json::to_string(&result) {
                    | Ok(json_str) => {
                        match CString::new(json_str) {
                            | Ok(c_str) => c_str.into_raw(),
                            | Err(_) => ptr::null_mut(),
                        }
                    },
                    | Err(_) => ptr::null_mut(),
                };
            }

            let json_str = match unsafe {

                CStr::from_ptr(json_ptr).to_str()
            } {
                | Ok(s) => s,
                | Err(_) => return ptr::null_mut(),
            };

            let input : Result<SpecialFunc1Input, _> = serde_json::from_str(json_str);

            let ffi_result = match input {
                | Ok(input_data) => {

                    let result = special_module::$wrapped_fn(input_data.x);

                    FfiResult {
                        ok : Some(result),
                        err : None::<String>,
                    }
                },
                | Err(e) => {
                    FfiResult {
                        ok : None,
                        err : Some(format!(
                            "JSON deserialization error: {}",
                            e
                        )),
                    }
                },
            };

            match serde_json::to_string(&ffi_result) {
                | Ok(json_str) => {
                    match CString::new(json_str) {
                        | Ok(c_str) => c_str.into_raw(),
                        | Err(_) => ptr::null_mut(),
                    }
                },
                | Err(_) => ptr::null_mut(),
            }
        }
    };
}

macro_rules! impl_special_fn_two_args {
    ($fn_name:ident, $wrapped_fn:ident, $note:expr_2021) => {
        #[deprecated(since = "0.1.6", note = $note)]
        #[unsafe(no_mangle)]
        /// Implements a FFI function that takes a JSON string representing two f64 arguments,
        /// calls a wrapped special function that operates on them, and returns the f64 result as a JSON string.
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        pub unsafe extern "C" fn $fn_name(json_ptr : *const c_char) -> *mut c_char {

            if json_ptr.is_null() {

                let result = FfiResult {
                    ok : None::<f64>,
                    err : Some("Null pointer passed to function".to_string()),
                };

                return match serde_json::to_string(&result) {
                    | Ok(json_str) => {
                        match CString::new(json_str) {
                            | Ok(c_str) => c_str.into_raw(),
                            | Err(_) => ptr::null_mut(),
                        }
                    },
                    | Err(_) => ptr::null_mut(),
                };
            }

            let json_str = match unsafe {

                CStr::from_ptr(json_ptr).to_str()
            } {
                | Ok(s) => s,
                | Err(_) => return ptr::null_mut(),
            };

            let input : Result<SpecialFunc2Input, _> = serde_json::from_str(json_str);

            let ffi_result = match input {
                | Ok(input_data) => {

                    let result = special_module::$wrapped_fn(
                        input_data.a,
                        input_data.b,
                    );

                    FfiResult {
                        ok : Some(result),
                        err : None::<String>,
                    }
                },
                | Err(e) => {
                    FfiResult {
                        ok : None,
                        err : Some(format!(
                            "JSON deserialization error: {}",
                            e
                        )),
                    }
                },
            };

            match serde_json::to_string(&ffi_result) {
                | Ok(json_str) => {
                    match CString::new(json_str) {
                        | Ok(c_str) => c_str.into_raw(),
                        | Err(_) => ptr::null_mut(),
                    }
                },
                | Err(_) => ptr::null_mut(),
            }
        }
    };
}

impl_special_fn_one_arg!(
    special_gamma,
    gamma_numerical,
    "Please use rssn_spec_gamma \
     instead."
);

impl_special_fn_one_arg!(
    special_ln_gamma,
    ln_gamma_numerical,
    "Please use rssn_spec_ln_gamma \
     instead."
);

impl_special_fn_one_arg!(
    special_erf,
    erf_numerical,
    "Please use rssn_spec_erf instead."
);

impl_special_fn_one_arg!(
    special_erfc,
    erfc_numerical,
    "Please use rssn_spec_erfc \
     instead."
);

impl_special_fn_two_args!(
    special_beta,
    beta_numerical,
    "Please use rssn_spec_beta \
     instead."
);

impl_special_fn_two_args!(
    special_ln_beta,
    ln_beta_numerical,
    "Please use rssn_spec_ln_beta \
     instead."
);

macro_rules! impl_rssn_special_fn_one_arg {
    ($fn_name:ident, $wrapped_fn:ident) => {
        /// Computes the value of the `$wrapped_fn` special function for a single `f64` argument.
        ///
        /// # Arguments
        /// * `x` - The input value for the special function.
        /// * `result` - A pointer to store the computed result.
        ///
        /// # Returns
        /// 0 on success, -1 on error (e.g., null pointer).
        ///
        /// # Safety
        /// The `result` pointer must point to a valid `f64` location.
        #[unsafe(no_mangle)]
        #[deprecated(
            since = "0.1.19",
            note = "This module is deprecated \
                    and please use special \
                    ffi api modules instead."
        )]
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        pub unsafe extern "C" fn $fn_name(
            x : f64,
            result : *mut f64,
        ) -> i32 { unsafe {

            if result.is_null() {

                update_last_error(format!(
                    "Null pointer passed for 'result' to {}",
                    stringify!($fn_name)
                ));

                return -1;
            }

             {

                *result = special_module::$wrapped_fn(x)
            };

            0
        }}
    };
}

macro_rules! impl_rssn_special_fn_two_args {
    ($fn_name:ident, $wrapped_fn:ident) => {
        #[unsafe(no_mangle)]
        #[deprecated(
            since = "0.1.19",
            note = "This module is deprecated \
                    and please use special \
                    ffi api modules instead."
        )]
        #[allow(clippy::too_long_first_doc_paragraph)]
        /// Implements a FFI function that takes two f64 arguments `a` and `b`,
        /// calls a wrapped special function that operates on them, and stores the f64 result
        /// in the `result` pointer. Returns 0 on success, -1 on error.
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        pub unsafe extern "C" fn $fn_name(
            a : f64,
            b : f64,
            result : *mut f64,
        ) -> i32 { unsafe {

            if result.is_null() {

                update_last_error(format!(
                    "Null pointer passed for 'result' to {}",
                    stringify!($fn_name)
                ));

                return -1;
            }

             {

                *result = special_module::$wrapped_fn(a, b)
            };

            0
        }}
    };
}

impl_rssn_special_fn_one_arg!(
    rssn_spec_gamma,
    gamma_numerical
);

impl_rssn_special_fn_one_arg!(
    rssn_spec_ln_gamma,
    ln_gamma_numerical
);

impl_rssn_special_fn_one_arg!(
    rssn_spec_erf,
    erf_numerical
);

impl_rssn_special_fn_one_arg!(
    rssn_spec_erfc,
    erfc_numerical
);

impl_rssn_special_fn_two_args!(
    rssn_spec_beta,
    beta_numerical
);

impl_rssn_special_fn_two_args!(
    rssn_spec_ln_beta,
    ln_beta_numerical
);

use num_complex::Complex;

use crate::numerical::transforms::fft;
use crate::numerical::transforms::ifft;

#[derive(Serialize, Deserialize)]

struct TransformsInput {
    data: Vec<Complex<f64>>,
}

/// Computes the Fast Fourier Transform (FFT) of a sequence of complex numbers.
#[deprecated(
    since = "0.1.6",
    note = "Please use rssn_fft \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn transforms_fft(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        let result = FfiResult {
            ok: None::<Vec<Complex<f64>>>,
            err: Some(
                "Null pointer passed \
                 to transforms_fft"
                    .to_string(),
            ),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let input: Result<
        TransformsInput,
        _,
    > = serde_json::from_str(json_str);

    let ffi_result = match input {
        | Ok(mut input_data) => {

            fft(&mut input_data.data);

            FfiResult {
                ok: Some(
                    input_data.data,
                ),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(format!(
                "JSON deserialization \
                 error: {e}"
            )),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

/// Computes the Inverse Fast Fourier Transform (IFFT) of a sequence of complex numbers.
#[deprecated(
    since = "0.1.6",
    note = "Please use rssn_ifft \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn transforms_ifft(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        let result = FfiResult {
            ok: None::<Vec<Complex<f64>>>,
            err: Some(
                "Null pointer passed \
                 to transforms_ifft"
                    .to_string(),
            ),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let input: Result<
        TransformsInput,
        _,
    > = serde_json::from_str(json_str);

    let ffi_result = match input {
        | Ok(mut input_data) => {

            ifft(&mut input_data.data);

            FfiResult {
                ok: Some(
                    input_data.data,
                ),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(format!(
                "JSON deserialization \
                 error: {e}"
            )),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

use crate::numerical::transforms::fft_slice;
use crate::numerical::transforms::ifft_slice;

/// Computes the Fast Fourier Transform (FFT) of a sequence of complex numbers in-place.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_fft(
    data: *mut Complex<f64>,
    len: usize,
) -> i32 {

    unsafe {

        if data.is_null() {

            update_last_error(
                "Null pointer passed \
                 to rssn_fft"
                    .to_string(),
            );

            return -1;
        }

        let complex_slice = {

            std::slice::from_raw_parts_mut(
            data, len,
        )
        };

        fft_slice(complex_slice);

        0
    }
}

/// Computes the Inverse Fast Fourier Transform (IFFT) of a sequence of complex numbers in-place.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_ifft(
    data: *mut Complex<f64>,
    len: usize,
) -> i32 {

    unsafe {

        if data.is_null() {

            update_last_error(
                "Null pointer passed \
                 to rssn_ifft"
                    .to_string(),
            );

            return -1;
        }

        let complex_slice = {

            std::slice::from_raw_parts_mut(
            data, len,
        )
        };

        ifft_slice(complex_slice);

        0
    }
}

#[derive(Deserialize)]

struct PolyInput {
    expr: Expr,
    var: String,
}

#[derive(Deserialize)]

struct PolyDivInput {
    n: Expr,
    d: Expr,
    var: String,
}

#[derive(Deserialize)]

struct PolyFromCoeffsInput {
    coeffs: Vec<Expr>,
    var: String,
}

/// Checks if an expression is a polynomial with respect to a given variable.
///
/// This function is deprecated.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_poly_is_polynomial \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn poly_is_polynomial(
    json_ptr: *const c_char
) -> bool {

    if json_ptr.is_null() {

        return false;
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => return false,
    };

    let input: Result<PolyInput, _> =
        serde_json::from_str(json_str);

    match input {
        | Ok(poly_input) => {
            poly_module::is_polynomial(
                &poly_input.expr,
                &poly_input.var,
            )
        },
        | Err(_) => false,
    }
}

/// Computes the degree of a polynomial expression with respect to a given variable.
///
/// This function is deprecated.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_poly_degree instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn poly_degree(
    json_ptr: *const c_char
) -> i64 {

    if json_ptr.is_null() {

        return -1;
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => return -1,
    };

    let input: Result<PolyInput, _> =
        serde_json::from_str(json_str);

    match input {
        | Ok(poly_input) => {
            poly_module::polynomial_degree(
                &poly_input.expr,
                &poly_input.var,
            )
        },
        | Err(_) => -1,
    }
}

/// Computes the leading coefficient of a polynomial expression with respect to a given variable.
///
/// This function is deprecated.
#[deprecated(
    since = "0.1.6",
    note = "Please use rssn_poly_leading_coefficient instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn poly_leading_coefficient(
    handle: *mut Expr,
    var_ptr: *const c_char,
) -> *mut Expr {

    if handle.is_null()
        || var_ptr.is_null()
    {

        return ptr::null_mut();
    }

    let expr = unsafe {

        &*handle
    };

    let var = match unsafe {

        CStr::from_ptr(var_ptr).to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let result_expr = poly_module::leading_coefficient(expr, var);

    Arc::into_raw(Arc::new(
        result_expr,
    ))
    .cast_mut()
}

/// Performs polynomial long division on two expressions with respect to a given variable.
///
/// This function is deprecated.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_poly_long_division \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn poly_long_division(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        let result = FfiResult {
            ok: None::<MatrixPair>,
            err: Some(
                "Null pointer passed \
                 to poly_long_division"
                    .to_string(),
            ),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let input: Result<PolyDivInput, _> =
        serde_json::from_str(json_str);

    let ffi_result = match input {
        | Ok(div_input) => {

            let (q, r) = poly_module::polynomial_long_division(
                &div_input.n,
                &div_input.d,
                &div_input.var,
            );

            FfiResult {
                ok: Some(MatrixPair {
                    p1: q,
                    p2: r,
                }),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(format!(
                "JSON deserialization \
                 error: {e}"
            )),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

/// Converts a polynomial expression to a vector of its coefficients.
///
/// This function is deprecated.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_poly_to_coeffs \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn poly_to_coeffs_vec(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        let result = FfiResult {
            ok: None::<Vec<Expr>>,
            err: Some(
                "Null pointer passed \
                 to poly_to_coeffs_vec"
                    .to_string(),
            ),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let input: Result<PolyInput, _> =
        serde_json::from_str(json_str);

    let ffi_result = match input {
        | Ok(poly_input) => {

            let coeffs = poly_module::to_polynomial_coeffs_vec(
                &poly_input.expr,
                &poly_input.var,
            );

            FfiResult {
                ok: Some(coeffs),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(format!(
                "JSON deserialization \
                 error: {e}"
            )),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

/// Creates a polynomial expression from a vector of its coefficients.
///
/// This function is deprecated.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_poly_from_coeffs \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn poly_from_coeffs_vec(
    json_ptr: *const c_char
) -> *mut Expr {

    if json_ptr.is_null() {

        return ptr::null_mut();
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let input: Result<
        PolyFromCoeffsInput,
        _,
    > = serde_json::from_str(json_str);

    match input {
        | Ok(input_data) => {

            let result_expr =
                from_coeffs_to_expr(
                    &input_data.coeffs,
                    &input_data.var,
                );

            Arc::into_raw(Arc::new(
                result_expr,
            ))
            .cast_mut()
        },
        | Err(_) => ptr::null_mut(),
    }
}

/// Checks if an expression is a polynomial with respect to a given variable.
///
/// # Arguments
/// * `expr_handle` - A handle to the expression to check.
/// * `var_ptr` - A pointer to a null-terminated C string representing the variable.
/// * `result` - A pointer to store the boolean result (true if it's a polynomial, false otherwise).
///
/// # Returns
/// 0 on success, -1 on error (e.g., null pointer, invalid handle).
///
/// # Safety
/// The `var_ptr` must point to a valid null-terminated C string.
/// The `result` pointer must point to a valid boolean location.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_poly_is_polynomial(
    expr_handle: usize,
    var_ptr: *const c_char,
    result: *mut bool,
) -> i32 {

    unsafe {

        if result.is_null()
            || var_ptr.is_null()
        {

            update_last_error(
            "Null pointer passed to \
             rssn_poly_is_polynomial"
                .to_string(),
        );

            return -1;
        }

        let var = match CStr::from_ptr(
            var_ptr,
        )
        .to_str()
        {
            | Ok(s) => s,
            | Err(e) => {

                update_last_error(
                    format!(
                "Invalid UTF-8 in \
                 var_ptr: {e}"
            ),
                );

                return -1;
            },
        };

        match HANDLE_MANAGER
            .get(expr_handle)
        {
            | Some(expr) => {

                {

                    *result = poly_module::is_polynomial(&expr, var);
                };

                0
            },
            | None => {

                update_last_error(format!(
                "Invalid handle passed to rssn_poly_is_polynomial: {expr_handle}"
            ));

                -1
            },
        }
    }
}

/// Computes the degree of a polynomial expression with respect to a given variable.
///
/// # Arguments
/// * `expr_handle` - A handle to the expression to check.
/// * `var_ptr` - A pointer to a null-terminated C string representing the variable.
/// * `result` - A pointer to store the degree of the polynomial.
///
/// # Returns
/// 0 on success, -1 on error (e.g., null pointer, invalid handle).
///
/// # Safety
/// The `var_ptr` must point to a valid null-terminated C string.
/// The `result` pointer must point to a valid `i64` location.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_poly_degree(
    expr_handle: usize,
    var_ptr: *const c_char,
    result: *mut i64,
) -> i32 {

    unsafe {

        if result.is_null()
            || var_ptr.is_null()
        {

            update_last_error(
                "Null pointer passed \
                 to rssn_poly_degree"
                    .to_string(),
            );

            return -1;
        }

        let var = match CStr::from_ptr(
            var_ptr,
        )
        .to_str()
        {
            | Ok(s) => s,
            | Err(e) => {

                update_last_error(
                    format!(
                "Invalid UTF-8 in \
                 var_ptr: {e}"
            ),
                );

                return -1;
            },
        };

        match HANDLE_MANAGER
            .get(expr_handle)
        {
            | Some(expr) => {

                {

                    *result = poly_module::polynomial_degree(&expr, var);
                };

                0
            },
            | None => {

                update_last_error(format!(
                "Invalid handle passed to rssn_poly_degree: {expr_handle}"
            ));

                -1
            },
        }
    }
}

/// Performs polynomial long division on two expressions with respect to a given variable.
///
/// # Arguments
/// * `n_handle` - A handle to the numerator expression.
/// * `d_handle` - A handle to the denominator expression.
/// * `var_ptr` - A pointer to a null-terminated C string representing the variable.
/// * `q_handle` - A pointer to store the handle of the quotient expression.
/// * `r_handle` - A pointer to store the handle of the remainder expression.
///
/// # Returns
/// 0 on success, -1 on error (e.g., null pointer, invalid handle).
///
/// # Safety
/// The `var_ptr` must point to a valid null-terminated C string.
/// The `q_handle` and `r_handle` pointers must point to valid `usize` locations.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_poly_long_division(
    n_handle: usize,
    d_handle: usize,
    var_ptr: *const c_char,
    q_handle: *mut usize,
    r_handle: *mut usize,
) -> i32 {

    unsafe {

        if q_handle.is_null()
            || r_handle.is_null()
            || var_ptr.is_null()
        {

            update_last_error(
            "Null pointer passed to \
             rssn_poly_long_division"
                .to_string(),
        );

            return -1;
        }

        let var = match CStr::from_ptr(
            var_ptr,
        )
        .to_str()
        {
            | Ok(s) => s,
            | Err(e) => {

                update_last_error(
                    format!(
                "Invalid UTF-8 in \
                 var_ptr: {e}"
            ),
                );

                return -1;
            },
        };

        match (
            HANDLE_MANAGER
                .get(n_handle),
            HANDLE_MANAGER
                .get(d_handle),
        ) {
            | (Some(n), Some(d)) => {

                let (q, r) = poly_module::polynomial_long_division(&n, &d, var);

                {

                    *q_handle =
                        HANDLE_MANAGER
                            .insert(q);

                    *r_handle =
                        HANDLE_MANAGER
                            .insert(r);
                }

                0
            },
            | _ => {

                update_last_error("Invalid handle passed to rssn_poly_long_division".to_string());

                -1
            },
        }
    }
}

#[derive(Deserialize)]

struct StatsDataInput {
    data: Vec<f64>,
}

#[derive(Deserialize)]

struct StatsTwoDataInput {
    data1: Vec<f64>,
    data2: Vec<f64>,
}

#[derive(Deserialize)]

struct PercentileInput {
    data: Vec<f64>,
    p: f64,
}

#[derive(Deserialize)]

struct RegressionInput {
    data: Vec<(f64, f64)>,
}

#[derive(Serialize)]

struct RegressionResult {
    slope: f64,
    intercept: f64,
}

macro_rules! impl_stats_fn_single_data {
    ($fn_name:ident, $wrapped_fn:ident) => {
        /// Computes `$wrapped_fn` for a single vector of `f64` data passed as a JSON string.
        #[unsafe(no_mangle)]
        #[deprecated(
            since = "0.1.19",
            note = "This module is deprecated \
                    and please use special \
                    ffi api modules instead."
        )]
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        pub unsafe extern "C" fn $fn_name(json_ptr : *const c_char) -> *mut c_char {

            if json_ptr.is_null() {

                let result = FfiResult {
                    ok : None::<f64>,
                    err : Some("Null pointer passed to function".to_string()),
                };

                return match serde_json::to_string(&result) {
                    | Ok(json_str) => {
                        match CString::new(json_str) {
                            | Ok(c_str) => c_str.into_raw(),
                            | Err(_) => ptr::null_mut(),
                        }
                    },
                    | Err(_) => ptr::null_mut(),
                };
            }

            let json_str = match unsafe {

                CStr::from_ptr(json_ptr).to_str()
            } {
                | Ok(s) => s,
                | Err(_) => return ptr::null_mut(),
            };

            let input : Result<StatsDataInput, _> = serde_json::from_str(json_str);

            let ffi_result = match input {
                | Ok(mut input_data) => {

                    let result = stats_module::$wrapped_fn(&mut input_data.data);

                    FfiResult {
                        ok : Some(result),
                        err : None::<String>,
                    }
                },
                | Err(e) => {
                    FfiResult {
                        ok : None,
                        err : Some(format!(
                            "JSON deserialization error: {}",
                            e
                        )),
                    }
                },
            };

            match serde_json::to_string(&ffi_result) {
                | Ok(json_str) => {
                    match CString::new(json_str) {
                        | Ok(c_str) => c_str.into_raw(),
                        | Err(_) => ptr::null_mut(),
                    }
                },
                | Err(_) => ptr::null_mut(),
            }
        }
    };
}

impl_stats_fn_single_data!(
    stats_mean,
    mean
);

impl_stats_fn_single_data!(
    stats_variance,
    variance
);

impl_stats_fn_single_data!(
    stats_std_dev,
    std_dev
);

impl_stats_fn_single_data!(
    stats_median,
    median
);

impl_stats_fn_single_data!(
    stats_min,
    min
);

impl_stats_fn_single_data!(
    stats_max,
    max
);

impl_stats_fn_single_data!(
    stats_skewness,
    skewness
);

impl_stats_fn_single_data!(
    stats_kurtosis,
    kurtosis
);

impl_stats_fn_single_data!(
    stats_shannon_entropy,
    shannon_entropy
);

/// Computes the percentile of a dataset.
///
/// This function is deprecated.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_stats_percentile \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn stats_percentile(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        let result = FfiResult {
            ok: None::<f64>,
            err: Some(
                "Null pointer passed \
                 to stats_percentile"
                    .to_string(),
            ),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let input: Result<
        PercentileInput,
        _,
    > = serde_json::from_str(json_str);

    let ffi_result = match input {
        | Ok(mut input_data) => {

            let result = stats_module::percentile(
                &mut input_data.data,
                input_data.p,
            );

            FfiResult {
                ok: Some(result),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(format!(
                "JSON deserialization \
                 error: {e}"
            )),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

macro_rules! impl_stats_fn_two_data {
    ($fn_name:ident, $wrapped_fn:ident) => {
        /// Computes `$wrapped_fn` for two vectors of `f64` data passed as a JSON string.
        #[unsafe(no_mangle)]
        #[deprecated(
            since = "0.1.19",
            note = "This module is deprecated \
                    and please use special \
                    ffi api modules instead."
        )]
        /// # Safety
        ///
        /// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
        /// The caller must ensure:
        /// 1. All pointer arguments are valid and point to initialized memory.
        /// 2. The memory layout of passed structures matches the expected C-ABI layout.
        /// 3. Any pointers returned by this function are managed according to the API's ownership rules.
        pub unsafe extern "C" fn $fn_name(json_ptr : *const c_char) -> *mut c_char {


            if json_ptr.is_null() {

                let result = FfiResult {
                    ok : None::<f64>,
                    err : Some("Null pointer passed to function".to_string()),
                };

                return match serde_json::to_string(&result) {
                    | Ok(json_str) => {
                        match CString::new(json_str) {
                            | Ok(c_str) => c_str.into_raw(),
                            | Err(_) => ptr::null_mut(),
                        }
                    },
                    | Err(_) => ptr::null_mut(),
                };
            }

            let json_str = match unsafe {

                CStr::from_ptr(json_ptr).to_str()
            } {
                | Ok(s) => s,
                | Err(_) => return ptr::null_mut(),
            };

            let input : Result<StatsTwoDataInput, _> = serde_json::from_str(json_str);

            let ffi_result = match input {
                | Ok(input_data) => {

                    let result = stats_module::$wrapped_fn(
                        &input_data.data1,
                        &input_data.data2,
                    );

                    FfiResult {
                        ok : Some(result),
                        err : None::<String>,
                    }
                },
                | Err(e) => {
                    FfiResult {
                        ok : None,
                        err : Some(format!(
                            "JSON deserialization error: {}",
                            e
                        )),
                    }
                },
            };

            match serde_json::to_string(&ffi_result) {
                | Ok(json_str) => {
                    match CString::new(json_str) {
                        | Ok(c_str) => c_str.into_raw(),
                        | Err(_) => ptr::null_mut(),
                    }
                },
                | Err(_) => ptr::null_mut(),
            }
        }
    };
}

impl_stats_fn_two_data!(
    stats_covariance,
    covariance
);

impl_stats_fn_two_data!(
    stats_correlation,
    correlation
);

/// Performs a simple linear regression on a dataset.
///
/// This function is deprecated.
#[deprecated(
    since = "0.1.6",
    note = "Please use rssn_stats_simple_linear_regression instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn stats_simple_linear_regression(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        let result = FfiResult {
            ok : None::<RegressionResult>,
            err : Some("Null pointer passed to stats_simple_linear_regression".to_string()),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let input: Result<
        RegressionInput,
        _,
    > = serde_json::from_str(json_str);

    let ffi_result = match input {
        | Ok(input_data) => {

            let (slope, intercept) = simple_linear_regression(&input_data.data);

            let result =
                RegressionResult {
                    slope,
                    intercept,
                };

            FfiResult {
                ok: Some(result),
                err: None::<String>,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(format!(
                "JSON deserialization \
                 error: {e}"
            )),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

/// Computes the mean of a slice of f64 values.
///
/// # Arguments
/// * `data` - A pointer to the first element of the data slice.
/// * `len` - The length of the data slice.
/// * `result` - A pointer to store the computed mean.
///
/// # Returns
/// 0 on success, -1 on error (e.g., null pointer).
///
/// # Safety
/// The `data` pointer must point to a valid array of `len` f64 elements.
/// The `result` pointer must point to a valid f64 location.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_stats_mean(
    data: *const f64,
    len: usize,
    result: *mut f64,
) -> i32 {

    unsafe {

        if data.is_null()
            || result.is_null()
        {

            update_last_error(
                "Null pointer passed \
                 to rssn_stats_mean"
                    .to_string(),
            );

            return -1;
        }

        let data_slice = {

            std::slice::from_raw_parts(
                data, len,
            )
        };

        {

            *result =
                stats_module::mean(
                    data_slice,
                );
        };

        0
    }
}

/// Computes the variance of a slice of f64 values.
///
/// # Arguments
/// * `data` - A pointer to the first element of the data slice.
/// * `len` - The length of the data slice.
/// * `result` - A pointer to store the computed variance.
///
/// # Returns
/// 0 on success, -1 on error (e.g., null pointer).
///
/// # Safety
/// The `data` pointer must point to a valid array of `len` f64 elements.
/// The `result` pointer must point to a valid f64 location.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_stats_variance(
    data: *const f64,
    len: usize,
    result: *mut f64,
) -> i32 {

    unsafe {

        if data.is_null()
            || result.is_null()
        {

            update_last_error(
            "Null pointer passed to \
             rssn_stats_variance"
                .to_string(),
        );

            return -1;
        }

        let data_slice = {

            std::slice::from_raw_parts(
                data, len,
            )
        };

        {

            *result =
                stats_module::variance(
                    data_slice,
                );
        };

        0
    }
}

/// Computes the standard deviation of a slice of f64 values.
///
/// # Arguments
/// * `data` - A pointer to the first element of the data slice.
/// * `len` - The length of the data slice.
/// * `result` - A pointer to store the computed standard deviation.
///
/// # Returns
/// 0 on success, -1 on error (e.g., null pointer).
///
/// # Safety
/// The `data` pointer must point to a valid array of `len` f64 elements.
/// The `result` pointer must point to a valid f64 location.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_stats_std_dev(
    data: *const f64,
    len: usize,
    result: *mut f64,
) -> i32 {

    unsafe {

        if data.is_null()
            || result.is_null()
        {

            update_last_error(
                "Null pointer passed \
                 to rssn_stats_std_dev"
                    .to_string(),
            );

            return -1;
        }

        let data_slice = {

            std::slice::from_raw_parts(
                data, len,
            )
        };

        {

            *result =
                stats_module::std_dev(
                    data_slice,
                );
        };

        0
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Computes the covariance between two datasets.
///
/// # Safety
/// * `d1` and `d2` must be valid pointers to arrays of length `len`.
/// * `result` must be a valid pointer to store the f64 result.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_stats_covariance(
    d1: *const f64,
    d2: *const f64,
    len: usize,
    result: *mut f64,
) -> i32 {

    unsafe {

        if d1.is_null()
            || d2.is_null()
            || result.is_null()
        {

            update_last_error(
            "Null pointer passed to \
             rssn_stats_covariance"
                .to_string(),
        );

            return -1;
        }

        let s1 = {

            std::slice::from_raw_parts(
                d1, len,
            )
        };

        let s2 = {

            std::slice::from_raw_parts(
                d2, len,
            )
        };

        {

            *result =
            stats_module::covariance(
                s1, s2,
            );
        };

        0
    }
}

use crate::numerical::calculus::gradient;
use crate::numerical::integrate::QuadratureMethod;
use crate::numerical::integrate::quadrature;
use crate::numerical::interpolate::{
    self as interp_module,
};
use crate::numerical::stats::simple_linear_regression;
use crate::numerical::stats::{
    self as stats_module,
};
use crate::physics::physics_sm::solve_advection_diffusion_1d;
use crate::symbolic::calculus::definite_integrate;
use crate::symbolic::calculus::differentiate;
use crate::symbolic::calculus::integrate;
use crate::symbolic::calculus::limit;
use crate::symbolic::calculus::substitute;
use crate::symbolic::matrix::add_matrices;
use crate::symbolic::matrix::characteristic_polynomial;
use crate::symbolic::matrix::determinant;
use crate::symbolic::matrix::eigen_decomposition;
use crate::symbolic::matrix::identity_matrix;
use crate::symbolic::matrix::inverse_matrix;
use crate::symbolic::matrix::lu_decomposition;
use crate::symbolic::matrix::mul_matrices;
use crate::symbolic::matrix::null_space;
use crate::symbolic::matrix::rref;
use crate::symbolic::matrix::scalar_mul_matrix;
use crate::symbolic::matrix::sub_matrices;
use crate::symbolic::matrix::trace;
use crate::symbolic::matrix::transpose_matrix;
use crate::symbolic::polynomial::from_coeffs_to_expr;
use crate::symbolic::polynomial::{
    self as poly_module,
};
use crate::symbolic::solve::solve;

/// Differentiates an `Expr` and returns a handle to the new, derivative expression.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_expr_differentiate \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn expr_differentiate(
    handle: *mut Expr,
    var_ptr: *const c_char,
) -> *mut Expr {

    if handle.is_null()
        || var_ptr.is_null()
    {

        return ptr::null_mut();
    }

    let expr = unsafe {

        &*handle
    };

    let var = match unsafe {

        CStr::from_ptr(var_ptr).to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let derivative_expr =
        differentiate(expr, var);

    Arc::into_raw(Arc::new(
        derivative_expr,
    ))
    .cast_mut()
}

/// Substitutes a variable in an `Expr` with another `Expr` and returns a handle to the new expression.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_expr_substitute \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn expr_substitute(
    handle: *mut Expr,
    var_ptr: *const c_char,
    replacement_handle: *mut Expr,
) -> *mut Expr {

    if handle.is_null()
        || var_ptr.is_null()
        || replacement_handle.is_null()
    {

        return ptr::null_mut();
    }

    let expr = unsafe {

        &*handle
    };

    let var = match unsafe {

        CStr::from_ptr(var_ptr).to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let replacement = unsafe {

        &*replacement_handle
    };

    let substituted_expr = substitute(
        expr,
        var,
        replacement,
    );

    Arc::into_raw(Arc::new(
        substituted_expr,
    ))
    .cast_mut()
}

/// Computes the indefinite integral of an `Expr` and returns a handle to the new expression.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_expr_integrate \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn expr_integrate(
    handle: *mut Expr,
    var_ptr: *const c_char,
) -> *mut Expr {

    if handle.is_null()
        || var_ptr.is_null()
    {

        return ptr::null_mut();
    }

    let expr = unsafe {

        &*handle
    };

    let var = match unsafe {

        CStr::from_ptr(var_ptr).to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let integral_expr = integrate(
        expr, var, None, None,
    );

    Arc::into_raw(Arc::new(
        integral_expr,
    ))
    .cast_mut()
}

/// Computes the definite integral of an `Expr` and returns a handle to the new expression.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_definite_integrate \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn expr_definite_integrate(
    handle: *mut Expr,
    var_ptr: *const c_char,
    lower_handle: *mut Expr,
    upper_handle: *mut Expr,
) -> *mut Expr {

    if handle.is_null()
        || var_ptr.is_null()
        || lower_handle.is_null()
        || upper_handle.is_null()
    {

        return ptr::null_mut();
    }

    let expr = unsafe {

        &*handle
    };

    let var = match unsafe {

        CStr::from_ptr(var_ptr).to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let lower = unsafe {

        &*lower_handle
    };

    let upper = unsafe {

        &*upper_handle
    };

    let integral_expr =
        definite_integrate(
            expr, var, lower, upper,
        );

    Arc::into_raw(Arc::new(
        integral_expr,
    ))
    .cast_mut()
}

/// Computes the limit of an `Expr` and returns a handle to the new expression.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_expr_limit instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn expr_limit(
    handle: *mut Expr,
    var_ptr: *const c_char,
    to_handle: *mut Expr,
) -> *mut Expr {

    if handle.is_null()
        || var_ptr.is_null()
        || to_handle.is_null()
    {

        return ptr::null_mut();
    }

    let expr = unsafe {

        &*handle
    };

    let var = match unsafe {

        CStr::from_ptr(var_ptr).to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let to = unsafe {

        &*to_handle
    };

    let limit_expr =
        limit(expr, var, to);

    Arc::into_raw(Arc::new(limit_expr))
        .cast_mut()
}

/// Solves an equation for a given variable and returns the solutions as a JSON string.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_expr_solve instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn expr_solve(
    handle: *mut Expr,
    var_ptr: *const c_char,
) -> *mut c_char {

    if handle.is_null()
        || var_ptr.is_null()
    {

        let result = FfiResult {
            ok: None::<Vec<Expr>>,
            err: Some(
                "Null pointer passed \
                 to expr_solve"
                    .to_string(),
            ),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let expr = unsafe {

        &*handle
    };

    let var = match unsafe {

        CStr::from_ptr(var_ptr).to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let solutions = solve(expr, var);

    let ffi_result = FfiResult {
        ok: Some(solutions),
        err: None::<String>,
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

/// Adds two matrices and returns a handle to the new matrix expression.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_matrix_add instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn matrix_add(
    h1: *mut Expr,
    h2: *mut Expr,
) -> *mut Expr {

    if h1.is_null() || h2.is_null() {

        return ptr::null_mut();
    }

    let m1 = unsafe {

        &*h1
    };

    let m2 = unsafe {

        &*h2
    };

    Arc::into_raw(Arc::new(
        add_matrices(m1, m2),
    ))
    .cast_mut()
}

/// Subtracts the second matrix from the first and returns a handle to the new matrix expression.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_matrix_sub instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn matrix_sub(
    h1: *mut Expr,
    h2: *mut Expr,
) -> *mut Expr {

    if h1.is_null() || h2.is_null() {

        return ptr::null_mut();
    }

    let m1 = unsafe {

        &*h1
    };

    let m2 = unsafe {

        &*h2
    };

    Arc::into_raw(Arc::new(
        sub_matrices(m1, m2),
    ))
    .cast_mut()
}

/// Multiplies two matrices and returns a handle to the new matrix expression.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_matrix_mul instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn matrix_mul(
    h1: *mut Expr,
    h2: *mut Expr,
) -> *mut Expr {

    if h1.is_null() || h2.is_null() {

        return ptr::null_mut();
    }

    let m1 = unsafe {

        &*h1
    };

    let m2 = unsafe {

        &*h2
    };

    Arc::into_raw(Arc::new(
        mul_matrices(m1, m2),
    ))
    .cast_mut()
}

/// Transposes a matrix and returns a handle to the new matrix expression.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_matrix_transpose \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn matrix_transpose(
    handle: *mut Expr
) -> *mut Expr {

    if handle.is_null() {

        return ptr::null_mut();
    }

    let m = unsafe {

        &*handle
    };

    Arc::into_raw(Arc::new(
        transpose_matrix(m),
    ))
    .cast_mut()
}

/// Computes the determinant of a matrix and returns a handle to the resulting expression.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_matrix_determinant \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn matrix_determinant(
    handle: *mut Expr
) -> *mut Expr {

    if handle.is_null() {

        return ptr::null_mut();
    }

    let m = unsafe {

        &*handle
    };

    Arc::into_raw(Arc::new(
        determinant(m),
    ))
    .cast_mut()
}

/// Inverts a matrix and returns a handle to the new matrix expression.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_matrix_inverse \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn matrix_inverse(
    handle: *mut Expr
) -> *mut Expr {

    if handle.is_null() {

        return ptr::null_mut();
    }

    let m = unsafe {

        &*handle
    };

    Arc::into_raw(Arc::new(
        inverse_matrix(m),
    ))
    .cast_mut()
}

/// Creates an identity matrix of a given size and returns a handle to it.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_matrix_identity \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn matrix_identity(
    size: usize
) -> *mut Expr {

    Arc::into_raw(Arc::new(
        identity_matrix(size),
    ))
    .cast_mut()
}

/// Multiplies a matrix by a scalar and returns a handle to the new matrix expression.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_matrix_scalar_mul \
            instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn matrix_scalar_mul(
    scalar_handle: *mut Expr,
    matrix_handle: *mut Expr,
) -> *mut Expr {

    if scalar_handle.is_null()
        || matrix_handle.is_null()
    {

        return ptr::null_mut();
    }

    let scalar = unsafe {

        &*scalar_handle
    };

    let matrix = unsafe {

        &*matrix_handle
    };

    Arc::into_raw(Arc::new(
        scalar_mul_matrix(
            scalar,
            matrix,
        ),
    ))
    .cast_mut()
}

/// Computes the trace of a matrix and returns the result as a JSON string.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_matrix_trace instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn matrix_trace(
    handle: *mut Expr
) -> *mut c_char {

    if handle.is_null() {

        let result = FfiResult {
            ok: None::<Expr>,
            err: Some(
                "Null pointer passed \
                 to matrix_trace"
                    .to_string(),
            ),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let m = unsafe {

        &*handle
    };

    let result = trace(m);

    let ffi_result = match result {
        | Ok(value) => {
            FfiResult {
                ok: Some(value),
                err: None,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

/// Computes the characteristic polynomial of a matrix and returns the result as a JSON string.
#[deprecated(
    since = "0.1.6",
    note = "Please use rssn_matrix_characteristic_polynomial instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn matrix_characteristic_polynomial(
    handle: *mut Expr,
    var_ptr: *const c_char,
) -> *mut c_char {

    if handle.is_null()
        || var_ptr.is_null()
    {

        let result = FfiResult {
            ok : None::<Expr>,
            err : Some("Null pointer passed to matrix_characteristic_polynomial".to_string()),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let m = unsafe {

        &*handle
    };

    let var = match unsafe {

        CStr::from_ptr(var_ptr).to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let result =
        characteristic_polynomial(
            m, var,
        );

    let ffi_result = match result {
        | Ok(value) => {
            FfiResult {
                ok: Some(value),
                err: None,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

/// Computes the Reduced Row Echelon Form (RREF) of a matrix and returns the result as a JSON string.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_matrix_rref instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn matrix_rref(
    handle: *mut Expr
) -> *mut c_char {

    if handle.is_null() {

        let result = FfiResult {
            ok: None::<Expr>,
            err: Some(
                "Null pointer passed \
                 to matrix_rref"
                    .to_string(),
            ),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let m = unsafe {

        &*handle
    };

    let result = rref(m);

    let ffi_result = match result {
        | Ok(value) => {
            FfiResult {
                ok: Some(value),
                err: None,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

/// Computes the null space of a matrix and returns the result as a JSON string.
#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_null_space instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn matrix_null_space(
    handle: *mut Expr
) -> *mut c_char {

    if handle.is_null() {

        let result = FfiResult {
            ok: None::<Expr>,
            err: Some(
                "Null pointer passed \
                 to matrix_null_space"
                    .to_string(),
            ),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let m = unsafe {

        &*handle
    };

    let result = null_space(m);

    let ffi_result = match result {
        | Ok(value) => {
            FfiResult {
                ok: Some(value),
                err: None,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

#[derive(Serialize)]

struct MatrixPair {
    p1: Expr,
    p2: Expr,
}

/// Computes the LU decomposition of a matrix and returns the L and U matrices as a JSON string.
#[deprecated(
    since = "0.1.6",
    note = "Please use rssn_matrix_lu_decomposition instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn matrix_lu_decomposition(
    handle: *mut Expr
) -> *mut c_char {

    if handle.is_null() {

        let result = FfiResult {
            ok : None::<MatrixPair>,
            err : Some("Null pointer passed to matrix_lu_decomposition".to_string()),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let m = unsafe {

        &*handle
    };

    let result = lu_decomposition(m);

    let ffi_result = match result {
        | Ok((l, u)) => {
            FfiResult {
                ok: Some(MatrixPair {
                    p1: l,
                    p2: u,
                }),
                err: None,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

/// Computes the eigenvalue decomposition of a matrix and returns the eigenvalues and eigenvectors as a JSON string.
#[deprecated(
    since = "0.1.6",
    note = "Please use rssn_matrix_eigen_decomposition instead."
)]
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn matrix_eigen_decomposition(
    handle: *mut Expr
) -> *mut c_char {

    if handle.is_null() {

        let result = FfiResult {
            ok : None::<MatrixPair>,
            err : Some("Null pointer passed to matrix_eigen_decomposition".to_string()),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let m = unsafe {

        &*handle
    };

    let result = eigen_decomposition(m);

    let ffi_result =
        match result {
            | Ok((
                eigenvalues,
                eigenvectors,
            )) => FfiResult {
                ok: Some(MatrixPair {
                    p1: eigenvalues,
                    p2: eigenvectors,
                }),
                err: None,
            },
            | Err(e) => {
                FfiResult {
                    ok: None,
                    err: Some(e),
                }
            },
        };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

#[derive(Serialize, Deserialize)]

struct GradientInput {
    expr: Expr,
    vars: Vec<String>,
    point: Vec<f64>,
}

#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_numerical_gradient \
            instead."
)]
#[unsafe(no_mangle)]

/// Computes the numerical gradient of an expression.
///
/// This function is deprecated.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn numerical_gradient(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        let result = FfiResult {
            ok: None::<Vec<f64>>,
            err: Some(
                "Null pointer passed \
                 to numerical_gradient"
                    .to_string(),
            ),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let input: Result<
        GradientInput,
        _,
    > = serde_json::from_str(json_str);

    let ffi_result = match input {
        | Ok(grad_input) => {

            let vars_as_str: Vec<&str> =
                grad_input
                    .vars
                    .iter()
                    .map(std::string::String::as_str)
                    .collect();

            let grad_result = gradient(
                &grad_input.expr,
                &vars_as_str,
                &grad_input.point,
            );

            match grad_result {
                | Ok(grad_vec) => {
                    FfiResult {
                        ok: Some(
                            grad_vec,
                        ),
                        err: None,
                    }
                },
                | Err(e) => {
                    FfiResult {
                        ok: None,
                        err: Some(e),
                    }
                },
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(format!(
                "JSON deserialization \
                 error: {e}"
            )),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

#[derive(Deserialize)]

enum FfiQuadratureMethod {
    Trapezoidal,
    Simpson,
}

#[derive(Deserialize)]

struct IntegrationInput {
    expr: Expr,
    var: String,
    start: f64,
    end: f64,
    n_steps: usize,
    method: FfiQuadratureMethod,
}

#[deprecated(
    since = "0.1.6",
    note = "Please use \
            rssn_numerical_integrate \
            instead."
)]
#[unsafe(no_mangle)]

/// Performs numerical integration of an expression.
///
/// This function is deprecated.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn numerical_integrate(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        let result = FfiResult {
            ok : None::<f64>,
            err : Some("Null pointer passed to numerical_integrate".to_string()),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let input: Result<
        IntegrationInput,
        _,
    > = serde_json::from_str(json_str);

    let ffi_result = match input {
        | Ok(int_input) => {

            let method = match int_input.method {
                | FfiQuadratureMethod::Trapezoidal => QuadratureMethod::Trapezoidal,
                | FfiQuadratureMethod::Simpson => QuadratureMethod::Simpson,
            };

            let int_result = quadrature(
                &int_input.expr,
                &int_input.var,
                (
                    int_input.start,
                    int_input.end,
                ),
                int_input.n_steps,
                &method,
            );

            match int_result {
                | Ok(val) => {
                    FfiResult {
                        ok: Some(val),
                        err: None,
                    }
                },
                | Err(e) => {
                    FfiResult {
                        ok: None,
                        err: Some(e),
                    }
                },
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(format!(
                "JSON deserialization \
                 error: {e}"
            )),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

#[derive(Deserialize)]

struct AdvectionDiffusion1DInput {
    initial_condition: Vec<f64>,
    dx: f64,
    c: f64,
    d: f64,
    dt: f64,
    steps: usize,
}

#[deprecated(
    since = "0.1.6",
    note = "Please use rssn_physics_solve_advection_diffusion_1d instead."
)]
#[unsafe(no_mangle)]

/// Solves the 1D advection-diffusion equation.
///
/// This function is deprecated.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn physics_solve_advection_diffusion_1d(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        let result = FfiResult {
            ok : None::<Vec<f64>>,
            err : Some("Null pointer passed to physics_solve_advection_diffusion_1d".to_string()),
        };

        return match serde_json::to_string(&result) {
            | Ok(json_str) => {
                match CString::new(json_str) {
                    | Ok(c_str) => c_str.into_raw(),
                    | Err(_) => ptr::null_mut(),
                }
            },
            | Err(_) => ptr::null_mut(),
        };
    }

    let json_str = match unsafe {

        CStr::from_ptr(json_ptr)
            .to_str()
    } {
        | Ok(s) => s,
        | Err(_) => {
            return ptr::null_mut();
        },
    };

    let input: Result<
        AdvectionDiffusion1DInput,
        _,
    > = serde_json::from_str(json_str);

    let ffi_result = match input {
        | Ok(sim_input) => {

            let result_vec = solve_advection_diffusion_1d(
                &sim_input.initial_condition,
                sim_input.dx,
                sim_input.c,
                sim_input.d,
                sim_input.dt,
                sim_input.steps,
            );

            FfiResult {
                ok: Some(result_vec),
                err: None,
            }
        },
        | Err(e) => {
            FfiResult {
                ok: None,
                err: Some(format!(
                "JSON deserialization \
                 error: {e}"
            )),
            }
        },
    };

    match serde_json::to_string(
        &ffi_result,
    ) {
        | Ok(json_str) => {
            match CString::new(json_str)
            {
                | Ok(c_str) => {
                    c_str.into_raw()
                },
                | Err(_) => {
                    ptr::null_mut()
                },
            }
        },
        | Err(_) => ptr::null_mut(),
    }
}

/// A point in 2D space for FFI communication.
#[repr(C)]

pub struct FfiPoint {
    x: f64,
    y: f64,
}

/// Computes a Lagrange interpolating polynomial from a set of points.
/// Returns a handle to the resulting polynomial expression.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_interp_lagrange(
    points_ptr: *const FfiPoint,
    num_points: usize,
    result_handle: *mut usize,
) -> i32 {

    unsafe {

        if points_ptr.is_null()
            || result_handle.is_null()
        {

            update_last_error(
            "Null pointer passed to \
             rssn_interp_lagrange"
                .to_string(),
        );

            return -1;
        }

        let points_slice = {

            std::slice::from_raw_parts(
                points_ptr,
                num_points,
            )
        };

        let points_vec: Vec<(
            f64,
            f64,
        )> = points_slice
            .iter()
            .map(|p| (p.x, p.y))
            .collect();

        match interp_module::lagrange_interpolation(&points_vec) {
        | Ok(poly) => {

            let expr_coeffs = poly
                .coeffs
                .into_iter()
                .map(Expr::Constant)
                .collect();

            let poly_expr = Expr::Polynomial(expr_coeffs);

             {

                *result_handle = HANDLE_MANAGER.insert(poly_expr);
            };

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
    }
}

/// Evaluates a point on a Bezier curve defined by control points.
#[unsafe(no_mangle)]
#[allow(clippy::indexing_slicing)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_interp_bezier_curve(
    points_ptr: *const FfiPoint,
    num_points: usize,
    t: f64,
    result_ptr: *mut FfiPoint,
) -> i32 {

    unsafe {

        if points_ptr.is_null()
            || result_ptr.is_null()
        {

            update_last_error(
            "Null pointer passed to \
             rssn_interp_bezier_curve"
                .to_string(),
        );

            return -1;
        }

        let points_slice = {

            std::slice::from_raw_parts(
                points_ptr,
                num_points,
            )
        };

        let control_points: Vec<
            Vec<f64>,
        > = points_slice
            .iter()
            .map(|p| vec![p.x, p.y])
            .collect();

        let result_vec =
            interp_module::bezier_curve(
                &control_points,
                t,
            );

        if result_vec.len() >= 2 {

            {

                (*result_ptr).x =
                    result_vec[0];

                (*result_ptr).y =
                    result_vec[1];
            }

            0
        } else {

            update_last_error(
                "Bezier curve \
                 evaluation returned \
                 an invalid point"
                    .to_string(),
            );

            -1
        }
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Computes the numerical definite integral of an expression using handles.
///
/// # Safety
/// * `var` must be a valid null-terminated C string.
/// * `result` must be a valid pointer to store the f64 result.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_numerical_integrate(
    expr_h: usize,
    var: *const c_char,
    start: f64,
    end: f64,
    n_steps: usize,
    method: u32,
    result: *mut f64,
) -> i32 {

    unsafe {

        if var.is_null()
            || result.is_null()
        {

            update_last_error(
            "Null pointer passed to \
             rssn_numerical_integrate"
                .to_string(),
        );

            return -1;
        }

        let var_str =
            match CStr::from_ptr(var)
                .to_str()
            {
                | Ok(s) => s,
                | Err(e) => {

                    update_last_error(
                        format!(
                "Invalid UTF-8 in \
                 var: {e}"
            ),
                    );

                    return -1;
                },
            };

        let quad_method = match method {
        | 0 => QuadratureMethod::Trapezoidal,
        | 1 => QuadratureMethod::Simpson,
        | _ => {

            update_last_error("Invalid quadrature method specified".to_string());

            return -1;
        },
    };

        match HANDLE_MANAGER.get(expr_h)
        {
            | Some(expr) => {
                match quadrature(
                    &expr,
                    var_str,
                    (start, end),
                    n_steps,
                    &quad_method,
                ) {
                    | Ok(val) => {

                        {

                            *result =
                                val;
                        };

                        0
                    },
                    | Err(e) => {

                        update_last_error(
                        e,
                    );

                        -1
                    },
                }
            },
            | None => {

                update_last_error(format!(
                "Invalid handle passed to rssn_numerical_integrate: {expr_h}"
            ));

                -1
            },
        }
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Subtracts two matrices using handles.
///
/// # Safety
/// * `result_h` must be a valid pointer to store the resulting handle.

pub unsafe extern "C" fn rssn_matrix_sub(
    h1: usize,
    h2: usize,
    result_h: *mut usize,
) -> i32 {

    unsafe {

        if result_h.is_null() {

            update_last_error(
                "Null pointer passed \
                 to rssn_matrix_sub"
                    .to_string(),
            );

            return -1;
        }

        match (
            HANDLE_MANAGER.get(h1),
            HANDLE_MANAGER.get(h2),
        ) {
            | (Some(m1), Some(m2)) => {

                let result =
                    sub_matrices(
                        &m1, &m2,
                    );

                {

                    *result_h =
                        HANDLE_MANAGER
                            .insert(
                                result,
                            );
                };

                0
            },
            | _ => {

                update_last_error("Invalid handle passed to rssn_matrix_sub".to_string());

                -1
            },
        }
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Multiplies two matrices using handles.
///
/// # Safety
/// * `result_h` must be a valid pointer to store the resulting handle.

pub unsafe extern "C" fn rssn_matrix_mul(
    h1: usize,
    h2: usize,
    result_h: *mut usize,
) -> i32 {

    unsafe {

        if result_h.is_null() {

            update_last_error(
                "Null pointer passed \
                 to rssn_matrix_mul"
                    .to_string(),
            );

            return -1;
        }

        match (
            HANDLE_MANAGER.get(h1),
            HANDLE_MANAGER.get(h2),
        ) {
            | (Some(m1), Some(m2)) => {

                let result =
                    mul_matrices(
                        &m1, &m2,
                    );

                {

                    *result_h =
                        HANDLE_MANAGER
                            .insert(
                                result,
                            );
                };

                0
            },
            | _ => {

                update_last_error("Invalid handle passed to rssn_matrix_mul".to_string());

                -1
            },
        }
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Transposes a matrix using a handle.
///
/// # Safety
/// * `result_h` must be a valid pointer to store the resulting handle.

pub unsafe extern "C" fn rssn_matrix_transpose(
    h: usize,
    result_h: *mut usize,
) -> i32 {

    unsafe {

        if result_h.is_null() {

            update_last_error(
            "Null pointer passed to \
             rssn_matrix_transpose"
                .to_string(),
        );

            return -1;
        }

        match HANDLE_MANAGER.get(h) {
            | Some(m) => {

                let result =
                    transpose_matrix(
                        &m,
                    );

                {

                    *result_h =
                        HANDLE_MANAGER
                            .insert(
                                result,
                            );
                };

                0
            },
            | _ => {

                update_last_error("Invalid handle passed to rssn_matrix_transpose".to_string());

                -1
            },
        }
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Computes the determinant of a matrix using a handle.
///
/// # Safety
/// * `result_h` must be a valid pointer to store the resulting handle to the determinant expression.

pub unsafe extern "C" fn rssn_matrix_determinant(
    h: usize,
    result_h: *mut usize,
) -> i32 {

    unsafe {

        if result_h.is_null() {

            update_last_error(
            "Null pointer passed to \
             rssn_matrix_determinant"
                .to_string(),
        );

            return -1;
        }

        match HANDLE_MANAGER.get(h) {
            | Some(m) => {

                let result =
                    determinant(&m);

                {

                    *result_h =
                        HANDLE_MANAGER
                            .insert(
                                result,
                            );
                };

                0
            },
            | _ => {

                update_last_error("Invalid handle passed to rssn_matrix_determinant".to_string());

                -1
            },
        }
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Computes the inverse of a matrix using a handle.
///
/// # Safety
/// * `result_h` must be a valid pointer to store the resulting handle.

pub unsafe extern "C" fn rssn_matrix_inverse(
    h: usize,
    result_h: *mut usize,
) -> i32 {

    unsafe {

        if result_h.is_null() {

            update_last_error(
            "Null pointer passed to \
             rssn_matrix_inverse"
                .to_string(),
        );

            return -1;
        }

        match HANDLE_MANAGER.get(h) {
            | Some(m) => {

                let result =
                    inverse_matrix(&m);

                {

                    *result_h =
                        HANDLE_MANAGER
                            .insert(
                                result,
                            );
                };

                0
            },
            | _ => {

                update_last_error("Invalid handle passed to rssn_matrix_inverse".to_string());

                -1
            },
        }
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Creates an identity matrix of the specified size.
///
/// # Safety
/// * `result_h` must be a valid pointer to store the resulting handle.

pub unsafe extern "C" fn rssn_matrix_identity(
    size: usize,
    result_h: *mut usize,
) -> i32 {

    unsafe {

        if result_h.is_null() {

            update_last_error(
            "Null pointer passed to \
             rssn_matrix_identity"
                .to_string(),
        );

            return -1;
        }

        let result =
            identity_matrix(size);

        {

            *result_h = HANDLE_MANAGER
                .insert(result);
        };

        0
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Multiplies a matrix by a scalar expression using handles.
///
/// # Safety
/// * `result_h` must be a valid pointer to store the resulting handle.

pub unsafe extern "C" fn rssn_matrix_scalar_mul(
    scalar_h: usize,
    matrix_h: usize,
    result_h: *mut usize,
) -> i32 {

    unsafe {

        if result_h.is_null() {

            update_last_error(
            "Null pointer passed to \
             rssn_matrix_scalar_mul"
                .to_string(),
        );

            return -1;
        }

        match (
            HANDLE_MANAGER
                .get(scalar_h),
            HANDLE_MANAGER
                .get(matrix_h),
        ) {
            | (
                Some(scalar),
                Some(matrix),
            ) => {

                let result =
                    scalar_mul_matrix(
                        &scalar,
                        &matrix,
                    );

                {

                    *result_h =
                        HANDLE_MANAGER
                            .insert(
                                result,
                            );
                };

                0
            },
            | _ => {

                update_last_error("Invalid handle passed to rssn_matrix_scalar_mul".to_string());

                -1
            },
        }
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Differentiates an expression with respect to a variable using handles.
///
/// # Safety
/// * `var` must be a valid null-terminated C string.
/// * `result_h` must be a valid pointer to store the resulting handle.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_calculus_differentiate(
    expr_h: usize,
    var: *const c_char,
    result_h: *mut usize,
) -> i32 {

    unsafe {

        if var.is_null()
            || result_h.is_null()
        {

            update_last_error("Null pointer passed to rssn_calculus_differentiate".to_string());

            return -1;
        }

        let var_str = {

            CStr::from_ptr(var)
                .to_str()
                .expect("var is empty")
        };

        match HANDLE_MANAGER.get(expr_h)
        {
            | Some(expr) => {

                let derivative =
                    differentiate(
                        &expr,
                        var_str,
                    );

                {

                    *result_h =
                        HANDLE_MANAGER
                            .insert(
                            derivative,
                        );
                };

                0
            },
            | None => {

                update_last_error(format!(
                "Invalid handle passed to rssn_calculus_differentiate: {expr_h}"
            ));

                -1
            },
        }
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Substitutes a variable with another expression in an expression using handles.
///
/// # Safety
/// * `var` must be a valid null-terminated C string.
/// * `result_h` must be a valid pointer to store the resulting handle.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_calculus_substitute(
    expr_h: usize,
    var: *const c_char,
    replacement_h: usize,
    result_h: *mut usize,
) -> i32 {

    unsafe {

        if var.is_null()
            || result_h.is_null()
        {

            update_last_error(
            "Null pointer passed to \
             rssn_calculus_substitute"
                .to_string(),
        );

            return -1;
        }

        let var_str = {

            CStr::from_ptr(var)
                .to_str()
                .expect("var is empty")
        };

        match (
            HANDLE_MANAGER.get(expr_h),
            HANDLE_MANAGER
                .get(replacement_h),
        ) {
            | (
                Some(expr),
                Some(rep),
            ) => {

                let result = substitute(
                    &expr,
                    var_str,
                    &rep,
                );

                {

                    *result_h =
                        HANDLE_MANAGER
                            .insert(
                                result,
                            );
                };

                0
            },
            | _ => {

                update_last_error("Invalid handle passed to rssn_calculus_substitute".to_string());

                -1
            },
        }
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Integrates an expression with respect to a variable using handles (indefinite).
///
/// # Safety
/// * `var` must be a valid null-terminated C string.
/// * `result_h` must be a valid pointer to store the resulting handle.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_calculus_integrate(
    expr_h: usize,
    var: *const c_char,
    result_h: *mut usize,
) -> i32 {

    unsafe {

        if var.is_null()
            || result_h.is_null()
        {

            update_last_error(
            "Null pointer passed to \
             rssn_calculus_integrate"
                .to_string(),
        );

            return -1;
        }

        let var_str = {

            CStr::from_ptr(var)
                .to_str()
                .expect("var is empty")
        };

        match HANDLE_MANAGER.get(expr_h)
        {
            | Some(expr) => {

                let integral =
                    integrate(
                        &expr,
                        var_str,
                        None,
                        None,
                    );

                {

                    *result_h =
                        HANDLE_MANAGER
                            .insert(
                            integral,
                        );
                };

                0
            },
            | None => {

                update_last_error(format!(
                "Invalid handle passed to rssn_calculus_integrate: {expr_h}"
            ));

                -1
            },
        }
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Computes the definite integral of an expression using handles.
///
/// # Safety
/// * `var` must be a valid null-terminated C string.
/// * `result_h` must be a valid pointer to store the resulting handle.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_calculus_definite_integrate(
    expr_h: usize,
    var: *const c_char,
    lower_h: usize,
    upper_h: usize,
    result_h: *mut usize,
) -> i32 {

    unsafe {

        if var.is_null()
            || result_h.is_null()
        {

            update_last_error("Null pointer passed to rssn_calculus_definite_integrate".to_string());

            return -1;
        }

        let var_str = {

            CStr::from_ptr(var)
                .to_str()
                .expect("var is empty")
        };

        match (
            HANDLE_MANAGER.get(expr_h),
            HANDLE_MANAGER.get(lower_h),
            HANDLE_MANAGER.get(upper_h),
        ) {
            | (
                Some(expr),
                Some(lower),
                Some(upper),
            ) => {

                let integral =
                    definite_integrate(
                        &expr,
                        var_str,
                        &lower,
                        &upper,
                    );

                {

                    *result_h =
                        HANDLE_MANAGER
                            .insert(
                            integral,
                        );
                };

                0
            },
            | _ => {

                update_last_error(
                "Invalid handle passed to rssn_calculus_definite_integrate".to_string(),
            );

                -1
            },
        }
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Computes the limit of an expression as a variable approaches a value using handles.
///
/// # Safety
/// * `var` must be a valid null-terminated C string.
/// * `result_h` must be a valid pointer to store the resulting handle.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_calculus_limit(
    expr_h: usize,
    var: *const c_char,
    to_h: usize,
    result_h: *mut usize,
) -> i32 {

    unsafe {

        if var.is_null()
            || result_h.is_null()
        {

            update_last_error(
            "Null pointer passed to \
             rssn_calculus_limit"
                .to_string(),
        );

            return -1;
        }

        let var_str = {

            CStr::from_ptr(var)
                .to_str()
                .expect("var is empty")
        };

        match (
            HANDLE_MANAGER.get(expr_h),
            HANDLE_MANAGER.get(to_h),
        ) {
            | (
                Some(expr),
                Some(to),
            ) => {

                let result = limit(
                    &expr,
                    var_str,
                    &to,
                );

                {

                    *result_h =
                        HANDLE_MANAGER
                            .insert(
                                result,
                            );
                };

                0
            },
            | _ => {

                update_last_error("Invalid handle passed to rssn_calculus_limit".to_string());

                -1
            },
        }
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Solves an equation or searches for zeros of an expression.
///
/// # Safety
/// * `var` must be a valid null-terminated C string.
/// * `result_h` must be a valid pointer to store the handle of the solution.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve(
    expr_h: usize,
    var: *const c_char,
    result_h: *mut usize,
) -> i32 {

    unsafe {

        let handle_error =
            |err_msg: String| {

                update_last_error(
                    err_msg,
                );

                -1
            };

        if var.is_null()
            || result_h.is_null()
        {

            return handle_error(
                "Null pointer passed \
                 to rssn_solve"
                    .to_string(),
            );
        }

        let var_str =
            match CStr::from_ptr(var)
                .to_str()
            {
                | Ok(s) => s,
                | Err(e) => {
                    return handle_error(
                        format!(
                    "Invalid UTF-8 in \
                     variable name: \
                     {e}"
                ),
                    );
                },
            };

        let expr = match HANDLE_MANAGER
            .get(expr_h)
        {
            | Some(e) => e,
            | None => {
                return handle_error(
                    format!(
                    "Invalid handle: \
                     {expr_h}"
                ),
                );
            },
        };

        let solutions = crate::numerical::testing::solve(&expr, var_str);

        if solutions.is_empty() {

            {

                *result_h = 0;
            };
        } else {

            // Return first solution for now, or could wrap in BinaryList
            {

                *result_h =
                    HANDLE_MANAGER
                        .insert(
                            solutions
                                [0]
                            .clone(),
                        );
            };
        }

        0
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Adds two matrices using handles.
///
/// # Safety
/// * `result_h` must be a valid pointer to store the resulting handle.

pub unsafe extern "C" fn rssn_matrix_add(
    h1: usize,
    h2: usize,
    result_h: *mut usize,
) -> i32 {

    unsafe {

        let handle_error =
            |err_msg: String| {

                update_last_error(
                    err_msg,
                );

                -1
            };

        if result_h.is_null() {

            return -1;
        }

        let m1 = match HANDLE_MANAGER
            .get(h1)
        {
            | Some(e) => e,
            | None => {
                return handle_error(
                    format!(
                    "Invalid handle \
                     h1: {h1}"
                ),
                );
            },
        };

        let m2 = match HANDLE_MANAGER
            .get(h2)
        {
            | Some(e) => e,
            | None => {
                return handle_error(
                    format!(
                    "Invalid handle \
                     h2: {h2}"
                ),
                );
            },
        };

        let res = crate::symbolic::matrix::add_matrices(&m1, &m2);

        {

            *result_h = HANDLE_MANAGER
                .insert(res);
        };

        0
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Computes the numerical gradient of an expression at a point.
///
/// # Safety
/// * `vars` must be a valid pointer to an array of `num_vars` null-terminated strings.
/// * `point` must be a valid pointer to an array of `point_len` f64 elements.
/// * `result_vec` must point to a buffer of size at least `point_len` f64 elements.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_numerical_gradient(
    expr_h: usize,
    vars: *const *const c_char,
    num_vars: usize,
    point: *const f64,
    point_len: usize,
    result_vec: *mut f64,
) -> i32 {

    unsafe {

        let handle_error =
            |err_msg: String| {

                update_last_error(
                    err_msg,
                );

                -1
            };

        if vars.is_null()
            || point.is_null()
            || result_vec.is_null()
        {

            return handle_error(
            "Null pointer passed to \
             rssn_numerical_gradient"
                .to_string(),
        );
        }

        let expr = match HANDLE_MANAGER
            .get(expr_h)
        {
            | Some(e) => e,
            | None => {
                return handle_error(
                    format!(
                    "Invalid handle: \
                     {expr_h}"
                ),
                );
            },
        };

        let mut vars_vec =
            Vec::with_capacity(
                num_vars,
            );

        for i in 0 .. num_vars {

            let v_ptr = {

                *vars.add(i)
            };

            if v_ptr.is_null() {

                return handle_error(
                    "Null pointer in \
                     vars array"
                        .to_string(),
                );
            }

            let v_str = match CStr::from_ptr(v_ptr).to_str() {
            | Ok(s) => s,
            | Err(e) => {
                return handle_error(format!(
                    "Invalid UTF-8 in var {i}: {e}"
                ))
            },
        };

            vars_vec.push(v_str);
        }

        let point_slice = {

            std::slice::from_raw_parts(
                point,
                point_len,
            )
        };

        match crate::numerical::vector_calculus::gradient(
        &expr,
        &vars_vec,
        point_slice,
    ) {
        | Ok(grad) => {

            if grad.len() != point_len {

                return handle_error(format!(
                    "Gradient length {} != point length {}",
                    grad.len(),
                    point_len
                ));
            }

             {

                std::ptr::copy_nonoverlapping(
                    grad.as_ptr(),
                    result_vec,
                    grad.len(),
                );
            }

            0
        },
        | Err(e) => handle_error(e),
    }
    }
}

#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// Solves the 1D advection-diffusion equation numerically.
///
/// # Safety
/// * `initial_cond` must be a valid pointer to an array of `len` f64 elements.
/// * `result_ptr` must point to a buffer of size at least `len` f64 elements.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_advection_diffusion_1d(
    initial_cond: *const f64,
    len: usize,
    dx: f64,
    c: f64,
    d: f64,
    dt: f64,
    steps: usize,
    result_ptr: *mut f64,
) -> i32 {

    unsafe {

        if initial_cond.is_null()
            || result_ptr.is_null()
        {

            update_last_error("Null pointer passed to rssn_physics_advection_diffusion_1d".to_string());

            return -1;
        }

        let init_slice = {

            std::slice::from_raw_parts(
                initial_cond,
                len,
            )
        };

        let res = crate::physics::physics_fdm::solve_advection_diffusion_1d(
        init_slice,
        dx,
        c,
        d,
        dt,
        steps,
    );

        {

            std::ptr::copy_nonoverlapping(
            res.as_ptr(),
            result_ptr,
            len,
        );
        }

        0
    }
}

/// Computes the greatest common divisor (GCD) of two numbers.
///
/// Returns 0 on success, -1 on error.
/// On error, call `rssn_get_last_error` to get the error message.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_nt_gcd(
    a: u64,
    b: u64,
    result: *mut u64,
) -> i32 {

    unsafe {

        if result.is_null() {

            update_last_error(
                "Null pointer passed \
                 for 'result' to \
                 rssn_nt_gcd"
                    .to_string(),
            );

            return -1;
        }

        {

            *result = nt::gcd(a, b);
        }

        0
    }
}

/// Checks if a number is prime using the Miller-Rabin test.
///
/// Returns 0 on success, -1 on error.
/// On error, call `rssn_get_last_error` to get the error message.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_nt_is_prime(
    n: u64,
    result: *mut bool,
) -> i32 {

    unsafe {

        if result.is_null() {

            update_last_error(
                "Null pointer passed \
                 for 'result' to \
                 rssn_nt_is_prime"
                    .to_string(),
            );

            return -1;
        }

        {

            *result =
            nt::is_prime_miller_rabin(
                n,
            );
        }

        0
    }
}

/// Computes modular exponentiation (base^exp % modulus).
///
/// Returns 0 on success, -1 on error.
/// On error, call `rssn_get_last_error` to get the error message.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_nt_mod_pow(
    base: u64,
    exp: u64,
    modulus: u64,
    result: *mut u64,
) -> i32 {

    unsafe {

        if result.is_null() {

            update_last_error(
                "Null pointer passed \
                 for 'result' to \
                 rssn_nt_mod_pow"
                    .to_string(),
            );

            return -1;
        }

        {

            *result = nt::mod_pow(
                u128::from(base),
                exp,
                modulus,
            );
        }

        0
    }
}

/// Computes the modular multiplicative inverse.
///
/// Returns 0 on success, -1 on error (e.g., if no inverse exists).
/// On error, call `rssn_get_last_error` to get the error message.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_nt_mod_inverse(
    a: i64,
    b: i64,
    result: *mut i64,
) -> i32 {

    unsafe {

        if result.is_null() {

            update_last_error(
                "Null pointer passed \
                 for 'result' to \
                 rssn_nt_mod_inverse"
                    .to_string(),
            );

            return -1;
        }

        match nt::mod_inverse(a, b) {
            | Some(inverse) => {

                {

                    *result = inverse;
                };

                0
            },
            | None => {

                update_last_error(
                    format!(
                "Modular inverse of \
                 {a} modulo {b} does \
                 not exist."
            ),
                );

                -1
            },
        }
    }
}

static PLUGIN_MANAGER:
    std::sync::LazyLock<
        Mutex<Option<PluginManager>>,
    > =
    std::sync::LazyLock::new(|| {

        Mutex::new(None)
    });

/// Initializes the plugin manager with a specified plugin directory.
///
/// This function must be called before any plugin operations are performed.
///
/// # Arguments
/// * `plugin_dir_ptr` - A null-terminated UTF-8 string for the plugin directory path.
///
/// # Returns
/// 0 on success, -1 on failure. On failure, an error message can be retrieved
/// with `rssn_get_last_error`.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_init_plugin_manager(
    plugin_dir_ptr: *const c_char
) -> i32 {

    unsafe {

        let handle_error =
            |err_msg: String| {

                update_last_error(
                    err_msg,
                );

                -1
            };

        let plugin_dir =
            match CStr::from_ptr(
                plugin_dir_ptr,
            )
            .to_str()
            {
                | Ok(s) => s,
                | Err(e) => {
                    return handle_error(
                        format!(
                    "Invalid UTF-8 in \
                     plugin_dir: {e}"
                ),
                    );
                },
            };

        match PluginManager::new(
            plugin_dir,
        ) {
            | Ok(manager) => {

                *PLUGIN_MANAGER.lock().expect("Plugin Manager error") = Some(manager);

                0
            },
            | Err(e) => {
                handle_error(format!(
                "Failed to initialize \
                 PluginManager: {e}"
            ))
            },
        }
    }
}

/// Executes a command on a loaded plugin.
///
/// # Arguments
/// * `plugin_name_ptr` - A null-terminated UTF-8 string representing the plugin's name.
/// * `command_ptr` - A null-terminated UTF-8 string for the command to execute.
/// * `args_handle` - A handle to the `Expr` object to be passed as an argument.
///
/// # Returns
/// A handle to the resulting `Expr` object on success, or 0 on failure.
/// On failure, an error message can be retrieved with `rssn_get_last_error`.
#[unsafe(no_mangle)]
#[deprecated(
    since = "0.1.19",
    note = "This module is deprecated \
            and please use special \
            ffi api modules instead."
)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_plugin_execute(
    plugin_name_ptr: *const c_char,
    command_ptr: *const c_char,
    args_handle: usize,
) -> usize {

    unsafe {

        let handle_error =
            |err_msg: String| {

                update_last_error(
                    err_msg,
                );

                0
            };

        let plugin_name =
            match CStr::from_ptr(
                plugin_name_ptr,
            )
            .to_str()
            {
                | Ok(s) => s,
                | Err(e) => {

                    return handle_error(format!("Invalid UTF-8 in plugin_name: {e}"));
                },
            };

        let command =
            match CStr::from_ptr(
                command_ptr,
            )
            .to_str()
            {
                | Ok(s) => s,
                | Err(e) => {

                    return handle_error(format!("Invalid UTF-8 in command: {e}"));
                },
            };

        let args_expr =
            match HANDLE_MANAGER
                .get(args_handle)
            {
                | Some(expr) => expr,
                | None => {

                    return handle_error(format!("Invalid handle for args: {args_handle}"));
                },
            };

        let result =
            match &*PLUGIN_MANAGER
                .lock()
                .expect(
                    "Plugin Manager \
                     Error",
                ) {
                | Some(pm) => {
                    pm.execute_plugin(
                        plugin_name,
                        command,
                        &args_expr,
                    )
                },
                | None => {

                    return handle_error(
                "Plugin manager not initialized. Call rssn_init_plugin_manager first.".to_string(),
            );
                },
            };

        match result {
            | Ok(result_expr) => {
                HANDLE_MANAGER
                    .insert(result_expr)
            },
            | Err(e) => {
                handle_error(format!(
                    "Plugin execution \
                     failed for \
                     '{plugin_name}': \
                     {e}"
                ))
            },
        }
    }
}
