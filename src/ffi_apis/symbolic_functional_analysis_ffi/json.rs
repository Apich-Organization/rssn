use std::os::raw::c_char;

use crate::ffi_apis::common::{from_json_string, to_json_string};
use crate::symbolic::core::Expr;
use crate::symbolic::functional_analysis::{HilbertSpace, inner_product, norm, gram_schmidt};

/// Constructs a Hilbert space from a JSON-encoded description.
///
/// The input string encodes a [`HilbertSpace`] specification (e.g., underlying
/// function space, inner product, and measure), which is deserialized and
/// returned in canonical internal form.
///
/// # Arguments
///
/// * `json_str` - C string pointer containing JSON for a `HilbertSpace` description.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `HilbertSpace`, or null if the
/// input cannot be deserialized.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer and
/// returns ownership of a heap-allocated C string that must be freed by the caller.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_hilbert_space_create(
    json_str: *const c_char
) -> *mut c_char {

    let space : HilbertSpace = match from_json_string(json_str) {
        | Some(s) => s,
        | None => return std::ptr::null_mut(),
    };

    to_json_string(&space)
}

/// Computes the inner product of two functions in a Hilbert space using JSON serialization.
///
/// Given a Hilbert space and two symbolic functions \(f\) and \(g\), this evaluates
/// the inner product \(\langle f, g \rangle\) according to the space's inner
/// product structure.
///
/// # Arguments
///
/// * `space_json` - C string pointer with JSON-encoded [`HilbertSpace`].
/// * `f_json` - C string pointer with JSON-encoded `Expr` for \(f\).
/// * `g_json` - C string pointer with JSON-encoded `Expr` for \(g\).
///
/// # Returns
///
/// A C string pointer containing JSON-encoded symbolic inner product (typically
/// an `Expr`), or null if any input cannot be deserialized or the computation fails.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw C string pointers and
/// returns ownership of a heap-allocated C string that must be freed by the caller.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_inner_product(
    space_json: *const c_char,
    f_json: *const c_char,
    g_json: *const c_char,
) -> *mut c_char {

    let space : HilbertSpace = match from_json_string(space_json) {
        | Some(s) => s,
        | None => return std::ptr::null_mut(),
    };

    let f: Expr = match from_json_string(
        f_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let g: Expr = match from_json_string(
        g_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let result =
        inner_product(&space, &f, &g);

    to_json_string(&result)
}

/// Computes the norm of a function in a Hilbert space using JSON serialization.
///
/// Given a Hilbert space and a symbolic function \(f\), this evaluates the norm
/// \(\|f\| = \sqrt{\langle f, f \rangle}\) induced by the inner product.
///
/// # Arguments
///
/// * `space_json` - C string pointer with JSON-encoded [`HilbertSpace`].
/// * `f_json` - C string pointer with JSON-encoded `Expr` for \(f\).
///
/// # Returns
///
/// A C string pointer containing JSON-encoded symbolic norm (typically an `Expr`),
/// or null if any input cannot be deserialized or the computation fails.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw C string pointers and
/// returns ownership of a heap-allocated C string that must be freed by the caller.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_norm(
    space_json: *const c_char,
    f_json: *const c_char,
) -> *mut c_char {

    let space : HilbertSpace = match from_json_string(space_json) {
        | Some(s) => s,
        | None => return std::ptr::null_mut(),
    };

    let f: Expr = match from_json_string(
        f_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let result = norm(&space, &f);

    to_json_string(&result)
}

/// Applies the Gram–Schmidt process to produce an orthonormal basis in a Hilbert space using JSON serialization.
///
/// Given a Hilbert space and a list of symbolic basis vectors, this performs the
/// Gram–Schmidt orthonormalization procedure to obtain an orthonormal basis.
///
/// # Arguments
///
/// * `space_json` - C string pointer with JSON-encoded [`HilbertSpace`].
/// * `basis_json` - C string pointer with JSON-encoded `Vec<Expr>` of basis vectors.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Vec<Expr>` for the orthonormal
/// basis, or null if any input cannot be deserialized or the computation fails.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw C string pointers and
/// returns ownership of a heap-allocated C string that must be freed by the caller.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_gram_schmidt(
    space_json: *const c_char,
    basis_json: *const c_char,
) -> *mut c_char {

    let space : HilbertSpace = match from_json_string(space_json) {
        | Some(s) => s,
        | None => return std::ptr::null_mut(),
    };

    let basis : Vec<Expr> = match from_json_string(basis_json) {
        | Some(b) => b,
        | None => return std::ptr::null_mut(),
    };

    let result =
        gram_schmidt(&space, &basis);

    to_json_string(&result)
}
