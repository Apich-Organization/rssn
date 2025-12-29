//! Handle-based FFI API for symbolic integral transforms.
//!
//! This module provides C-compatible FFI functions for Fourier, Laplace, and Z-transforms,
//! including their inverse transforms and related properties (time shift, frequency shift,
//! scaling, differentiation, and convolution theorems).

use std::os::raw::c_char;

use crate::ffi_apis::common::c_str_to_str;
use crate::symbolic::core::Expr;
use crate::symbolic::transforms;

/// Computes the symbolic Fourier transform of an expression.
///
/// # Safety
/// Caller must ensure `expr` is a valid pointer to an `Expr`.
/// `in_var` and `out_var` must be valid C strings or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_fourier_transform(
    expr: *const Expr,
    in_var: *const c_char,
    out_var: *const c_char,
) -> *mut Expr {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let in_v = c_str_to_str(in_var)
        .unwrap_or("t");

    let out_v = c_str_to_str(out_var)
        .unwrap_or("omega");

    Box::into_raw(Box::new(
        transforms::fourier_transform(
            &*expr,
            in_v,
            out_v,
        ),
    ))
}

/// Computes the symbolic inverse Fourier transform of an expression.
///
/// # Safety
/// Caller must ensure `expr` is a valid pointer to an `Expr`.
/// `in_var` and `out_var` must be valid C strings or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_inverse_fourier_transform(
    expr: *const Expr,
    in_var: *const c_char,
    out_var: *const c_char,
) -> *mut Expr {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let in_v = c_str_to_str(in_var)
        .unwrap_or("omega");

    let out_v = c_str_to_str(out_var)
        .unwrap_or("t");

    Box::into_raw(Box::new(
        transforms::inverse_fourier_transform(&*expr, in_v, out_v),
    ))
}

/// Computes the symbolic Laplace transform of an expression.
///
/// # Safety
/// Caller must ensure `expr` is a valid pointer to an `Expr`.
/// `in_var` and `out_var` must be valid C strings or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_laplace_transform(
    expr: *const Expr,
    in_var: *const c_char,
    out_var: *const c_char,
) -> *mut Expr {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let in_v = c_str_to_str(in_var)
        .unwrap_or("t");

    let out_v = c_str_to_str(out_var)
        .unwrap_or("s");

    Box::into_raw(Box::new(
        transforms::laplace_transform(
            &*expr,
            in_v,
            out_v,
        ),
    ))
}

/// Computes the symbolic inverse Laplace transform of an expression.
///
/// # Safety
/// Caller must ensure `expr` is a valid pointer to an `Expr`.
/// `in_var` and `out_var` must be valid C strings or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_inverse_laplace_transform(
    expr: *const Expr,
    in_var: *const c_char,
    out_var: *const c_char,
) -> *mut Expr {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let in_v = c_str_to_str(in_var)
        .unwrap_or("s");

    let out_v = c_str_to_str(out_var)
        .unwrap_or("t");

    Box::into_raw(Box::new(
        transforms::inverse_laplace_transform(&*expr, in_v, out_v),
    ))
}

/// Computes the symbolic Z-transform of an expression.
///
/// # Safety
/// Caller must ensure `expr` is a valid pointer to an `Expr`.
/// `in_var` and `out_var` must be valid C strings or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_z_transform(
    expr: *const Expr,
    in_var: *const c_char,
    out_var: *const c_char,
) -> *mut Expr {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let in_v = c_str_to_str(in_var)
        .unwrap_or("n");

    let out_v = c_str_to_str(out_var)
        .unwrap_or("z");

    Box::into_raw(Box::new(
        transforms::z_transform(
            &*expr,
            in_v,
            out_v,
        ),
    ))
}

/// Computes the symbolic inverse Z-transform of an expression.
///
/// # Safety
/// Caller must ensure `expr` is a valid pointer to an `Expr`.
/// `in_var` and `out_var` must be valid C strings or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_inverse_z_transform(
    expr: *const Expr,
    in_var: *const c_char,
    out_var: *const c_char,
) -> *mut Expr {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let in_v = c_str_to_str(in_var)
        .unwrap_or("z");

    let out_v = c_str_to_str(out_var)
        .unwrap_or("n");

    Box::into_raw(Box::new(
        transforms::inverse_z_transform(
            &*expr,
            in_v,
            out_v,
        ),
    ))
}

/// Applies the time shift property of the Fourier transform.
///
/// # Safety
/// Caller must ensure `f_omega` and `a` are valid pointers to an `Expr`.
/// `out_var` must be a valid C string or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_fourier_time_shift(
    f_omega: *const Expr,
    a: *const Expr,
    out_var: *const c_char,
) -> *mut Expr {

    if f_omega.is_null() || a.is_null()
    {

        return std::ptr::null_mut();
    }

    let out_v = c_str_to_str(out_var)
        .unwrap_or("omega");

    Box::into_raw(Box::new(
        transforms::fourier_time_shift(
            &*f_omega,
            &*a,
            out_v,
        ),
    ))
}

/// Applies the frequency shift property of the Fourier transform.
///
/// # Safety
/// Caller must ensure `f_omega` and `a` are valid pointers to an `Expr`.
/// `out_var` must be a valid C string or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_fourier_frequency_shift(
    f_omega: *const Expr,
    a: *const Expr,
    out_var: *const c_char,
) -> *mut Expr {

    if f_omega.is_null() || a.is_null()
    {

        return std::ptr::null_mut();
    }

    let out_v = c_str_to_str(out_var)
        .unwrap_or("omega");

    Box::into_raw(Box::new(
        transforms::fourier_frequency_shift(
            &*f_omega,
            &*a,
            out_v,
        ),
    ))
}

/// Applies the scaling property of the Fourier transform.
///
/// # Safety
/// Caller must ensure `f_omega` and `a` are valid pointers to an `Expr`.
/// `out_var` must be a valid C string or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_fourier_scaling(
    f_omega: *const Expr,
    a: *const Expr,
    out_var: *const c_char,
) -> *mut Expr {

    if f_omega.is_null() || a.is_null()
    {

        return std::ptr::null_mut();
    }

    let out_v = c_str_to_str(out_var)
        .unwrap_or("omega");

    Box::into_raw(Box::new(
        transforms::fourier_scaling(
            &*f_omega,
            &*a,
            out_v,
        ),
    ))
}

/// Applies the differentiation property of the Fourier transform.
///
/// # Safety
/// Caller must ensure `f_omega` is a valid pointer to an `Expr`.
/// `out_var` must be a valid C string or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_fourier_differentiation(
    f_omega: *const Expr,
    out_var: *const c_char,
) -> *mut Expr {

    if f_omega.is_null() {

        return std::ptr::null_mut();
    }

    let out_v = c_str_to_str(out_var)
        .unwrap_or("omega");

    Box::into_raw(Box::new(
        transforms::fourier_differentiation(&*f_omega, out_v),
    ))
}

/// Applies the time shift property of the Laplace transform.
///
/// # Safety
/// Caller must ensure `f_s` and `a` are valid pointers to an `Expr`.
/// `out_var` must be a valid C string or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_laplace_time_shift(
    f_s: *const Expr,
    a: *const Expr,
    out_var: *const c_char,
) -> *mut Expr {

    if f_s.is_null() || a.is_null() {

        return std::ptr::null_mut();
    }

    let out_v = c_str_to_str(out_var)
        .unwrap_or("s");

    Box::into_raw(Box::new(
        transforms::laplace_time_shift(
            &*f_s, &*a, out_v,
        ),
    ))
}

/// Applies the differentiation property of the Laplace transform.
///
/// # Safety
/// Caller must ensure `f_s` and `f_zero` are valid pointers to an `Expr`.
/// `out_var` must be a valid C string or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_laplace_differentiation(
    f_s: *const Expr,
    out_var: *const c_char,
    f_zero: *const Expr,
) -> *mut Expr {

    if f_s.is_null() || f_zero.is_null()
    {

        return std::ptr::null_mut();
    }

    let out_v = c_str_to_str(out_var)
        .unwrap_or("s");

    Box::into_raw(Box::new(
        transforms::laplace_differentiation(
            &*f_s,
            out_v,
            &*f_zero,
        ),
    ))
}

/// Computes the convolution of two expressions in the Fourier domain.
///
/// # Safety
/// Caller must ensure `f` and `g` are valid pointers to an `Expr`.
/// `in_var` and `out_var` must be valid C strings or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_convolution_fourier(
    f: *const Expr,
    g: *const Expr,
    in_var: *const c_char,
    out_var: *const c_char,
) -> *mut Expr {

    if f.is_null() || g.is_null() {

        return std::ptr::null_mut();
    }

    let in_v = c_str_to_str(in_var)
        .unwrap_or("t");

    let out_v = c_str_to_str(out_var)
        .unwrap_or("omega");

    Box::into_raw(Box::new(
        transforms::convolution_fourier(
            &*f, &*g, in_v, out_v,
        ),
    ))
}

/// Computes the convolution of two expressions in the Laplace domain.
///
/// # Safety
/// Caller must ensure `f` and `g` are valid pointers to an `Expr`.
/// `in_var` and `out_var` must be valid C strings or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_convolution_laplace(
    f: *const Expr,
    g: *const Expr,
    in_var: *const c_char,
    out_var: *const c_char,
) -> *mut Expr {

    if f.is_null() || g.is_null() {

        return std::ptr::null_mut();
    }

    let in_v = c_str_to_str(in_var)
        .unwrap_or("t");

    let out_v = c_str_to_str(out_var)
        .unwrap_or("s");

    Box::into_raw(Box::new(
        transforms::convolution_laplace(
            &*f, &*g, in_v, out_v,
        ),
    ))
}

/// Applies the frequency shift property of the Laplace transform.
///
/// # Safety
/// Caller must ensure `f_s` and `a` are valid pointers to an `Expr`.
/// `out_var` must be a valid C string or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_laplace_frequency_shift(
    f_s: *const Expr,
    a: *const Expr,
    out_var: *const c_char,
) -> *mut Expr {

    if f_s.is_null() || a.is_null() {

        return std::ptr::null_mut();
    }

    let out_v = c_str_to_str(out_var)
        .unwrap_or("s");

    Box::into_raw(Box::new(
        transforms::laplace_frequency_shift(&*f_s, &*a, out_v),
    ))
}

/// Applies the scaling property of the Laplace transform.
///
/// # Safety
/// Caller must ensure `f_s` and `a` are valid pointers to an `Expr`.
/// `out_var` must be a valid C string or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_laplace_scaling(
    f_s: *const Expr,
    a: *const Expr,
    out_var: *const c_char,
) -> *mut Expr {

    if f_s.is_null() || a.is_null() {

        return std::ptr::null_mut();
    }

    let out_v = c_str_to_str(out_var)
        .unwrap_or("s");

    Box::into_raw(Box::new(
        transforms::laplace_scaling(
            &*f_s, &*a, out_v,
        ),
    ))
}

/// Applies the integration property of the Laplace transform.
///
/// # Safety
/// Caller must ensure `f_s` is a valid pointer to an `Expr`.
/// `out_var` must be a valid C string or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_laplace_integration(
    f_s: *const Expr,
    out_var: *const c_char,
) -> *mut Expr {

    if f_s.is_null() {

        return std::ptr::null_mut();
    }

    let out_v = c_str_to_str(out_var)
        .unwrap_or("s");

    Box::into_raw(Box::new(
        transforms::laplace_integration(
            &*f_s, out_v,
        ),
    ))
}

/// Applies the time shift property of the Z-transform.
///
/// # Safety
/// Caller must ensure `f_z` and `k` are valid pointers to an `Expr`.
/// `out_var` must be a valid C string or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_z_time_shift(
    f_z: *const Expr,
    k: *const Expr,
    out_var: *const c_char,
) -> *mut Expr {

    if f_z.is_null() || k.is_null() {

        return std::ptr::null_mut();
    }

    let out_v = c_str_to_str(out_var)
        .unwrap_or("z");

    Box::into_raw(Box::new(
        transforms::z_time_shift(
            &*f_z, &*k, out_v,
        ),
    ))
}

/// Applies the scaling property of the Z-transform.
///
/// # Safety
/// Caller must ensure `f_z` and `a` are valid pointers to an `Expr`.
/// `out_var` must be a valid C string or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_z_scaling(
    f_z: *const Expr,
    a: *const Expr,
    out_var: *const c_char,
) -> *mut Expr {

    if f_z.is_null() || a.is_null() {

        return std::ptr::null_mut();
    }

    let out_v = c_str_to_str(out_var)
        .unwrap_or("z");

    Box::into_raw(Box::new(
        transforms::z_scaling(
            &*f_z, &*a, out_v,
        ),
    ))
}

/// Applies the differentiation property of the Z-transform.
///
/// # Safety
/// Caller must ensure `f_z` is a valid pointer to an `Expr`.
/// `out_var` must be a valid C string or null (defaults apply).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_z_differentiation(
    f_z: *const Expr,
    out_var: *const c_char,
) -> *mut Expr {

    if f_z.is_null() {

        return std::ptr::null_mut();
    }

    let out_v = c_str_to_str(out_var)
        .unwrap_or("z");

    Box::into_raw(Box::new(
        transforms::z_differentiation(
            &*f_z, out_v,
        ),
    ))
}

// --- ExprList Support for Partial Fraction Decomposition ---

/// A C-compatible wrapper for a `Vec<Expr>`, used for returning lists of expressions via FFI.

pub struct ExprList(pub Vec<Expr>);

/// Computes the partial fraction decomposition of an expression.

///

/// Takes a raw pointer to an `Expr` (expression) and a C-style string (variable).

/// Returns a raw pointer to an `ExprList` representing the decomposition.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_partial_fraction_decomposition(
    expr: *const Expr,
    var: *const c_char,
) -> *mut ExprList {

    if expr.is_null() {

        return std::ptr::null_mut();
    }

    let v = c_str_to_str(var)
        .unwrap_or("x");

    if let Some(res) = transforms::partial_fraction_decomposition(&*expr, v) {

        Box::into_raw(Box::new(ExprList(
            res,
        )))
    } else {

        std::ptr::null_mut()
    }
}

/// Returns the length of an `ExprList`.

///

/// Takes a raw pointer to an `ExprList`.

/// Returns a `usize` representing its length.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub const unsafe extern "C" fn rssn_expr_list_len(
    list: *const ExprList
) -> usize {

    if list.is_null() {

        return 0;
    }

    (*list).0.len()
}

/// Returns a specific element from an `ExprList`.

///

/// Takes a raw pointer to an `ExprList` and a `usize` index.

/// Returns a raw pointer to an `Expr` at that index.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_expr_list_get(
    list: *const ExprList,
    index: usize,
) -> *mut Expr {

    if list.is_null() {

        return std::ptr::null_mut();
    }

    if let Some(item) =
        (&(*list).0).get(index)
    {

        Box::into_raw(Box::new(
            item.clone(),
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Frees the memory allocated for an `ExprList`.

///

/// Takes a raw mutable pointer to an `ExprList`.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_expr_list_free(
    list: *mut ExprList
) {

    if !list.is_null() {

        drop(Box::from_raw(list));
    }
}
