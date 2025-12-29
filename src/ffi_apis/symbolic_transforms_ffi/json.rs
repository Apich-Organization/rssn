//! JSON-based FFI API for symbolic integral transforms.
//!
//! This module provides JSON string-based FFI functions for Fourier, Laplace, and Z-transforms,
//! allowing for language-agnostic integration through JSON serialization/deserialization.

use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::transforms;

// --- Fourier Transform ---

/// Computes the Fourier transform of an expression.

///

/// Takes JSON strings representing `Expr` (expression), `String` (input variable), and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the Fourier transform.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_fourier_transform(
    expr_json: *const c_char,
    in_var_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let in_var: Option<String> =
        from_json_string(in_var_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (
        Some(e),
        Some(iv),
        Some(ov),
    ) = (
        expr,
        in_var,
        out_var,
    ) {

        to_json_string(&transforms::fourier_transform(&e, &iv, &ov))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the inverse Fourier transform of an expression.

///

/// Takes JSON strings representing `Expr` (expression), `String` (input variable), and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the inverse Fourier transform.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_inverse_fourier_transform(
    expr_json: *const c_char,
    in_var_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let in_var: Option<String> =
        from_json_string(in_var_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (
        Some(e),
        Some(iv),
        Some(ov),
    ) = (
        expr,
        in_var,
        out_var,
    ) {

        to_json_string(&transforms::inverse_fourier_transform(&e, &iv, &ov))
    } else {

        std::ptr::null_mut()
    }
}

/// Applies the time shift property of the Fourier transform.

///

/// Takes JSON strings representing `Expr` (frequency domain expression), `Expr` (time shift amount `a`),

/// and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the transformed expression.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_fourier_time_shift(
    f_omega_json: *const c_char,
    a_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let f: Option<Expr> =
        from_json_string(f_omega_json);

    let a: Option<Expr> =
        from_json_string(a_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (
        Some(f),
        Some(a),
        Some(ov),
    ) = (f, a, out_var)
    {

        to_json_string(&transforms::fourier_time_shift(&f, &a, &ov))
    } else {

        std::ptr::null_mut()
    }
}

/// Applies the frequency shift property of the Fourier transform.

///

/// Takes JSON strings representing `Expr` (frequency domain expression), `Expr` (frequency shift amount `a`),

/// and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the transformed expression.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_fourier_frequency_shift(
    f_omega_json: *const c_char,
    a_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let f: Option<Expr> =
        from_json_string(f_omega_json);

    let a: Option<Expr> =
        from_json_string(a_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (
        Some(f),
        Some(a),
        Some(ov),
    ) = (f, a, out_var)
    {

        to_json_string(&transforms::fourier_frequency_shift(&f, &a, &ov))
    } else {

        std::ptr::null_mut()
    }
}

/// Applies the scaling property of the Fourier transform.

///

/// Takes JSON strings representing `Expr` (frequency domain expression), `Expr` (scaling factor `a`),

/// and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the transformed expression.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_fourier_scaling(
    f_omega_json: *const c_char,
    a_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let f: Option<Expr> =
        from_json_string(f_omega_json);

    let a: Option<Expr> =
        from_json_string(a_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (
        Some(f),
        Some(a),
        Some(ov),
    ) = (f, a, out_var)
    {

        to_json_string(&transforms::fourier_scaling(&f, &a, &ov))
    } else {

        std::ptr::null_mut()
    }
}

/// Applies the differentiation property of the Fourier transform.

///

/// Takes JSON strings representing `Expr` (frequency domain expression) and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the transformed expression.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_fourier_differentiation(
    f_omega_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let f: Option<Expr> =
        from_json_string(f_omega_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (Some(f), Some(ov)) =
        (f, out_var)
    {

        to_json_string(&transforms::fourier_differentiation(&f, &ov))
    } else {

        std::ptr::null_mut()
    }
}

// --- Laplace Transform ---

/// Computes the Laplace transform of an expression.

///

/// Takes JSON strings representing `Expr` (expression), `String` (input variable), and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the Laplace transform.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_laplace_transform(
    expr_json: *const c_char,
    in_var_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let in_var: Option<String> =
        from_json_string(in_var_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (
        Some(e),
        Some(iv),
        Some(ov),
    ) = (
        expr,
        in_var,
        out_var,
    ) {

        to_json_string(&transforms::laplace_transform(&e, &iv, &ov))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the inverse Laplace transform of an expression.

///

/// Takes JSON strings representing `Expr` (expression), `String` (input variable), and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the inverse Laplace transform.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_inverse_laplace_transform(
    expr_json: *const c_char,
    in_var_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let in_var: Option<String> =
        from_json_string(in_var_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (
        Some(e),
        Some(iv),
        Some(ov),
    ) = (
        expr,
        in_var,
        out_var,
    ) {

        to_json_string(&transforms::inverse_laplace_transform(&e, &iv, &ov))
    } else {

        std::ptr::null_mut()
    }
}

/// Applies the time shift property of the Laplace transform.

///

/// Takes JSON strings representing `Expr` (s-domain expression), `Expr` (time shift amount `a`),

/// and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the transformed expression.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_laplace_time_shift(
    f_s_json: *const c_char,
    a_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let f: Option<Expr> =
        from_json_string(f_s_json);

    let a: Option<Expr> =
        from_json_string(a_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (
        Some(f),
        Some(a),
        Some(ov),
    ) = (f, a, out_var)
    {

        to_json_string(&transforms::laplace_time_shift(&f, &a, &ov))
    } else {

        std::ptr::null_mut()
    }
}

/// Applies the frequency shift property of the Laplace transform.

///

/// Takes JSON strings representing `Expr` (s-domain expression), `Expr` (frequency shift amount `a`),

/// and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the transformed expression.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_laplace_frequency_shift(
    f_s_json: *const c_char,
    a_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let f: Option<Expr> =
        from_json_string(f_s_json);

    let a: Option<Expr> =
        from_json_string(a_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (
        Some(f),
        Some(a),
        Some(ov),
    ) = (f, a, out_var)
    {

        to_json_string(&transforms::laplace_frequency_shift(&f, &a, &ov))
    } else {

        std::ptr::null_mut()
    }
}

/// Applies the scaling property of the Laplace transform.

///

/// Takes JSON strings representing `Expr` (s-domain expression), `Expr` (scaling factor `a`),

/// and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the transformed expression.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_laplace_scaling(
    f_s_json: *const c_char,
    a_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let f: Option<Expr> =
        from_json_string(f_s_json);

    let a: Option<Expr> =
        from_json_string(a_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (
        Some(f),
        Some(a),
        Some(ov),
    ) = (f, a, out_var)
    {

        to_json_string(&transforms::laplace_scaling(&f, &a, &ov))
    } else {

        std::ptr::null_mut()
    }
}

/// Applies the differentiation property of the Laplace transform.

///

/// Takes JSON strings representing `Expr` (s-domain expression), `String` (output variable),

/// and `Expr` (`f(0)` - initial condition).

/// Returns a JSON string representing the `Expr` of the transformed expression.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_laplace_differentiation(
    f_s_json: *const c_char,
    out_var_json: *const c_char,
    f_zero_json: *const c_char,
) -> *mut c_char {

    let f: Option<Expr> =
        from_json_string(f_s_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    let f_zero: Option<Expr> =
        from_json_string(f_zero_json);

    if let (
        Some(f),
        Some(ov),
        Some(fz),
    ) = (f, out_var, f_zero)
    {

        to_json_string(&transforms::laplace_differentiation(&f, &ov, &fz))
    } else {

        std::ptr::null_mut()
    }
}

/// Applies the integration property of the Laplace transform.

///

/// Takes JSON strings representing `Expr` (s-domain expression) and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the transformed expression.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_laplace_integration(
    f_s_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let f: Option<Expr> =
        from_json_string(f_s_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (Some(f), Some(ov)) =
        (f, out_var)
    {

        to_json_string(&transforms::laplace_integration(&f, &ov))
    } else {

        std::ptr::null_mut()
    }
}

// --- Z-Transform ---

/// Computes the Z-transform of an expression.

///

/// Takes JSON strings representing `Expr` (expression), `String` (input variable), and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the Z-transform.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_z_transform(
    expr_json: *const c_char,
    in_var_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let in_var: Option<String> =
        from_json_string(in_var_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (
        Some(e),
        Some(iv),
        Some(ov),
    ) = (
        expr,
        in_var,
        out_var,
    ) {

        to_json_string(
            &transforms::z_transform(
                &e, &iv, &ov,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the inverse Z-transform of an expression.

///

/// Takes JSON strings representing `Expr` (expression), `String` (input variable), and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the inverse Z-transform.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_inverse_z_transform(
    expr_json: *const c_char,
    in_var_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let in_var: Option<String> =
        from_json_string(in_var_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (
        Some(e),
        Some(iv),
        Some(ov),
    ) = (
        expr,
        in_var,
        out_var,
    ) {

        to_json_string(&transforms::inverse_z_transform(&e, &iv, &ov))
    } else {

        std::ptr::null_mut()
    }
}

/// Applies the time shift property of the Z-transform.

///

/// Takes JSON strings representing `Expr` (z-domain expression), `Expr` (time shift amount `k`),

/// and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the transformed expression.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_z_time_shift(
    f_z_json: *const c_char,
    k_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let f: Option<Expr> =
        from_json_string(f_z_json);

    let k: Option<Expr> =
        from_json_string(k_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (
        Some(f),
        Some(k),
        Some(ov),
    ) = (f, k, out_var)
    {

        to_json_string(
            &transforms::z_time_shift(
                &f, &k, &ov,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Applies the scaling property of the Z-transform.

///

/// Takes JSON strings representing `Expr` (z-domain expression), `Expr` (scaling factor `a`),

/// and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the transformed expression.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_z_scaling(
    f_z_json: *const c_char,
    a_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let f: Option<Expr> =
        from_json_string(f_z_json);

    let a: Option<Expr> =
        from_json_string(a_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (
        Some(f),
        Some(a),
        Some(ov),
    ) = (f, a, out_var)
    {

        to_json_string(
            &transforms::z_scaling(
                &f, &a, &ov,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Applies the differentiation property of the Z-transform.

///

/// Takes JSON strings representing `Expr` (z-domain expression) and `String` (output variable).

/// Returns a JSON string representing the `Expr` of the transformed expression.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_z_differentiation(
    f_z_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let f: Option<Expr> =
        from_json_string(f_z_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (Some(f), Some(ov)) =
        (f, out_var)
    {

        to_json_string(&transforms::z_differentiation(&f, &ov))
    } else {

        std::ptr::null_mut()
    }
}

// --- Utils ---

/// Computes the convolution of two functions using the Fourier transform property.

///

/// Takes JSON strings representing `Expr` (two functions `f` and `g`),

/// and `String` (input variable), `String` (output variable).

/// Returns a JSON string representing the `Expr` of the convolution.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_convolution_fourier(
    f_json: *const c_char,
    g_json: *const c_char,
    in_var_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let f: Option<Expr> =
        from_json_string(f_json);

    let g: Option<Expr> =
        from_json_string(g_json);

    let in_var: Option<String> =
        from_json_string(in_var_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (
        Some(f),
        Some(g),
        Some(iv),
        Some(ov),
    ) = (
        f,
        g,
        in_var,
        out_var,
    ) {

        to_json_string(&transforms::convolution_fourier(&f, &g, &iv, &ov))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the convolution of two functions using the Laplace transform property.

///

/// Takes JSON strings representing `Expr` (two functions `f` and `g`),

/// and `String` (input variable), `String` (output variable).

/// Returns a JSON string representing the `Expr` of the convolution.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_convolution_laplace(
    f_json: *const c_char,
    g_json: *const c_char,
    in_var_json: *const c_char,
    out_var_json: *const c_char,
) -> *mut c_char {

    let f: Option<Expr> =
        from_json_string(f_json);

    let g: Option<Expr> =
        from_json_string(g_json);

    let in_var: Option<String> =
        from_json_string(in_var_json);

    let out_var: Option<String> =
        from_json_string(out_var_json);

    if let (
        Some(f),
        Some(g),
        Some(iv),
        Some(ov),
    ) = (
        f,
        g,
        in_var,
        out_var,
    ) {

        to_json_string(&transforms::convolution_laplace(&f, &g, &iv, &ov))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the partial fraction decomposition of an expression.

///

/// Takes JSON strings representing `Expr` (expression) and `String` (variable).

/// Returns a JSON string representing `Vec<Expr>` (partial fraction decomposition).

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_partial_fraction_decomposition(
    expr_json: *const c_char,
    var_json: *const c_char,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    let var: Option<String> =
        from_json_string(var_json);

    if let (Some(expr), Some(var)) =
        (expr, var)
    {

        if let Some(result) = transforms::partial_fraction_decomposition(&expr, &var) {

            return to_json_string(&result);
        }
    }

    std::ptr::null_mut()
}
