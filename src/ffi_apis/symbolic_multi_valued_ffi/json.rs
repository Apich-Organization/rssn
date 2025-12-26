use std::ffi::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::multi_valued::*;

/// Computes general multi-valued logarithm (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_general_log(
    z_json: *const c_char,
    k_json: *const c_char,
) -> *mut c_char {

    let z: Option<Expr> =
        from_json_string(z_json);

    let k: Option<Expr> =
        from_json_string(k_json);

    if let (
        Some(z_expr),
        Some(k_expr),
    ) = (z, k)
    {

        let result = general_log(
            &z_expr,
            &k_expr,
        );

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes general multi-valued square root (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_general_sqrt(
    z_json: *const c_char,
    k_json: *const c_char,
) -> *mut c_char {

    let z: Option<Expr> =
        from_json_string(z_json);

    let k: Option<Expr> =
        from_json_string(k_json);

    if let (
        Some(z_expr),
        Some(k_expr),
    ) = (z, k)
    {

        let result = general_sqrt(
            &z_expr,
            &k_expr,
        );

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes general multi-valued power (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_general_power(
    z_json: *const c_char,
    w_json: *const c_char,
    k_json: *const c_char,
) -> *mut c_char {

    let z: Option<Expr> =
        from_json_string(z_json);

    let w: Option<Expr> =
        from_json_string(w_json);

    let k: Option<Expr> =
        from_json_string(k_json);

    if let (
        Some(z_expr),
        Some(w_expr),
        Some(k_expr),
    ) = (z, w, k)
    {

        let result = general_power(
            &z_expr,
            &w_expr,
            &k_expr,
        );

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes general multi-valued n-th root (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_general_nth_root(
    z_json: *const c_char,
    n_json: *const c_char,
    k_json: *const c_char,
) -> *mut c_char {

    let z: Option<Expr> =
        from_json_string(z_json);

    let n: Option<Expr> =
        from_json_string(n_json);

    let k: Option<Expr> =
        from_json_string(k_json);

    if let (
        Some(z_expr),
        Some(n_expr),
        Some(k_expr),
    ) = (z, n, k)
    {

        let result = general_nth_root(
            &z_expr,
            &n_expr,
            &k_expr,
        );

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes general multi-valued arcsin (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_general_arcsin(
    z_json: *const c_char,
    k_json: *const c_char,
) -> *mut c_char {

    let z: Option<Expr> =
        from_json_string(z_json);

    let k: Option<Expr> =
        from_json_string(k_json);

    if let (
        Some(z_expr),
        Some(k_expr),
    ) = (z, k)
    {

        let result = general_arcsin(
            &z_expr,
            &k_expr,
        );

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes general multi-valued arccos (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_general_arccos(
    z_json: *const c_char,
    k_json: *const c_char,
    s_json: *const c_char,
) -> *mut c_char {

    let z: Option<Expr> =
        from_json_string(z_json);

    let k: Option<Expr> =
        from_json_string(k_json);

    let s: Option<Expr> =
        from_json_string(s_json);

    if let (
        Some(z_expr),
        Some(k_expr),
        Some(s_expr),
    ) = (z, k, s)
    {

        let result = general_arccos(
            &z_expr,
            &k_expr,
            &s_expr,
        );

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes general multi-valued arctan (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_general_arctan(
    z_json: *const c_char,
    k_json: *const c_char,
) -> *mut c_char {

    let z: Option<Expr> =
        from_json_string(z_json);

    let k: Option<Expr> =
        from_json_string(k_json);

    if let (
        Some(z_expr),
        Some(k_expr),
    ) = (z, k)
    {

        let result = general_arctan(
            &z_expr,
            &k_expr,
        );

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes argument (angle) of complex number (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_arg(
    z_json: *const c_char
) -> *mut c_char {

    let z: Option<Expr> =
        from_json_string(z_json);

    if let Some(z_expr) = z {

        let result = arg(&z_expr);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes absolute value (magnitude) of complex number (JSON)
#[no_mangle]

pub extern "C" fn rssn_json_abs(
    z_json: *const c_char
) -> *mut c_char {

    let z: Option<Expr> =
        from_json_string(z_json);

    if let Some(z_expr) = z {

        let result = abs(&z_expr);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}
