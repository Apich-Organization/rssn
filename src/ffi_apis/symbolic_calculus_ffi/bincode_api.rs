//! Bincode-based FFI API for symbolic calculus functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::symbolic::calculus;
use crate::symbolic::core::Expr;

/// Differentiates an expression using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_differentiate(
    expr_buf: BincodeBuffer,
    var: *const c_char,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (expr, var_str)
    { (Some(e), Some(v)) => {

        let result =
            calculus::differentiate(
                &e, v,
            );

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Integrates an expression using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_integrate(
    expr_buf: BincodeBuffer,
    var: *const c_char,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (expr, var_str)
    { (Some(e), Some(v)) => {

        let result =
            calculus::integrate(
                &e, v, None, None,
            );

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Computes definite integral using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_definite_integrate(
    expr_buf: BincodeBuffer,
    var: *const c_char,
    lower_buf: BincodeBuffer,
    upper_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let lower: Option<Expr> =
        from_bincode_buffer(&lower_buf);

    let upper: Option<Expr> =
        from_bincode_buffer(&upper_buf);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (
        expr,
        var_str,
        lower,
        upper,
    ) { (
        Some(e),
        Some(v),
        Some(l),
        Some(u),
    ) => {

        let result = calculus::definite_integrate(&e, v, &l, &u);

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Checks analytic using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_check_analytic(
    expr_buf: BincodeBuffer,
    var: *const c_char,
) -> bool {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (expr, var_str)
    { (Some(e), Some(v)) => {

        calculus::check_analytic(&e, v)
    } _ => {

        false
    }}
}

/// Computes limit using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_limit(
    expr_buf: BincodeBuffer,
    var: *const c_char,
    point_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let point: Option<Expr> =
        from_bincode_buffer(&point_buf);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (expr, var_str, point)
    { (Some(e), Some(v), Some(p)) => {

        let result =
            calculus::limit(&e, v, &p);

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Finds poles using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_find_poles(
    expr_buf: BincodeBuffer,
    var: *const c_char,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (expr, var_str)
    { (Some(e), Some(v)) => {

        let result =
            calculus::find_poles(&e, v);

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Calculates residue using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_calculate_residue(
    expr_buf: BincodeBuffer,
    var: *const c_char,
    pole_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let pole: Option<Expr> =
        from_bincode_buffer(&pole_buf);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (expr, var_str, pole)
    { (Some(e), Some(v), Some(p)) => {

        let result =
            calculus::calculate_residue(
                &e, v, &p,
            );

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Finds pole order using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_find_pole_order(
    expr_buf: BincodeBuffer,
    var: *const c_char,
    pole_buf: BincodeBuffer,
) -> usize {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let pole: Option<Expr> =
        from_bincode_buffer(&pole_buf);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (expr, var_str, pole)
    { (Some(e), Some(v), Some(p)) => {

        calculus::find_pole_order(
            &e, v, &p,
        )
    } _ => {

        0
    }}
}

/// Substitutes using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_substitute(
    expr_buf: BincodeBuffer,
    var: *const c_char,
    replacement_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let replacement: Option<Expr> =
        from_bincode_buffer(
            &replacement_buf,
        );

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (
        expr,
        var_str,
        replacement,
    ) { (Some(e), Some(v), Some(r)) => {

        let result =
            calculus::substitute(
                &e, v, &r,
            );

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Gets real and imaginary parts using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_get_real_imag_parts(
    expr_buf: BincodeBuffer
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    if let Some(e) = expr {

        let (re, im) = calculus::get_real_imag_parts(&e);

        to_bincode_buffer(&vec![re, im])
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes path integral using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_path_integrate(
    expr_buf: BincodeBuffer,
    var: *const c_char,
    contour_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let contour: Option<Expr> =
        from_bincode_buffer(
            &contour_buf,
        );

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (
        expr,
        var_str,
        contour,
    ) { (Some(e), Some(v), Some(c)) => {

        let result =
            calculus::path_integrate(
                &e, v, &c,
            );

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Evaluates at point using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_evaluate_at_point(
    expr_buf: BincodeBuffer,
    var: *const c_char,
    value_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let value: Option<Expr> =
        from_bincode_buffer(&value_buf);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .ok()
        }
    };

    match (expr, var_str, value)
    { (
        Some(e),
        Some(v),
        Some(val),
    ) => {

        let result =
            calculus::evaluate_at_point(
                &e, v, &val,
            );

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}
