use crate::symbolic::convergence::ConvergenceResult;
use crate::symbolic::core::Expr;
use crate::symbolic::series::{
    analytic_continuation, asymptotic_expansion, fourier_series, laurent_series, product,
    summation, taylor_series,
};

#[no_mangle]

pub extern "C" fn rssn_taylor_series_handle(
    expr: *const Expr,
    var: *const std::ffi::c_char,
    center: *const Expr,
    order: usize,
) -> *mut Expr {

    let expr_ref = unsafe {

        &*expr
    };

    let center_ref = unsafe {

        &*center
    };

    let var_str = unsafe {

        if var.is_null() {

            return std::ptr::null_mut();
        }

        std::ffi::CStr::from_ptr(var)
            .to_string_lossy()
            .into_owned()
    };

    let result = taylor_series(expr_ref, &var_str, center_ref, order);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub extern "C" fn rssn_laurent_series_handle(
    expr: *const Expr,
    var: *const std::ffi::c_char,
    center: *const Expr,
    order: usize,
) -> *mut Expr {

    let expr_ref = unsafe {

        &*expr
    };

    let center_ref = unsafe {

        &*center
    };

    let var_str = unsafe {

        if var.is_null() {

            return std::ptr::null_mut();
        }

        std::ffi::CStr::from_ptr(var)
            .to_string_lossy()
            .into_owned()
    };

    let result = laurent_series(expr_ref, &var_str, center_ref, order);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub extern "C" fn rssn_fourier_series_handle(
    expr: *const Expr,
    var: *const std::ffi::c_char,
    period: *const Expr,
    order: usize,
) -> *mut Expr {

    let expr_ref = unsafe {

        &*expr
    };

    let period_ref = unsafe {

        &*period
    };

    let var_str = unsafe {

        if var.is_null() {

            return std::ptr::null_mut();
        }

        std::ffi::CStr::from_ptr(var)
            .to_string_lossy()
            .into_owned()
    };

    let result = fourier_series(expr_ref, &var_str, period_ref, order);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub extern "C" fn rssn_summation_handle(
    expr: *const Expr,
    var: *const std::ffi::c_char,
    lower: *const Expr,
    upper: *const Expr,
) -> *mut Expr {

    let expr_ref = unsafe {

        &*expr
    };

    let lower_ref = unsafe {

        &*lower
    };

    let upper_ref = unsafe {

        &*upper
    };

    let var_str = unsafe {

        if var.is_null() {

            return std::ptr::null_mut();
        }

        std::ffi::CStr::from_ptr(var)
            .to_string_lossy()
            .into_owned()
    };

    let result = summation(expr_ref, &var_str, lower_ref, upper_ref);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub extern "C" fn rssn_product_handle(
    expr: *const Expr,
    var: *const std::ffi::c_char,
    lower: *const Expr,
    upper: *const Expr,
) -> *mut Expr {

    let expr_ref = unsafe {

        &*expr
    };

    let lower_ref = unsafe {

        &*lower
    };

    let upper_ref = unsafe {

        &*upper
    };

    let var_str = unsafe {

        if var.is_null() {

            return std::ptr::null_mut();
        }

        std::ffi::CStr::from_ptr(var)
            .to_string_lossy()
            .into_owned()
    };

    let result = product(expr_ref, &var_str, lower_ref, upper_ref);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub extern "C" fn rssn_series_analyze_convergence_handle(
    series: *const Expr,
    var: *const std::ffi::c_char,
) -> *mut ConvergenceResult {

    let series_ref = unsafe {

        &*series
    };

    let var_str = unsafe {

        if var.is_null() {

            return std::ptr::null_mut();
        }

        std::ffi::CStr::from_ptr(var)
            .to_string_lossy()
            .into_owned()
    };

    let result = crate::symbolic::convergence::analyze_convergence(series_ref, &var_str);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub extern "C" fn rssn_asymptotic_expansion_handle(
    expr: *const Expr,
    var: *const std::ffi::c_char,
    point: *const Expr,
    order: usize,
) -> *mut Expr {

    let expr_ref = unsafe {

        &*expr
    };

    let point_ref = unsafe {

        &*point
    };

    let var_str = unsafe {

        if var.is_null() {

            return std::ptr::null_mut();
        }

        std::ffi::CStr::from_ptr(var)
            .to_string_lossy()
            .into_owned()
    };

    let result = asymptotic_expansion(expr_ref, &var_str, point_ref, order);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub extern "C" fn rssn_analytic_continuation_handle(
    expr: *const Expr,
    var: *const std::ffi::c_char,
    orig_center: *const Expr,
    new_center: *const Expr,
    order: usize,
) -> *mut Expr {

    let expr_ref = unsafe {

        &*expr
    };

    let orig_center_ref = unsafe {

        &*orig_center
    };

    let new_center_ref = unsafe {

        &*new_center
    };

    let var_str = unsafe {

        if var.is_null() {

            return std::ptr::null_mut();
        }

        std::ffi::CStr::from_ptr(var)
            .to_string_lossy()
            .into_owned()
    };

    let result = analytic_continuation(expr_ref, &var_str, orig_center_ref, new_center_ref, order);

    Box::into_raw(Box::new(result))
}
