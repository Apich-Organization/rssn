use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::series::{
    taylor_series, laurent_series, fourier_series, summation, product,
    analyze_convergence, asymptotic_expansion, analytic_continuation
};
use std::ffi::c_char;

#[no_mangle]
pub extern "C" fn rssn_json_taylor_series(
    expr_json: *const c_char,
    var_json: *const c_char,
    center_json: *const c_char,
    order_json: *const c_char
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let var: Option<String> = from_json_string(var_json);
    let center: Option<Expr> = from_json_string(center_json);
    let order: Option<usize> = from_json_string(order_json);
    
    if let (Some(e), Some(v), Some(c), Some(o)) = (expr, var, center, order) {
        let result = taylor_series(&e, &v, &c, o);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub extern "C" fn rssn_json_laurent_series(
    expr_json: *const c_char,
    var_json: *const c_char,
    center_json: *const c_char,
    order_json: *const c_char
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let var: Option<String> = from_json_string(var_json);
    let center: Option<Expr> = from_json_string(center_json);
    let order: Option<usize> = from_json_string(order_json);
    
    if let (Some(e), Some(v), Some(c), Some(o)) = (expr, var, center, order) {
        let result = laurent_series(&e, &v, &c, o);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub extern "C" fn rssn_json_fourier_series(
    expr_json: *const c_char,
    var_json: *const c_char,
    period_json: *const c_char,
    order_json: *const c_char
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let var: Option<String> = from_json_string(var_json);
    let period: Option<Expr> = from_json_string(period_json);
    let order: Option<usize> = from_json_string(order_json);
    
    if let (Some(e), Some(v), Some(p), Some(o)) = (expr, var, period, order) {
        let result = fourier_series(&e, &v, &p, o);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub extern "C" fn rssn_json_summation(
    expr_json: *const c_char,
    var_json: *const c_char,
    lower_json: *const c_char,
    upper_json: *const c_char
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let var: Option<String> = from_json_string(var_json);
    let lower: Option<Expr> = from_json_string(lower_json);
    let upper: Option<Expr> = from_json_string(upper_json);
    
    if let (Some(e), Some(v), Some(l), Some(u)) = (expr, var, lower, upper) {
        let result = summation(&e, &v, &l, &u);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub extern "C" fn rssn_json_product(
    expr_json: *const c_char,
    var_json: *const c_char,
    lower_json: *const c_char,
    upper_json: *const c_char
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let var: Option<String> = from_json_string(var_json);
    let lower: Option<Expr> = from_json_string(lower_json);
    let upper: Option<Expr> = from_json_string(upper_json);
    
    if let (Some(e), Some(v), Some(l), Some(u)) = (expr, var, lower, upper) {
        let result = product(&e, &v, &l, &u);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub extern "C" fn rssn_series_json_analyze_convergence(
    series_json: *const c_char,
    var_json: *const c_char
) -> *mut c_char {
    let series: Option<Expr> = from_json_string(series_json);
    let var: Option<String> = from_json_string(var_json);
    
    if let (Some(s), Some(v)) = (series, var) {
        let result = analyze_convergence(&s, &v);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub extern "C" fn rssn_json_asymptotic_expansion(
    expr_json: *const c_char,
    var_json: *const c_char,
    point_json: *const c_char,
    order_json: *const c_char
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let var: Option<String> = from_json_string(var_json);
    let point: Option<Expr> = from_json_string(point_json);
    let order: Option<usize> = from_json_string(order_json);
    
    if let (Some(e), Some(v), Some(p), Some(o)) = (expr, var, point, order) {
        let result = asymptotic_expansion(&e, &v, &p, o);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub extern "C" fn rssn_json_analytic_continuation(
    expr_json: *const c_char,
    var_json: *const c_char,
    orig_center_json: *const c_char,
    new_center_json: *const c_char,
    order_json: *const c_char
) -> *mut c_char {
    let expr: Option<Expr> = from_json_string(expr_json);
    let var: Option<String> = from_json_string(var_json);
    let orig_center: Option<Expr> = from_json_string(orig_center_json);
    let new_center: Option<Expr> = from_json_string(new_center_json);
    let order: Option<usize> = from_json_string(order_json);
    
    if let (Some(e), Some(v), Some(oc), Some(nc), Some(o)) = (expr, var, orig_center, new_center, order) {
        let result = analytic_continuation(&e, &v, &oc, &nc, o);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}
