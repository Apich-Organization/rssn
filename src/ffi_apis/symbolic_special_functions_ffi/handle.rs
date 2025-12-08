use crate::symbolic::core::Expr;
use crate::symbolic::special_functions;
use crate::ffi_apis::common::*;
use std::os::raw::c_char;

#[no_mangle]
pub unsafe extern "C" fn rssn_gamma(arg: *const Expr) -> *mut Expr {
    if arg.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(special_functions::gamma((*arg).clone())))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_beta(a: *const Expr, b: *const Expr) -> *mut Expr {
    if a.is_null() || b.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(special_functions::beta((*a).clone(), (*b).clone())))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_erf(arg: *const Expr) -> *mut Expr {
    if arg.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(special_functions::erf((*arg).clone())))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_erfc(arg: *const Expr) -> *mut Expr {
    if arg.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(special_functions::erfc((*arg).clone())))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_erfi(arg: *const Expr) -> *mut Expr {
    if arg.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(special_functions::erfi((*arg).clone())))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_zeta(arg: *const Expr) -> *mut Expr {
    if arg.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(special_functions::zeta((*arg).clone())))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_digamma(arg: *const Expr) -> *mut Expr {
    if arg.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(special_functions::digamma((*arg).clone())))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_bessel_j(order: *const Expr, arg: *const Expr) -> *mut Expr {
    if order.is_null() || arg.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(special_functions::bessel_j((*order).clone(), (*arg).clone())))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_bessel_y(order: *const Expr, arg: *const Expr) -> *mut Expr {
    if order.is_null() || arg.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(special_functions::bessel_y((*order).clone(), (*arg).clone())))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_legendre_p(degree: *const Expr, arg: *const Expr) -> *mut Expr {
    if degree.is_null() || arg.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(special_functions::legendre_p((*degree).clone(), (*arg).clone())))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_laguerre_l(degree: *const Expr, arg: *const Expr) -> *mut Expr {
    if degree.is_null() || arg.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(special_functions::laguerre_l((*degree).clone(), (*arg).clone())))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_hermite_h(degree: *const Expr, arg: *const Expr) -> *mut Expr {
    if degree.is_null() || arg.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(special_functions::hermite_h((*degree).clone(), (*arg).clone())))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_bessel_differential_equation(
    y: *const Expr, 
    x: *const Expr, 
    n: *const Expr
) -> *mut Expr {
     if y.is_null() || x.is_null() || n.is_null() { return std::ptr::null_mut(); }
     Box::into_raw(Box::new(special_functions::bessel_differential_equation(&*y, &*x, &*n)))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_legendre_differential_equation(
    y: *const Expr, 
    x: *const Expr, 
    n: *const Expr
) -> *mut Expr {
     if y.is_null() || x.is_null() || n.is_null() { return std::ptr::null_mut(); }
     Box::into_raw(Box::new(special_functions::legendre_differential_equation(&*y, &*x, &*n)))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_laguerre_differential_equation(
    y: *const Expr, 
    x: *const Expr, 
    n: *const Expr
) -> *mut Expr {
     if y.is_null() || x.is_null() || n.is_null() { return std::ptr::null_mut(); }
     Box::into_raw(Box::new(special_functions::laguerre_differential_equation(&*y, &*x, &*n)))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_hermite_differential_equation(
    y: *const Expr, 
    x: *const Expr, 
    n: *const Expr
) -> *mut Expr {
     if y.is_null() || x.is_null() || n.is_null() { return std::ptr::null_mut(); }
     Box::into_raw(Box::new(special_functions::hermite_differential_equation(&*y, &*x, &*n)))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_legendre_rodrigues_formula(
    n: *const Expr,
    x: *const Expr
) -> *mut Expr {
    if n.is_null() || x.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(special_functions::legendre_rodrigues_formula(&*n, &*x)))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_hermite_rodrigues_formula(
    n: *const Expr,
    x: *const Expr
) -> *mut Expr {
    if n.is_null() || x.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(special_functions::hermite_rodrigues_formula(&*n, &*x)))
}
