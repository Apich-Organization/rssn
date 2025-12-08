use crate::symbolic::core::Expr;
use crate::symbolic::special_functions;
use crate::ffi_apis::common::*;
use std::os::raw::c_char;

#[no_mangle]
pub unsafe extern "C" fn rssn_json_gamma(arg_json: *const c_char) -> *mut c_char {
    let arg: Option<Expr> = from_json_string(arg_json);
    if let Some(a) = arg {
        to_json_string(&special_functions::gamma(a))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_beta(a_json: *const c_char, b_json: *const c_char) -> *mut c_char {
    let a: Option<Expr> = from_json_string(a_json);
    let b: Option<Expr> = from_json_string(b_json);
    if let (Some(val_a), Some(val_b)) = (a, b) {
        to_json_string(&special_functions::beta(val_a, val_b))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_erf(arg_json: *const c_char) -> *mut c_char {
    let arg: Option<Expr> = from_json_string(arg_json);
    if let Some(a) = arg {
        to_json_string(&special_functions::erf(a))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_erfc(arg_json: *const c_char) -> *mut c_char {
    let arg: Option<Expr> = from_json_string(arg_json);
    if let Some(a) = arg {
        to_json_string(&special_functions::erfc(a))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_erfi(arg_json: *const c_char) -> *mut c_char {
    let arg: Option<Expr> = from_json_string(arg_json);
    if let Some(a) = arg {
        to_json_string(&special_functions::erfi(a))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_zeta(arg_json: *const c_char) -> *mut c_char {
    let arg: Option<Expr> = from_json_string(arg_json);
    if let Some(a) = arg {
        to_json_string(&special_functions::zeta(a))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_digamma(arg_json: *const c_char) -> *mut c_char {
    let arg: Option<Expr> = from_json_string(arg_json);
    if let Some(a) = arg {
        to_json_string(&special_functions::digamma(a))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_bessel_j(order_json: *const c_char, arg_json: *const c_char) -> *mut c_char {
    let order: Option<Expr> = from_json_string(order_json);
    let arg: Option<Expr> = from_json_string(arg_json);
    if let (Some(o), Some(a)) = (order, arg) {
        to_json_string(&special_functions::bessel_j(o, a))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_bessel_y(order_json: *const c_char, arg_json: *const c_char) -> *mut c_char {
    let order: Option<Expr> = from_json_string(order_json);
    let arg: Option<Expr> = from_json_string(arg_json);
    if let (Some(o), Some(a)) = (order, arg) {
        to_json_string(&special_functions::bessel_y(o, a))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_legendre_p(degree_json: *const c_char, arg_json: *const c_char) -> *mut c_char {
    let degree: Option<Expr> = from_json_string(degree_json);
    let arg: Option<Expr> = from_json_string(arg_json);
    if let (Some(d), Some(a)) = (degree, arg) {
        to_json_string(&special_functions::legendre_p(d, a))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_laguerre_l(degree_json: *const c_char, arg_json: *const c_char) -> *mut c_char {
    let degree: Option<Expr> = from_json_string(degree_json);
    let arg: Option<Expr> = from_json_string(arg_json);
    if let (Some(d), Some(a)) = (degree, arg) {
        to_json_string(&special_functions::laguerre_l(d, a))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_hermite_h(degree_json: *const c_char, arg_json: *const c_char) -> *mut c_char {
    let degree: Option<Expr> = from_json_string(degree_json);
    let arg: Option<Expr> = from_json_string(arg_json);
    if let (Some(d), Some(a)) = (degree, arg) {
        to_json_string(&special_functions::hermite_h(d, a))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_bessel_differential_equation(y_json: *const c_char, x_json: *const c_char, n_json: *const c_char) -> *mut c_char {
    let y: Option<Expr> = from_json_string(y_json);
    let x: Option<Expr> = from_json_string(x_json);
    let n: Option<Expr> = from_json_string(n_json);
    if let (Some(y), Some(x), Some(n)) = (y, x, n) {
        to_json_string(&special_functions::bessel_differential_equation(&y, &x, &n))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_legendre_differential_equation(y_json: *const c_char, x_json: *const c_char, n_json: *const c_char) -> *mut c_char {
    let y: Option<Expr> = from_json_string(y_json);
    let x: Option<Expr> = from_json_string(x_json);
    let n: Option<Expr> = from_json_string(n_json);
    if let (Some(y), Some(x), Some(n)) = (y, x, n) {
        to_json_string(&special_functions::legendre_differential_equation(&y, &x, &n))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_laguerre_differential_equation(y_json: *const c_char, x_json: *const c_char, n_json: *const c_char) -> *mut c_char {
    let y: Option<Expr> = from_json_string(y_json);
    let x: Option<Expr> = from_json_string(x_json);
    let n: Option<Expr> = from_json_string(n_json);
    if let (Some(y), Some(x), Some(n)) = (y, x, n) {
        to_json_string(&special_functions::laguerre_differential_equation(&y, &x, &n))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_hermite_differential_equation(y_json: *const c_char, x_json: *const c_char, n_json: *const c_char) -> *mut c_char {
    let y: Option<Expr> = from_json_string(y_json);
    let x: Option<Expr> = from_json_string(x_json);
    let n: Option<Expr> = from_json_string(n_json);
    if let (Some(y), Some(x), Some(n)) = (y, x, n) {
        to_json_string(&special_functions::hermite_differential_equation(&y, &x, &n))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_legendre_rodrigues_formula(n_json: *const c_char, x_json: *const c_char) -> *mut c_char {
    let n: Option<Expr> = from_json_string(n_json);
    let x: Option<Expr> = from_json_string(x_json);
    if let (Some(n), Some(x)) = (n, x) {
        to_json_string(&special_functions::legendre_rodrigues_formula(&n, &x))
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_hermite_rodrigues_formula(n_json: *const c_char, x_json: *const c_char) -> *mut c_char {
    let n: Option<Expr> = from_json_string(n_json);
    let x: Option<Expr> = from_json_string(x_json);
    if let (Some(n), Some(x)) = (n, x) {
        to_json_string(&special_functions::hermite_rodrigues_formula(&n, &x))
    } else {
        std::ptr::null_mut()
    }
}
