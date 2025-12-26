//! Handle-based FFI API for numerical multi-valued functions.

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::multi_valued;
use crate::symbolic::core::Expr;
use num_complex::Complex;

/// Finds a root of a complex function using Newton's method.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_mv_newton_method_complex(
    f_ptr: *const Expr,
    f_prime_ptr: *const Expr,
    start_re: f64,
    start_im: f64,
    tolerance: f64,
    max_iter: usize,
    res_re: *mut f64,
    res_im: *mut f64,
) -> i32 {

    if f_ptr.is_null()
        || f_prime_ptr.is_null()
        || res_re.is_null()
        || res_im.is_null()
    {

        update_last_error(
            "Null pointer passed to \
             rssn_num_mv_newton_method_complex"
                .to_string(),
        );

        return -1;
    }

    let f = &*f_ptr;

    let f_prime = &*f_prime_ptr;

    let start_point = Complex::new(
        start_re, start_im,
    );

    match multi_valued::newton_method_complex(
        f,
        f_prime,
        start_point,
        tolerance,
        max_iter,
    ) {
        | Some(root) => {

            *res_re = root.re;

            *res_im = root.im;

            0
        },
        | None => {

            update_last_error(
                "Newton's method failed to converge"
                    .to_string(),
            );

            -1
        },
    }
}

/// Computes the k-th branch of the complex logarithm.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_mv_complex_log_k(
    re: f64,
    im: f64,
    k: i32,
    res_re: *mut f64,
    res_im: *mut f64,
) {

    let z = Complex::new(re, im);

    let res =
        multi_valued::complex_log_k(
            z, k,
        );

    *res_re = res.re;

    *res_im = res.im;
}

/// Computes the k-th branch of the complex square root.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_mv_complex_sqrt_k(
    re: f64,
    im: f64,
    k: i32,
    res_re: *mut f64,
    res_im: *mut f64,
) {

    let z = Complex::new(re, im);

    let res =
        multi_valued::complex_sqrt_k(
            z, k,
        );

    *res_re = res.re;

    *res_im = res.im;
}

/// Computes the k-th branch of the complex power z^w.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_mv_complex_pow_k(
    z_re: f64,
    z_im: f64,
    w_re: f64,
    w_im: f64,
    k: i32,
    res_re: *mut f64,
    res_im: *mut f64,
) {

    let z = Complex::new(z_re, z_im);

    let w = Complex::new(w_re, w_im);

    let res =
        multi_valued::complex_pow_k(
            z, w, k,
        );

    *res_re = res.re;

    *res_im = res.im;
}

/// Computes the k-th branch of the complex n-th root.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_mv_complex_nth_root_k(
    re: f64,
    im: f64,
    n: u32,
    k: i32,
    res_re: *mut f64,
    res_im: *mut f64,
) {

    let z = Complex::new(re, im);

    let res = multi_valued::complex_nth_root_k(z, n, k);

    *res_re = res.re;

    *res_im = res.im;
}

/// Computes the k-th branch of the complex arcsine.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_mv_complex_arcsin_k(
    re: f64,
    im: f64,
    k: i32,
    res_re: *mut f64,
    res_im: *mut f64,
) {

    let z = Complex::new(re, im);

    let res =
        multi_valued::complex_arcsin_k(
            z, k,
        );

    *res_re = res.re;

    *res_im = res.im;
}

/// Computes the k-th branch of the complex arccosine.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_mv_complex_arccos_k(
    re: f64,
    im: f64,
    k: i32,
    s: i32,
    res_re: *mut f64,
    res_im: *mut f64,
) {

    let z = Complex::new(re, im);

    let res =
        multi_valued::complex_arccos_k(
            z, k, s,
        );

    *res_re = res.re;

    *res_im = res.im;
}

/// Computes the k-th branch of the complex arctangent.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_mv_complex_arctan_k(
    re: f64,
    im: f64,
    k: i32,
    res_re: *mut f64,
    res_im: *mut f64,
) {

    let z = Complex::new(re, im);

    let res =
        multi_valued::complex_arctan_k(
            z, k,
        );

    *res_re = res.re;

    *res_im = res.im;
}
