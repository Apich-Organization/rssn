use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::multi_valued::*;

/// Computes general multi-valued logarithm (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_general_log(z_buf: BincodeBuffer, k_buf: BincodeBuffer) -> BincodeBuffer {
    let z: Option<Expr> = from_bincode_buffer(&z_buf);
    let k: Option<Expr> = from_bincode_buffer(&k_buf);

    if let (Some(z_expr), Some(k_expr)) = (z, k) {
        let result = general_log(&z_expr, &k_expr);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

/// Computes general multi-valued square root (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_general_sqrt(z_buf: BincodeBuffer, k_buf: BincodeBuffer) -> BincodeBuffer {
    let z: Option<Expr> = from_bincode_buffer(&z_buf);
    let k: Option<Expr> = from_bincode_buffer(&k_buf);

    if let (Some(z_expr), Some(k_expr)) = (z, k) {
        let result = general_sqrt(&z_expr, &k_expr);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

/// Computes general multi-valued power (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_general_power(
    z_buf: BincodeBuffer,
    w_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
) -> BincodeBuffer {
    let z: Option<Expr> = from_bincode_buffer(&z_buf);
    let w: Option<Expr> = from_bincode_buffer(&w_buf);
    let k: Option<Expr> = from_bincode_buffer(&k_buf);

    if let (Some(z_expr), Some(w_expr), Some(k_expr)) = (z, w, k) {
        let result = general_power(&z_expr, &w_expr, &k_expr);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

/// Computes general multi-valued n-th root (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_general_nth_root(
    z_buf: BincodeBuffer,
    n_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
) -> BincodeBuffer {
    let z: Option<Expr> = from_bincode_buffer(&z_buf);
    let n: Option<Expr> = from_bincode_buffer(&n_buf);
    let k: Option<Expr> = from_bincode_buffer(&k_buf);

    if let (Some(z_expr), Some(n_expr), Some(k_expr)) = (z, n, k) {
        let result = general_nth_root(&z_expr, &n_expr, &k_expr);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

/// Computes general multi-valued arcsin (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_general_arcsin(z_buf: BincodeBuffer, k_buf: BincodeBuffer) -> BincodeBuffer {
    let z: Option<Expr> = from_bincode_buffer(&z_buf);
    let k: Option<Expr> = from_bincode_buffer(&k_buf);

    if let (Some(z_expr), Some(k_expr)) = (z, k) {
        let result = general_arcsin(&z_expr, &k_expr);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

/// Computes general multi-valued arccos (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_general_arccos(
    z_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
    s_buf: BincodeBuffer,
) -> BincodeBuffer {
    let z: Option<Expr> = from_bincode_buffer(&z_buf);
    let k: Option<Expr> = from_bincode_buffer(&k_buf);
    let s: Option<Expr> = from_bincode_buffer(&s_buf);

    if let (Some(z_expr), Some(k_expr), Some(s_expr)) = (z, k, s) {
        let result = general_arccos(&z_expr, &k_expr, &s_expr);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

/// Computes general multi-valued arctan (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_general_arctan(z_buf: BincodeBuffer, k_buf: BincodeBuffer) -> BincodeBuffer {
    let z: Option<Expr> = from_bincode_buffer(&z_buf);
    let k: Option<Expr> = from_bincode_buffer(&k_buf);

    if let (Some(z_expr), Some(k_expr)) = (z, k) {
        let result = general_arctan(&z_expr, &k_expr);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

/// Computes argument (angle) of complex number (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_arg(z_buf: BincodeBuffer) -> BincodeBuffer {
    let z: Option<Expr> = from_bincode_buffer(&z_buf);

    if let Some(z_expr) = z {
        let result = arg(&z_expr);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

/// Computes absolute value (magnitude) of complex number (Bincode)
#[no_mangle]
pub extern "C" fn rssn_bincode_abs(z_buf: BincodeBuffer) -> BincodeBuffer {
    let z: Option<Expr> = from_bincode_buffer(&z_buf);

    if let Some(z_expr) = z {
        let result = abs(&z_expr);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}
