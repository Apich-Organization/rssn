use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::symbolic::core::Expr;
use crate::symbolic::multi_valued::abs;
use crate::symbolic::multi_valued::arg;
use crate::symbolic::multi_valued::general_arccos;
use crate::symbolic::multi_valued::general_arcsin;
use crate::symbolic::multi_valued::general_arctan;
use crate::symbolic::multi_valued::general_log;
use crate::symbolic::multi_valued::general_nth_root;
use crate::symbolic::multi_valued::general_power;
use crate::symbolic::multi_valued::general_sqrt;

/// Computes general multi-valued logarithm (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_general_log(
    z_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
) -> BincodeBuffer {

    let z: Option<Expr> =
        from_bincode_buffer(&z_buf);

    let k: Option<Expr> =
        from_bincode_buffer(&k_buf);

    match (z, k)
    { (
        Some(z_expr),
        Some(k_expr),
    ) => {

        let result = general_log(
            &z_expr,
            &k_expr,
        );

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Computes general multi-valued square root (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_general_sqrt(
    z_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
) -> BincodeBuffer {

    let z: Option<Expr> =
        from_bincode_buffer(&z_buf);

    let k: Option<Expr> =
        from_bincode_buffer(&k_buf);

    match (z, k)
    { (
        Some(z_expr),
        Some(k_expr),
    ) => {

        let result = general_sqrt(
            &z_expr,
            &k_expr,
        );

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Computes general multi-valued power (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_general_power(
    z_buf: BincodeBuffer,
    w_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
) -> BincodeBuffer {

    let z: Option<Expr> =
        from_bincode_buffer(&z_buf);

    let w: Option<Expr> =
        from_bincode_buffer(&w_buf);

    let k: Option<Expr> =
        from_bincode_buffer(&k_buf);

    match (z, w, k)
    { (
        Some(z_expr),
        Some(w_expr),
        Some(k_expr),
    ) => {

        let result = general_power(
            &z_expr,
            &w_expr,
            &k_expr,
        );

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Computes general multi-valued n-th root (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_general_nth_root(
    z_buf: BincodeBuffer,
    n_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
) -> BincodeBuffer {

    let z: Option<Expr> =
        from_bincode_buffer(&z_buf);

    let n: Option<Expr> =
        from_bincode_buffer(&n_buf);

    let k: Option<Expr> =
        from_bincode_buffer(&k_buf);

    match (z, n, k)
    { (
        Some(z_expr),
        Some(n_expr),
        Some(k_expr),
    ) => {

        let result = general_nth_root(
            &z_expr,
            &n_expr,
            &k_expr,
        );

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Computes general multi-valued arcsin (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_general_arcsin(
    z_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
) -> BincodeBuffer {

    let z: Option<Expr> =
        from_bincode_buffer(&z_buf);

    let k: Option<Expr> =
        from_bincode_buffer(&k_buf);

    match (z, k)
    { (
        Some(z_expr),
        Some(k_expr),
    ) => {

        let result = general_arcsin(
            &z_expr,
            &k_expr,
        );

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Computes general multi-valued arccos (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_general_arccos(
    z_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
    s_buf: BincodeBuffer,
) -> BincodeBuffer {

    let z: Option<Expr> =
        from_bincode_buffer(&z_buf);

    let k: Option<Expr> =
        from_bincode_buffer(&k_buf);

    let s: Option<Expr> =
        from_bincode_buffer(&s_buf);

    match (z, k, s)
    { (
        Some(z_expr),
        Some(k_expr),
        Some(s_expr),
    ) => {

        let result = general_arccos(
            &z_expr,
            &k_expr,
            &s_expr,
        );

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Computes general multi-valued arctan (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_general_arctan(
    z_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
) -> BincodeBuffer {

    let z: Option<Expr> =
        from_bincode_buffer(&z_buf);

    let k: Option<Expr> =
        from_bincode_buffer(&k_buf);

    match (z, k)
    { (
        Some(z_expr),
        Some(k_expr),
    ) => {

        let result = general_arctan(
            &z_expr,
            &k_expr,
        );

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Computes argument (angle) of complex number (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_arg(
    z_buf: BincodeBuffer
) -> BincodeBuffer {

    let z: Option<Expr> =
        from_bincode_buffer(&z_buf);

    if let Some(z_expr) = z {

        let result = arg(&z_expr);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes absolute value (magnitude) of complex number (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_abs(
    z_buf: BincodeBuffer
) -> BincodeBuffer {

    let z: Option<Expr> =
        from_bincode_buffer(&z_buf);

    if let Some(z_expr) = z {

        let result = abs(&z_expr);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}
