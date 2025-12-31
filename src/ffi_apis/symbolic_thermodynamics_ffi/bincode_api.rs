//! Bincode-based FFI API for thermodynamics functions.

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::symbolic::core::Expr;
use crate::symbolic::thermodynamics;

/// Calculates ideal gas Law using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_ideal_gas_law(
    p_buf: BincodeBuffer,
    v_buf: BincodeBuffer,
    n_buf: BincodeBuffer,
    r_buf: BincodeBuffer,
    t_buf: BincodeBuffer,
) -> BincodeBuffer {

    let p: Option<Expr> =
        from_bincode_buffer(&p_buf);

    let v: Option<Expr> =
        from_bincode_buffer(&v_buf);

    let n: Option<Expr> =
        from_bincode_buffer(&n_buf);

    let r: Option<Expr> =
        from_bincode_buffer(&r_buf);

    let t: Option<Expr> =
        from_bincode_buffer(&t_buf);

    match (p, v, n, r, t)
    { (
        Some(p),
        Some(v),
        Some(n),
        Some(r),
        Some(t),
    ) => {

        to_bincode_buffer(&thermodynamics::ideal_gas_law(&p, &v, &n, &r, &t))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Calculates Gibbs Free Energy using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_gibbs_free_energy(
    h_buf: BincodeBuffer,
    t_buf: BincodeBuffer,
    s_buf: BincodeBuffer,
) -> BincodeBuffer {

    let h: Option<Expr> =
        from_bincode_buffer(&h_buf);

    let t: Option<Expr> =
        from_bincode_buffer(&t_buf);

    let s: Option<Expr> =
        from_bincode_buffer(&s_buf);

    match (h, t, s)
    { (Some(h), Some(t), Some(s)) => {

        to_bincode_buffer(&thermodynamics::gibbs_free_energy(&h, &t, &s))
    } _ => {

        BincodeBuffer::empty()
    }}
}
