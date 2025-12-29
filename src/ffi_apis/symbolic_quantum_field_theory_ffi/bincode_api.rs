use crate::ffi_apis::common::{BincodeBuffer, from_bincode_buffer, to_bincode_buffer};
use crate::symbolic::core::Expr;
use crate::symbolic::quantum_field_theory;

/// Computes a propagator using Bincode.
#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_qft_propagator(
    p_buf: BincodeBuffer,
    m_buf: BincodeBuffer,
    is_fermion: bool,
) -> BincodeBuffer {

    let p: Option<Expr> =
        from_bincode_buffer(&p_buf);

    let m: Option<Expr> =
        from_bincode_buffer(&m_buf);

    if let (Some(p), Some(m)) = (p, m) {

        to_bincode_buffer(&quantum_field_theory::propagator(&p, &m, is_fermion))
    } else {

        BincodeBuffer::empty()
    }
}

/// Lagrangian density for a scalar field using Bincode.
#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_scalar_field_lagrangian(
    phi_buf: BincodeBuffer,
    m_buf: BincodeBuffer,
) -> BincodeBuffer {

    let phi: Option<Expr> =
        from_bincode_buffer(&phi_buf);

    let m: Option<Expr> =
        from_bincode_buffer(&m_buf);

    if let (Some(phi), Some(m)) =
        (phi, m)
    {

        to_bincode_buffer(&quantum_field_theory::scalar_field_lagrangian(&phi, &m))
    } else {

        BincodeBuffer::empty()
    }
}
