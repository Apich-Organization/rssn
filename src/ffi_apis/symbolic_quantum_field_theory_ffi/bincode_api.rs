use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::symbolic::core::Expr;
use crate::symbolic::quantum_field_theory;

/// Computes a propagator using Bincode.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_qft_propagator(
    p_buf: BincodeBuffer,
    m_buf: BincodeBuffer,
    is_fermion: bool,
) -> BincodeBuffer {

    let p: Option<Expr> =
        from_bincode_buffer(&p_buf);

    let m: Option<Expr> =
        from_bincode_buffer(&m_buf);

    match (p, m) { (Some(p), Some(m)) => {

        to_bincode_buffer(&quantum_field_theory::propagator(&p, &m, is_fermion))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Lagrangian density for a scalar field using Bincode.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_scalar_field_lagrangian(
    phi_buf: BincodeBuffer,
    m_buf: BincodeBuffer,
) -> BincodeBuffer {

    let phi: Option<Expr> =
        from_bincode_buffer(&phi_buf);

    let m: Option<Expr> =
        from_bincode_buffer(&m_buf);

    match (phi, m)
    { (Some(phi), Some(m)) => {

        to_bincode_buffer(&quantum_field_theory::scalar_field_lagrangian(&phi, &m))
    } _ => {

        BincodeBuffer::empty()
    }}
}
