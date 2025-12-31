//! Bincode-based FFI API for electromagnetism functions.

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::symbolic::core::Expr;
use crate::symbolic::electromagnetism;
use crate::symbolic::vector::Vector;

/// Calculates Lorentz force using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_lorentz_force(
    charge_buf: BincodeBuffer,
    e_field_buf: BincodeBuffer,
    velocity_buf: BincodeBuffer,
    b_field_buf: BincodeBuffer,
) -> BincodeBuffer {

    let charge: Option<Expr> =
        from_bincode_buffer(
            &charge_buf,
        );

    let e_field: Option<Vector> =
        from_bincode_buffer(
            &e_field_buf,
        );

    let velocity: Option<Vector> =
        from_bincode_buffer(
            &velocity_buf,
        );

    let b_field: Option<Vector> =
        from_bincode_buffer(
            &b_field_buf,
        );

    match (
        charge,
        e_field,
        velocity,
        b_field,
    ) { (
        Some(q),
        Some(e),
        Some(v),
        Some(b),
    ) => {

        to_bincode_buffer(&electromagnetism::lorentz_force(&q, &e, &v, &b))
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Calculates energy density using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_electromagnetic_energy_density(
    e_field_buf: BincodeBuffer,
    b_field_buf: BincodeBuffer,
) -> BincodeBuffer {

    let e_field: Option<Vector> =
        from_bincode_buffer(
            &e_field_buf,
        );

    let b_field: Option<Vector> =
        from_bincode_buffer(
            &b_field_buf,
        );

    match (e_field, b_field)
    { (Some(e), Some(b)) => {

        to_bincode_buffer(&electromagnetism::energy_density(&e, &b))
    } _ => {

        BincodeBuffer::empty()
    }}
}
