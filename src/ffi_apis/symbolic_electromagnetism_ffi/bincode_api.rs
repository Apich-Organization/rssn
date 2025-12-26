//! Bincode-based FFI API for electromagnetism functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::electromagnetism;
use crate::symbolic::vector::Vector;

/// Calculates Lorentz force using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_lorentz_force(
    charge_buf : BincodeBuffer,
    e_field_buf : BincodeBuffer,
    velocity_buf : BincodeBuffer,
    b_field_buf : BincodeBuffer,
) -> BincodeBuffer {

    let charge : Option<Expr> = from_bincode_buffer(&charge_buf);

    let e_field : Option<Vector> = from_bincode_buffer(&e_field_buf);

    let velocity : Option<Vector> = from_bincode_buffer(&velocity_buf);

    let b_field : Option<Vector> = from_bincode_buffer(&b_field_buf);

    if let (Some(q), Some(e), Some(v), Some(b)) = (
        charge,
        e_field,
        velocity,
        b_field,
    ) {

        to_bincode_buffer(&electromagnetism::lorentz_force(&q, &e, &v, &b))
    } else {

        BincodeBuffer::empty()
    }
}

/// Calculates energy density using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_electromagnetic_energy_density(
    e_field_buf : BincodeBuffer,
    b_field_buf : BincodeBuffer,
) -> BincodeBuffer {

    let e_field : Option<Vector> = from_bincode_buffer(&e_field_buf);

    let b_field : Option<Vector> = from_bincode_buffer(&b_field_buf);

    if let (Some(e), Some(b)) = (e_field, b_field) {

        to_bincode_buffer(&electromagnetism::energy_density(&e, &b))
    } else {

        BincodeBuffer::empty()
    }
}
