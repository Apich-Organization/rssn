//! JSON-based FFI API for electromagnetism functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::electromagnetism;
use crate::symbolic::vector::Vector;

/// Calculates Lorentz force using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_lorentz_force(
    charge_json: *const c_char,
    e_field_json: *const c_char,
    velocity_json: *const c_char,
    b_field_json: *const c_char,
) -> *mut c_char {

    let charge: Option<Expr> =
        from_json_string(charge_json);

    let e_field: Option<Vector> =
        from_json_string(e_field_json);

    let velocity: Option<Vector> =
        from_json_string(velocity_json);

    let b_field: Option<Vector> =
        from_json_string(b_field_json);

    if let (
        Some(q),
        Some(e),
        Some(v),
        Some(b),
    ) = (
        charge,
        e_field,
        velocity,
        b_field,
    ) {

        to_json_string(&electromagnetism::lorentz_force(&q, &e, &v, &b))
    } else {

        std::ptr::null_mut()
    }
}

/// Calculates energy density using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_electromagnetic_energy_density(
    e_field_json: *const c_char,
    b_field_json: *const c_char,
) -> *mut c_char {

    let e_field: Option<Vector> =
        from_json_string(e_field_json);

    let b_field: Option<Vector> =
        from_json_string(b_field_json);

    if let (Some(e), Some(b)) =
        (e_field, b_field)
    {

        to_json_string(&electromagnetism::energy_density(&e, &b))
    } else {

        std::ptr::null_mut()
    }
}
