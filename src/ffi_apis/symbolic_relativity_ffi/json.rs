//! JSON-based FFI API for relativity functions.

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::relativity;
use std::os::raw::c_char;

/// Calculates Lorentz factor using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_lorentz_factor(velocity_json: *const c_char) -> *mut c_char {

    let velocity: Option<Expr> = from_json_string(velocity_json);

    if let Some(v) = velocity {

        to_json_string(&relativity::lorentz_factor(&v))
    } else {

        std::ptr::null_mut()
    }
}

/// Calculates mass-energy equivalence using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_mass_energy_equivalence(mass_json: *const c_char) -> *mut c_char {

    let mass: Option<Expr> = from_json_string(mass_json);

    if let Some(m) = mass {

        to_json_string(&relativity::mass_energy_equivalence(&m))
    } else {

        std::ptr::null_mut()
    }
}

/// Calculates Schwarzschild radius using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_schwarzschild_radius(mass_json: *const c_char) -> *mut c_char {

    let mass: Option<Expr> = from_json_string(mass_json);

    if let Some(m) = mass {

        to_json_string(&relativity::schwarzschild_radius(&m))
    } else {

        std::ptr::null_mut()
    }
}
