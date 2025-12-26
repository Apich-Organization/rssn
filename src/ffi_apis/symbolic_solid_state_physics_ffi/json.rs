//! JSON-based FFI API for solid-state physics functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::solid_state_physics;

/// Computes the density of states for a 3D electron gas using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_density_of_states_3d(
    energy_json : *const c_char,
    mass_json : *const c_char,
    volume_json : *const c_char,
) -> *mut c_char {

    let energy : Option<Expr> =
        from_json_string(energy_json);

    let mass : Option<Expr> =
        from_json_string(mass_json);

    let volume : Option<Expr> =
        from_json_string(volume_json);

    if let (
        Some(energy),
        Some(mass),
        Some(volume),
    ) = (energy, mass, volume)
    {

        to_json_string(
            &solid_state_physics::density_of_states_3d(
                &energy, &mass, &volume,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Computes Fermi energy for a 3D electron gas using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_fermi_energy_3d(
    concentration_json : *const c_char,
    mass_json : *const c_char,
) -> *mut c_char {

    let concentration : Option<Expr> =
        from_json_string(
            concentration_json,
        );

    let mass : Option<Expr> =
        from_json_string(mass_json);

    if let (
        Some(concentration),
        Some(mass),
    ) = (concentration, mass)
    {

        to_json_string(
            &solid_state_physics::fermi_energy_3d(
                &concentration,
                &mass,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Computes Drude conductivity using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_drude_conductivity(
    n_json : *const c_char,
    e_json : *const c_char,
    tau_json : *const c_char,
    mass_json : *const c_char,
) -> *mut c_char {

    let n : Option<Expr> =
        from_json_string(n_json);

    let e : Option<Expr> =
        from_json_string(e_json);

    let tau : Option<Expr> =
        from_json_string(tau_json);

    let mass : Option<Expr> =
        from_json_string(mass_json);

    if let (
        Some(n),
        Some(e),
        Some(tau),
        Some(mass),
    ) = (n, e, tau, mass)
    {

        to_json_string(
            &solid_state_physics::drude_conductivity(
                &n, &e, &tau, &mass,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}
