//! Bincode-based FFI API for solid-state physics functions.

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::solid_state_physics;

/// Computes the density of states for a 3D electron gas using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_density_of_states_3d(
    energy_buf : BincodeBuffer,
    mass_buf : BincodeBuffer,
    volume_buf : BincodeBuffer,
) -> BincodeBuffer {

    let energy : Option<Expr> = from_bincode_buffer(&energy_buf);

    let mass : Option<Expr> = from_bincode_buffer(&mass_buf);

    let volume : Option<Expr> = from_bincode_buffer(&volume_buf);

    if let (Some(energy), Some(mass), Some(volume)) = (energy, mass, volume) {

        to_bincode_buffer(
            &solid_state_physics::density_of_states_3d(
                &energy,
                &mass,
                &volume,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes Fermi energy for a 3D electron gas using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_fermi_energy_3d(
    concentration_buf : BincodeBuffer,
    mass_buf : BincodeBuffer,
) -> BincodeBuffer {

    let concentration : Option<Expr> = from_bincode_buffer(&concentration_buf);

    let mass : Option<Expr> = from_bincode_buffer(&mass_buf);

    if let (Some(concentration), Some(mass)) = (concentration, mass) {

        to_bincode_buffer(
            &solid_state_physics::fermi_energy_3d(
                &concentration,
                &mass,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes Drude conductivity using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_drude_conductivity(
    n_buf : BincodeBuffer,
    e_buf : BincodeBuffer,
    tau_buf : BincodeBuffer,
    mass_buf : BincodeBuffer,
) -> BincodeBuffer {

    let n : Option<Expr> = from_bincode_buffer(&n_buf);

    let e : Option<Expr> = from_bincode_buffer(&e_buf);

    let tau : Option<Expr> = from_bincode_buffer(&tau_buf);

    let mass : Option<Expr> = from_bincode_buffer(&mass_buf);

    if let (Some(n), Some(e), Some(tau), Some(mass)) = (n, e, tau, mass) {

        to_bincode_buffer(&solid_state_physics::drude_conductivity(&n, &e, &tau, &mass))
    } else {

        BincodeBuffer::empty()
    }
}
