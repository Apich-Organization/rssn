//! Bincode-based FFI API for solid-state physics functions.

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::symbolic::core::Expr;
use crate::symbolic::solid_state_physics;

/// Computes the density of states for a 3D electron gas using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_density_of_states_3d(
    energy_buf: BincodeBuffer,
    mass_buf: BincodeBuffer,
    volume_buf: BincodeBuffer,
) -> BincodeBuffer {

    let energy: Option<Expr> =
        from_bincode_buffer(
            &energy_buf,
        );

    let mass: Option<Expr> =
        from_bincode_buffer(&mass_buf);

    let volume: Option<Expr> =
        from_bincode_buffer(
            &volume_buf,
        );

    match (energy, mass, volume)
    { (
        Some(energy),
        Some(mass),
        Some(volume),
    ) => {

        to_bincode_buffer(
            &solid_state_physics::density_of_states_3d(
                &energy,
                &mass,
                &volume,
            ),
        )
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Computes Fermi energy for a 3D electron gas using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_fermi_energy_3d(
    concentration_buf: BincodeBuffer,
    mass_buf: BincodeBuffer,
) -> BincodeBuffer {

    let concentration: Option<Expr> =
        from_bincode_buffer(
            &concentration_buf,
        );

    let mass: Option<Expr> =
        from_bincode_buffer(&mass_buf);

    match (concentration, mass)
    { (
        Some(concentration),
        Some(mass),
    ) => {

        to_bincode_buffer(
            &solid_state_physics::fermi_energy_3d(
                &concentration,
                &mass,
            ),
        )
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Computes Drude conductivity using Bincode.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_drude_conductivity(
    n_buf: BincodeBuffer,
    e_buf: BincodeBuffer,
    tau_buf: BincodeBuffer,
    mass_buf: BincodeBuffer,
) -> BincodeBuffer {

    let n: Option<Expr> =
        from_bincode_buffer(&n_buf);

    let e: Option<Expr> =
        from_bincode_buffer(&e_buf);

    let tau: Option<Expr> =
        from_bincode_buffer(&tau_buf);

    let mass: Option<Expr> =
        from_bincode_buffer(&mass_buf);

    match (n, e, tau, mass)
    { (
        Some(n),
        Some(e),
        Some(tau),
        Some(mass),
    ) => {

        to_bincode_buffer(&solid_state_physics::drude_conductivity(&n, &e, &tau, &mass))
    } _ => {

        BincodeBuffer::empty()
    }}
}
