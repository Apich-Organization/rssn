//! Handle-based FFI API for numerical MD functions.

use crate::numerical::physics_md;

// ============================================================================
// Physical Constants
// ============================================================================

/// Returns Boltzmann constant in SI units.
#[unsafe(no_mangle)]

pub const extern "C" fn rssn_num_md_boltzmann_constant_si(
) -> f64 {

    physics_md::BOLTZMANN_CONSTANT_SI
}

/// Returns Avogadro's number.
#[unsafe(no_mangle)]

pub const extern "C" fn rssn_num_md_avogadro_number(
) -> f64 {

    physics_md::AVOGADRO_NUMBER
}

/// Returns temperature unit for argon in reduced units.
#[unsafe(no_mangle)]

pub const extern "C" fn rssn_num_md_temperature_unit_argon(
) -> f64 {

    physics_md::TEMPERATURE_UNIT_ARGON
}

// ============================================================================
// System Properties
// ============================================================================

/// Calculates CFL number for MD stability.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_md_cfl_check(
    velocity: f64,
    dt: f64,
    sigma: f64,
) -> f64 {

    velocity * dt / sigma
}

// ============================================================================
// Periodic Boundary Conditions
// ============================================================================

/// Applies minimum image convention to a 1D distance.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_md_minimum_image_1d(
    dx: f64,
    box_length: f64,
) -> f64 {

    let mut d = dx % box_length;

    if d > box_length / 2.0 {

        d -= box_length;
    } else if d < -box_length / 2.0 {

        d += box_length;
    }

    d
}

/// Applies periodic boundary condition to a 1D position.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_num_md_apply_pbc_1d(
    x: f64,
    box_length: f64,
) -> f64 {

    let mut wrapped = x % box_length;

    if wrapped < 0.0 {

        wrapped += box_length;
    }

    wrapped
}
