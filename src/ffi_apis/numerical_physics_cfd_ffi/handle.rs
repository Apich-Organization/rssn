//! Handle-based FFI API for numerical CFD functions.

use crate::numerical::physics_cfd;

// ============================================================================
// Fluid Properties
// ============================================================================

/// Returns kinematic viscosity of air at 20°C.
#[no_mangle]
pub extern "C" fn rssn_num_cfd_air_kinematic_viscosity() -> f64 {
    physics_cfd::FluidProperties::air().kinematic_viscosity()
}

/// Returns kinematic viscosity of water at 20°C.
#[no_mangle]
pub extern "C" fn rssn_num_cfd_water_kinematic_viscosity() -> f64 {
    physics_cfd::FluidProperties::water().kinematic_viscosity()
}

/// Returns Prandtl number of air.
#[no_mangle]
pub extern "C" fn rssn_num_cfd_air_prandtl_number() -> f64 {
    physics_cfd::FluidProperties::air().prandtl_number()
}

/// Returns Prandtl number of water.
#[no_mangle]
pub extern "C" fn rssn_num_cfd_water_prandtl_number() -> f64 {
    physics_cfd::FluidProperties::water().prandtl_number()
}

// ============================================================================
// Dimensionless Numbers
// ============================================================================

/// Calculates Reynolds number.
#[no_mangle]
pub extern "C" fn rssn_num_cfd_reynolds_number(
    velocity: f64,
    length: f64,
    kinematic_viscosity: f64,
) -> f64 {
    physics_cfd::reynolds_number(velocity, length, kinematic_viscosity)
}

/// Calculates Mach number.
#[no_mangle]
pub extern "C" fn rssn_num_cfd_mach_number(velocity: f64, speed_of_sound: f64) -> f64 {
    physics_cfd::mach_number(velocity, speed_of_sound)
}

/// Calculates Froude number.
#[no_mangle]
pub extern "C" fn rssn_num_cfd_froude_number(velocity: f64, length: f64, gravity: f64) -> f64 {
    physics_cfd::froude_number(velocity, length, gravity)
}

/// Calculates CFL number.
#[no_mangle]
pub extern "C" fn rssn_num_cfd_cfl_number(velocity: f64, dt: f64, dx: f64) -> f64 {
    physics_cfd::cfl_number(velocity, dt, dx)
}

/// Checks CFL stability.
#[no_mangle]
pub extern "C" fn rssn_num_cfd_check_cfl_stability(
    velocity: f64,
    dt: f64,
    dx: f64,
    max_cfl: f64,
) -> bool {
    physics_cfd::check_cfl_stability(velocity, dt, dx, max_cfl)
}

/// Calculates diffusion number.
#[no_mangle]
pub extern "C" fn rssn_num_cfd_diffusion_number(alpha: f64, dt: f64, dx: f64) -> f64 {
    physics_cfd::diffusion_number(alpha, dt, dx)
}
