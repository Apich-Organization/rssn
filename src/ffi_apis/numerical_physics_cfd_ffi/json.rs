//! JSON-based FFI API for numerical CFD functions.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::physics_cfd;

// ============================================================================
// Input/Output structs
// ============================================================================

#[derive(Deserialize)]

struct FluidPropertiesInput {
    density: f64,
    dynamic_viscosity: f64,
    thermal_conductivity: f64,
    specific_heat: f64,
}

#[derive(Serialize)]

struct FluidPropertiesOutput {
    kinematic_viscosity: f64,
    thermal_diffusivity: f64,
    prandtl_number: f64,
}

#[derive(Deserialize)]

struct ReynoldsInput {
    velocity: f64,
    length: f64,
    kinematic_viscosity: f64,
}

#[derive(Deserialize)]

struct CflInput {
    velocity: f64,
    dt: f64,
    dx: f64,
}

#[derive(Deserialize)]

struct Advection1DInput {
    u0: Vec<f64>,
    c: f64,
    dx: f64,
    dt: f64,
    num_steps: usize,
}

#[derive(Deserialize)]

struct Diffusion1DInput {
    u0: Vec<f64>,
    alpha: f64,
    dx: f64,
    dt: f64,
    num_steps: usize,
}

#[derive(Deserialize)]

struct AdvectionDiffusion1DInput {
    u0: Vec<f64>,
    c: f64,
    alpha: f64,
    dx: f64,
    dt: f64,
    num_steps: usize,
}

#[derive(Deserialize)]

struct Burgers1DInput {
    u0: Vec<f64>,
    nu: f64,
    dx: f64,
    dt: f64,
    num_steps: usize,
}

// ============================================================================
// Fluid Properties
// ============================================================================

/// Computes derived fluid properties from fundamental properties using JSON serialization.
///
/// This function calculates kinematic viscosity (ν = μ/ρ), thermal diffusivity (α = `k/(ρc_p)`),
/// and Prandtl number (Pr = ν/α) from input material properties.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `density`: Fluid density ρ (kg/m³)
///   - `dynamic_viscosity`: Dynamic viscosity μ (Pa·s)
///   - `thermal_conductivity`: Thermal conductivity k (W/(m·K))
///   - `specific_heat`: Specific heat capacity `c_p` (J/(kg·K))
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult` with:
/// - `kinematic_viscosity`: ν (m²/s)
/// - `thermal_diffusivity`: α (m²/s)
/// - `prandtl_number`: Pr (dimensionless)
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_cfd_fluid_properties_json(
    input: *const c_char
) -> *mut c_char {

    let input : FluidPropertiesInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<FluidPropertiesOutput, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let fluid = physics_cfd::FluidProperties::new(
        input.density,
        input.dynamic_viscosity,
        input.thermal_conductivity,
        input.specific_heat,
    );

    let output =
        FluidPropertiesOutput {
            kinematic_viscosity: fluid
                .kinematic_viscosity(),
            thermal_diffusivity: fluid
                .thermal_diffusivity(),
            prandtl_number: fluid
                .prandtl_number(),
        };

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(output),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Returns standard air properties at sea level and 15°C using JSON serialization.
///
/// Provides reference fluid properties for air (ρ ≈ 1.225 kg/m³, μ ≈ 1.81×10⁻⁵ Pa·s).
///
/// # Arguments
///
/// * `_input` - Unused parameter for API consistency (can be null)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult` with standard air properties:
/// - `kinematic_viscosity`: ν (m²/s)
/// - `thermal_diffusivity`: α (m²/s)
/// - `prandtl_number`: Pr (dimensionless, ≈ 0.71 for air)
///
/// # Safety
///
/// This function is unsafe because it returns a raw pointer that the caller must free.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_cfd_air_properties_json(
    _input: *const c_char
) -> *mut c_char {

    let fluid = physics_cfd::FluidProperties::air();

    let output =
        FluidPropertiesOutput {
            kinematic_viscosity: fluid
                .kinematic_viscosity(),
            thermal_diffusivity: fluid
                .thermal_diffusivity(),
            prandtl_number: fluid
                .prandtl_number(),
        };

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(output),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Returns standard water properties at 20°C using JSON serialization.
///
/// Provides reference fluid properties for water (ρ ≈ 998 kg/m³, μ ≈ 1.0×10⁻³ Pa·s).
///
/// # Arguments
///
/// * `_input` - Unused parameter for API consistency (can be null)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult` with standard water properties:
/// - `kinematic_viscosity`: ν (m²/s)
/// - `thermal_diffusivity`: α (m²/s)
/// - `prandtl_number`: Pr (dimensionless, ≈ 7 for water)
///
/// # Safety
///
/// This function is unsafe because it returns a raw pointer that the caller must free.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_cfd_water_properties_json(
    _input: *const c_char
) -> *mut c_char {

    let fluid = physics_cfd::FluidProperties::water();

    let output =
        FluidPropertiesOutput {
            kinematic_viscosity: fluid
                .kinematic_viscosity(),
            thermal_diffusivity: fluid
                .thermal_diffusivity(),
            prandtl_number: fluid
                .prandtl_number(),
        };

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(output),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

// ============================================================================
// Dimensionless Numbers
// ============================================================================

/// Computes the Reynolds number for fluid flow using JSON serialization.
///
/// The Reynolds number characterizes the flow regime (laminar vs. turbulent):
/// Re = (velocity × length) / `kinematic_viscosity`
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `velocity`: Flow velocity (m/s)
///   - `length`: Characteristic length scale (m)
///   - `kinematic_viscosity`: Kinematic viscosity ν (m²/s)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the computed Reynolds number (dimensionless).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_cfd_reynolds_number_json(
    input: *const c_char
) -> *mut c_char {

    let input : ReynoldsInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let re =
        physics_cfd::reynolds_number(
            input.velocity,
            input.length,
            input.kinematic_viscosity,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(re),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the Courant-Friedrichs-Lewy (CFL) number using JSON serialization.
///
/// The CFL condition ensures numerical stability in explicit time-stepping schemes:
/// CFL = (velocity × dt) / dx. For stability, CFL ≤ 1 is typically required.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `velocity`: Flow velocity or wave speed (m/s)
///   - `dt`: Time step size (s)
///   - `dx`: Spatial grid spacing (m)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the computed CFL number (dimensionless).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_cfd_cfl_number_json(
    input: *const c_char
) -> *mut c_char {

    let input : CflInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let cfl = physics_cfd::cfl_number(
        input.velocity,
        input.dt,
        input.dx,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(cfl),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

// ============================================================================
// 1D Solvers
// ============================================================================

/// Solves the 1D advection equation using JSON serialization.
///
/// The advection equation describes transport by a velocity field: ∂u/∂t + c ∂u/∂x = 0.
/// Uses an explicit finite difference scheme to time-step the solution.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `u0`: Initial condition array (spatial distribution)
///   - `c`: Advection velocity (constant, m/s)
///   - `dx`: Spatial grid spacing (m)
///   - `dt`: Time step size (s)
///   - `num_steps`: Number of time steps to compute
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<Vec<f64>>, String>` with
/// the solution history, where each inner vector represents the spatial distribution at one time step.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_cfd_solve_advection_1d_json(
    input: *const c_char
) -> *mut c_char {

    let input : Advection1DInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<Vec<f64>>, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let results =
        physics_cfd::solve_advection_1d(
            &input.u0,
            input.c,
            input.dx,
            input.dt,
            input.num_steps,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(results),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Solves the 1D diffusion equation using JSON serialization.
///
/// The diffusion equation models heat conduction or mass diffusion: ∂u/∂t = α ∂²u/∂x².
/// Uses an explicit finite difference scheme for time integration.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `u0`: Initial condition array (spatial distribution)
///   - `alpha`: Diffusivity coefficient α (m²/s)
///   - `dx`: Spatial grid spacing (m)
///   - `dt`: Time step size (s)
///   - `num_steps`: Number of time steps to compute
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<Vec<f64>>, String>` with
/// the solution history, where each inner vector represents the spatial distribution at one time step.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_cfd_solve_diffusion_1d_json(
    input: *const c_char
) -> *mut c_char {

    let input : Diffusion1DInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<Vec<f64>>, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let results =
        physics_cfd::solve_diffusion_1d(
            &input.u0,
            input.alpha,
            input.dx,
            input.dt,
            input.num_steps,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(results),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Solves the 1D advection-diffusion equation using JSON serialization.
///
/// The advection-diffusion equation combines transport and diffusion:
/// ∂u/∂t + c ∂u/∂x = α ∂²u/∂x².
/// This models phenomena like pollutant dispersion in a moving fluid.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `u0`: Initial condition array (spatial distribution)
///   - `c`: Advection velocity (m/s)
///   - `alpha`: Diffusivity coefficient α (m²/s)
///   - `dx`: Spatial grid spacing (m)
///   - `dt`: Time step size (s)
///   - `num_steps`: Number of time steps to compute
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<Vec<f64>>, String>` with
/// the solution history, where each inner vector represents the spatial distribution at one time step.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_cfd_solve_advection_diffusion_1d_json(
    input: *const c_char
) -> *mut c_char {

    let input : AdvectionDiffusion1DInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<Vec<f64>>, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let results = physics_cfd::solve_advection_diffusion_1d(
        &input.u0,
        input.c,
        input.alpha,
        input.dx,
        input.dt,
        input.num_steps,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(results),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Solves the 1D Burgers equation using JSON serialization.
///
/// The Burgers equation is a nonlinear PDE combining convection and diffusion:
/// ∂u/∂t + u ∂u/∂x = ν ∂²u/∂x².
/// It models shock wave formation and is often used as a simplified model for fluid turbulence.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `u0`: Initial condition array (spatial distribution)
///   - `nu`: Kinematic viscosity ν (m²/s)
///   - `dx`: Spatial grid spacing (m)
///   - `dt`: Time step size (s)
///   - `num_steps`: Number of time steps to compute
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<Vec<f64>>, String>` with
/// the solution history, where each inner vector represents the spatial distribution at one time step.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_cfd_solve_burgers_1d_json(
    input: *const c_char
) -> *mut c_char {

    let input : Burgers1DInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<Vec<f64>>, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let results =
        physics_cfd::solve_burgers_1d(
            &input.u0,
            input.nu,
            input.dx,
            input.dt,
            input.num_steps,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(results),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}
