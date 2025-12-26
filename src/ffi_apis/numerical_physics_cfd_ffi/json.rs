//! JSON-based FFI API for numerical CFD functions.

use crate::ffi_apis::common::{from_json_string, to_c_string};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::physics_cfd;
use serde::{Deserialize, Serialize};
use std::os::raw::c_char;

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

#[no_mangle]

pub unsafe extern "C" fn rssn_num_cfd_fluid_properties_json(input: *const c_char) -> *mut c_char {

    let input: FluidPropertiesInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<FluidPropertiesOutput, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        }
    };

    let fluid = physics_cfd::FluidProperties::new(
        input.density,
        input.dynamic_viscosity,
        input.thermal_conductivity,
        input.specific_heat,
    );

    let output = FluidPropertiesOutput {
        kinematic_viscosity: fluid.kinematic_viscosity(),
        thermal_diffusivity: fluid.thermal_diffusivity(),
        prandtl_number: fluid.prandtl_number(),
    };

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(output),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_cfd_air_properties_json(_input: *const c_char) -> *mut c_char {

    let fluid = physics_cfd::FluidProperties::air();

    let output = FluidPropertiesOutput {
        kinematic_viscosity: fluid.kinematic_viscosity(),
        thermal_diffusivity: fluid.thermal_diffusivity(),
        prandtl_number: fluid.prandtl_number(),
    };

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(output),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_cfd_water_properties_json(_input: *const c_char) -> *mut c_char {

    let fluid = physics_cfd::FluidProperties::water();

    let output = FluidPropertiesOutput {
        kinematic_viscosity: fluid.kinematic_viscosity(),
        thermal_diffusivity: fluid.thermal_diffusivity(),
        prandtl_number: fluid.prandtl_number(),
    };

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(output),
            err: None::<String>,
        })
        .unwrap(),
    )
}

// ============================================================================
// Dimensionless Numbers
// ============================================================================

#[no_mangle]

pub unsafe extern "C" fn rssn_num_cfd_reynolds_number_json(input: *const c_char) -> *mut c_char {

    let input: ReynoldsInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        }
    };

    let re = physics_cfd::reynolds_number(
        input.velocity,
        input.length,
        input.kinematic_viscosity,
    );

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(re),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_cfd_cfl_number_json(input: *const c_char) -> *mut c_char {

    let input: CflInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        }
    };

    let cfl = physics_cfd::cfl_number(
        input.velocity,
        input.dt,
        input.dx,
    );

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(cfl),
            err: None::<String>,
        })
        .unwrap(),
    )
}

// ============================================================================
// 1D Solvers
// ============================================================================

#[no_mangle]

pub unsafe extern "C" fn rssn_num_cfd_solve_advection_1d_json(input: *const c_char) -> *mut c_char {

    let input: Advection1DInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<Vec<f64>>, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        }
    };

    let results = physics_cfd::solve_advection_1d(
        &input.u0,
        input.c,
        input.dx,
        input.dt,
        input.num_steps,
    );

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(results),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_cfd_solve_diffusion_1d_json(input: *const c_char) -> *mut c_char {

    let input: Diffusion1DInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<Vec<f64>>, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        }
    };

    let results = physics_cfd::solve_diffusion_1d(
        &input.u0,
        input.alpha,
        input.dx,
        input.dt,
        input.num_steps,
    );

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(results),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_cfd_solve_advection_diffusion_1d_json(
    input: *const c_char
) -> *mut c_char {

    let input: AdvectionDiffusion1DInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<Vec<f64>>, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        }
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
        serde_json::to_string(&FfiResult {
            ok: Some(results),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_cfd_solve_burgers_1d_json(input: *const c_char) -> *mut c_char {

    let input: Burgers1DInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<Vec<f64>>, String> {
                        ok: None,
                        err: Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        }
    };

    let results = physics_cfd::solve_burgers_1d(
        &input.u0,
        input.nu,
        input.dx,
        input.dt,
        input.num_steps,
    );

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(results),
            err: None::<String>,
        })
        .unwrap(),
    )
}
