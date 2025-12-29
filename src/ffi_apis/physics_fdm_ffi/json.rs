//! JSON-based FFI API for physics FDM functions.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_fdm::Dimensions;
use crate::physics::physics_fdm::FdmGrid;
use crate::physics::physics_fdm::{
    self,
};

#[derive(Deserialize)]

struct HeatEquationInput {
    width: usize,
    height: usize,
    alpha: f64,
    dx: f64,
    dy: f64,
    dt: f64,
    steps: usize,
    initial_temp: f64, /* simplified for JSON: constant initial temp except source */
}

#[derive(Deserialize)]

struct WaveEquationInput {
    width: usize,
    height: usize,
    c: f64,
    dx: f64,
    dy: f64,
    dt: f64,
    steps: usize,
}

#[derive(Deserialize)]

struct PoissonInput {
    width: usize,
    height: usize,
    source: Vec<f64>,
    dx: f64,
    dy: f64,
    omega: f64,
    max_iter: usize,
    tolerance: f64,
}

#[derive(Deserialize)]

struct BurgersInput {
    initial_u: Vec<f64>,
    dx: f64,
    nu: f64,
    dt: f64,
    steps: usize,
}

/// Solves the 2D heat equation using Finite Difference Method (FDM) via JSON serialization.
///
/// The heat equation ∂u/∂t = α∇²u is solved using explicit finite differences with
/// a square heat source at the center of the domain.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `width`: Grid width (number of cells in x-direction)
///   - `height`: Grid height (number of cells in y-direction)
///   - `alpha`: Thermal diffusivity coefficient α
///   - `dx`: Spatial step size in x-direction
///   - `dy`: Spatial step size in y-direction
///   - `dt`: Time step size
///   - `steps`: Number of time steps to simulate
///   - `initial_temp`: Temperature of the central heat source
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<FdmGrid<f64>, String>` with
/// the final temperature field grid.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_fdm_heat_json(
    input: *const c_char
) -> *mut c_char {

    let input : HeatEquationInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    FdmGrid<f64>,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let result = physics_fdm::solve_heat_equation_2d(
        input.width,
        input.height,
        input.alpha,
        input.dx,
        input.dy,
        input.dt,
        input.steps,
        |x, y| {
            if x > input.width / 3
                && x < 2 * input.width / 3
                && y > input.height / 3
                && y < 2 * input.height / 3
            {

                input.initial_temp
            } else {

                0.0
            }
        },
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                FdmGrid<f64>,
                String,
            >::ok(result),
        )
        .unwrap(),
    )
}

/// Solves the 2D wave equation using Finite Difference Method (FDM) via JSON serialization.
///
/// The wave equation ∂²u/∂t² = c²∇²u is solved using explicit finite differences with
/// a Gaussian initial condition centered at the grid midpoint.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `width`: Grid width (number of cells in x-direction)
///   - `height`: Grid height (number of cells in y-direction)
///   - `c`: Wave speed
///   - `dx`: Spatial step size in x-direction
///   - `dy`: Spatial step size in y-direction
///   - `dt`: Time step size
///   - `steps`: Number of time steps to simulate
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<FdmGrid<f64>, String>` with
/// the final wave field grid.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_fdm_wave_json(
    input: *const c_char
) -> *mut c_char {

    let input : WaveEquationInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    FdmGrid<f64>,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let result = physics_fdm::solve_wave_equation_2d(
        input.width,
        input.height,
        input.c,
        input.dx,
        input.dy,
        input.dt,
        input.steps,
        |x, y| {

            let dx_cen = x as f64 - (input.width / 2) as f64;

            let dy_cen = y as f64 - (input.height / 2) as f64;

            let dist2 = dx_cen.powi(2) + dy_cen.powi(2);

            (-dist2 / 20.0).exp()
        },
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                FdmGrid<f64>,
                String,
            >::ok(result),
        )
        .unwrap(),
    )
}

/// Solves the 2D Poisson equation using Finite Difference Method with SOR via JSON serialization.
///
/// The Poisson equation ∇²u = f is solved using Successive Over-Relaxation (SOR)
/// iteration to find the steady-state potential field given a source distribution.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `width`: Grid width (number of cells in x-direction)
///   - `height`: Grid height (number of cells in y-direction)
///   - `source`: Source term f (flattened 2D array)
///   - `dx`: Spatial step size in x-direction
///   - `dy`: Spatial step size in y-direction
///   - `omega`: SOR relaxation parameter (1 < ω < 2 for optimal convergence)
///   - `max_iter`: Maximum number of iterations
///   - `tolerance`: Convergence tolerance for residual norm
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<FdmGrid<f64>, String>` with
/// the solution grid u.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_fdm_poisson_json(
    input: *const c_char
) -> *mut c_char {

    let input : PoissonInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    FdmGrid<f64>,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let source_grid =
        FdmGrid::from_data(
            input.source,
            Dimensions::D2(
                input.width,
                input.height,
            ),
        );

    let result =
        physics_fdm::solve_poisson_2d(
            input.width,
            input.height,
            &source_grid,
            input.dx,
            input.dy,
            input.omega,
            input.max_iter,
            input.tolerance,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                FdmGrid<f64>,
                String,
            >::ok(result),
        )
        .unwrap(),
    )
}

/// Solves the 1D Burgers' equation using Finite Difference Method via JSON serialization.
///
/// Burgers' equation ∂u/∂t + u∂u/∂x = ν∂²u/∂x² combines nonlinear convection
/// with diffusion, modeling shock wave formation and viscous fluid flow.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `initial_u`: Initial velocity field
///   - `dx`: Spatial step size
///   - `nu`: Kinematic viscosity coefficient ν
///   - `dt`: Time step size
///   - `steps`: Number of time steps to simulate
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the final velocity field.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_fdm_burgers_json(
    input: *const c_char
) -> *mut c_char {

    let input : BurgersInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let result =
        physics_fdm::solve_burgers_1d(
            &input.initial_u,
            input.dx,
            input.nu,
            input.dt,
            input.steps,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                Vec<f64>,
                String,
            >::ok(result),
        )
        .unwrap(),
    )
}
