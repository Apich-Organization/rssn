//! JSON-based FFI API for physics MTM functions.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_mtm;

#[derive(Deserialize)]

struct Multigrid1DInput {
    n_interior: usize,
    f: Vec<f64>,
    num_cycles: usize,
}

#[derive(Deserialize)]

struct Multigrid2DInput {
    n: usize,
    f: Vec<f64>,
    num_cycles: usize,
}

/// Solves the 1D Poisson equation using Multigrid Method via JSON serialization.
///
/// The Poisson equation -d²u/dx² = f is solved using the multigrid method, which achieves
/// optimal O(N) complexity through hierarchical coarse-grid correction.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `n_interior`: Number of interior grid points
///   - `f`: Right-hand side source term
///   - `num_cycles`: Number of V-cycles or W-cycles to perform
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the solution vector u.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_mtm_solve_poisson_1d_json(
    input: *const c_char
) -> *mut c_char {

    let input : Multigrid1DInput = match from_json_string(input) {
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

    match physics_mtm::solve_poisson_1d_multigrid(
        input.n_interior,
        &input.f,
        input.num_cycles,
    ) {
        | Ok(res) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::ok(
                    res
                ))
                .unwrap(),
            )
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    e
                ))
                .unwrap(),
            )
        },
    }
}

/// Solves the 2D Poisson equation using Multigrid Method via JSON serialization.
///
/// The Poisson equation ∇²u = f is solved using the multigrid method, which achieves
/// optimal O(N) complexity through hierarchical coarse-grid correction.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `n`: Grid size (n×n interior points)
///   - `f`: Right-hand side source term (flattened 2D array)
///   - `num_cycles`: Number of V-cycles or W-cycles to perform
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the solution vector u (flattened 2D array).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_mtm_solve_poisson_2d_json(
    input: *const c_char
) -> *mut c_char {

    let input : Multigrid2DInput = match from_json_string(input) {
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

    match physics_mtm::solve_poisson_2d_multigrid(
        input.n,
        &input.f,
        input.num_cycles,
    ) {
        | Ok(res) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::ok(
                    res
                ))
                .unwrap(),
            )
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    e
                ))
                .unwrap(),
            )
        },
    }
}
