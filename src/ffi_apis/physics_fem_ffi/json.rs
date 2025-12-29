//! JSON-based FFI API for physics FEM functions.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_fem;

#[derive(Deserialize)]

struct Poisson1DInput {
    n_elements: usize,
    domain_length: f64,
}

/// Solves the 1D Poisson equation using Finite Element Method (FEM) via JSON serialization.
///
/// The Poisson equation -d²u/dx² = f(x) is solved using linear finite elements
/// with homogeneous Dirichlet boundary conditions (u = 0 at boundaries).
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `n_elements`: Number of finite elements in the mesh
///   - `domain_length`: Total length of the 1D domain
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the solution vector u at nodal points.
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

pub unsafe extern "C" fn rssn_physics_fem_solve_poisson_1d_json(
    input: *const c_char
) -> *mut c_char {

    let input : Poisson1DInput = match from_json_string(input) {
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

    match physics_fem::solve_poisson_1d(
        input.n_elements,
        input.domain_length,
        |_| 2.0,
    ) {
        | Ok(res) => {
            to_c_string(
                serde_json::to_string(
                    &FfiResult::<
                        Vec<f64>,
                        String,
                    >::ok(
                        res
                    ),
                )
                .unwrap(),
            )
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(
                    &FfiResult::<
                        Vec<f64>,
                        String,
                    >::err(
                        e
                    ),
                )
                .unwrap(),
            )
        },
    }
}
