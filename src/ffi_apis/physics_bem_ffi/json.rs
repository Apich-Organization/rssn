//! JSON-based FFI API for physics BEM functions.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_bem::BoundaryCondition;
use crate::physics::physics_bem::{
    self,
};

#[derive(Deserialize)]

struct Bem2DInput {
    points: Vec<(f64, f64)>,
    bcs: Vec<BemBoundaryCondition>,
}

#[derive(Deserialize, Serialize)]

enum BemBoundaryCondition {
    Potential(f64),
    Flux(f64),
}

#[derive(Serialize)]

struct Bem2DOutput {
    u: Vec<f64>,
    q: Vec<f64>,
}

/// Solves the 2D Laplace equation using Boundary Element Method (BEM) via JSON serialization.
///
/// The Laplace equation ∇²u = 0 is solved using BEM, where the domain is discretized
/// into boundary elements and the solution is represented by potential u and flux q
/// on the boundary.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `points`: Boundary points as (x, y) coordinate pairs
///   - `bcs`: Boundary conditions (Potential or Flux) at each point
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Bem2DOutput, String>` with
/// potential `u` and flux `q` values at boundary nodes.
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

pub unsafe extern "C" fn rssn_physics_bem_solve_laplace_2d_json(
    input: *const c_char
) -> *mut c_char {

    let input : Bem2DInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Bem2DOutput,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let bcs : Vec<BoundaryCondition<f64>> = input
        .bcs
        .into_iter()
        .map(|bc| {

            match bc {
                | BemBoundaryCondition::Potential(v) => BoundaryCondition::Potential(v),
                | BemBoundaryCondition::Flux(v) => BoundaryCondition::Flux(v),
            }
        })
        .collect();

    match physics_bem::solve_laplace_bem_2d(&input.points, &bcs) {
        | Ok((u, q)) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Bem2DOutput,
                    String,
                >::ok(
                    Bem2DOutput {
                        u,
                        q,
                    },
                ))
                .unwrap(),
            )
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Bem2DOutput,
                    String,
                >::err(
                    e
                ))
                .unwrap(),
            )
        },
    }
}
