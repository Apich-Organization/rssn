//! JSON-based FFI API for physics BEM functions.

use crate::ffi_apis::common::{
    from_json_string,
    to_c_string,
};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_bem::{
    self,
    BoundaryCondition,
};
use serde::{
    Deserialize,
    Serialize,
};
use std::os::raw::c_char;

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

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_bem_solve_laplace_2d_json(
    input: *const c_char
) -> *mut c_char {

    let input: Bem2DInput = match from_json_string(input) {
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

    let bcs: Vec<BoundaryCondition<f64>> = input
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
                    Bem2DOutput { u, q },
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
