//! JSON-based FFI API for physics MTM functions.

use crate::ffi_apis::common::{from_json_string, to_c_string};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_mtm;
use serde::{Deserialize, Serialize};
use std::os::raw::c_char;

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

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_mtm_solve_poisson_1d_json(
    input: *const c_char
) -> *mut c_char {

    let input: Multigrid1DInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        }
    };

    match physics_mtm::solve_poisson_1d_multigrid(
        input.n_interior,
        &input.f,
        input.num_cycles,
    ) {
        Ok(res) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::ok(
                    res
                ))
                .unwrap(),
            )
        }
        Err(e) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    e
                ))
                .unwrap(),
            )
        }
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_mtm_solve_poisson_2d_json(
    input: *const c_char
) -> *mut c_char {

    let input: Multigrid2DInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        }
    };

    match physics_mtm::solve_poisson_2d_multigrid(
        input.n,
        &input.f,
        input.num_cycles,
    ) {
        Ok(res) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::ok(
                    res
                ))
                .unwrap(),
            )
        }
        Err(e) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    e
                ))
                .unwrap(),
            )
        }
    }
}
