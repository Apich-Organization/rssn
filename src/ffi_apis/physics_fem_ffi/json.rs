//! JSON-based FFI API for physics FEM functions.

use crate::ffi_apis::common::{from_json_string, to_c_string};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_fem;
use serde::Deserialize;
use std::os::raw::c_char;

#[derive(Deserialize)]
struct Poisson1DInput {
    n_elements: usize,
    domain_length: f64,
}

#[no_mangle]
pub unsafe extern "C" fn rssn_physics_fem_solve_poisson_1d_json(
    input: *const c_char,
) -> *mut c_char {
    let input: Poisson1DInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vec<f64>, String>::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        }
    };

    match physics_fem::solve_poisson_1d(input.n_elements, input.domain_length, |_| 2.0) {
        Ok(res) => {
            to_c_string(serde_json::to_string(&FfiResult::<Vec<f64>, String>::ok(res)).unwrap())
        }
        Err(e) => {
            to_c_string(serde_json::to_string(&FfiResult::<Vec<f64>, String>::err(e)).unwrap())
        }
    }
}
