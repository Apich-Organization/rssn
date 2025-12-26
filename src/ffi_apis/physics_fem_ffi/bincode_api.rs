//! Bincode-based FFI API for physics FEM functions.

use crate::ffi_apis::common::{
    from_bincode_buffer,
    to_bincode_buffer,
    BincodeBuffer,
};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_fem;
use serde::Deserialize;

#[derive(Deserialize)]

struct Poisson1DInput {
    n_elements: usize,
    domain_length: f64,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_fem_solve_poisson_1d_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: Poisson1DInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                Vec<f64>,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
        },
    };

    match physics_fem::solve_poisson_1d(
        input.n_elements,
        input.domain_length,
        |_| 2.0,
    ) {
        | Ok(res) => {
            to_bincode_buffer(&FfiResult::<
                Vec<f64>,
                String,
            >::ok(
                res
            ))
        },
        | Err(e) => {
            to_bincode_buffer(&FfiResult::<
                Vec<f64>,
                String,
            >::err(
                e
            ))
        },
    }
}
