//! Bincode-based FFI API for physics MTM functions.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_mtm;

#[derive(Deserialize)]

struct Multigrid2DInput {
    n: usize,
    f: Vec<f64>,
    num_cycles: usize,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_mtm_solve_poisson_2d_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : Multigrid2DInput = match from_bincode_buffer(&buffer) {
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

    match physics_mtm::solve_poisson_2d_multigrid(
        input.n,
        &input.f,
        input.num_cycles,
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
