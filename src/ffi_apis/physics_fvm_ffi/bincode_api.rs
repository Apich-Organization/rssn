//! Bincode-based FFI API for physics FVM functions.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_fvm::SweState;
use crate::physics::physics_fvm::{
    self,
};

#[derive(Deserialize)]

struct SweInput {
    h: Vec<f64>,
    hu: Vec<f64>,
    dx: f64,
    dt: f64,
    steps: usize,
    g: f64,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_fvm_swe_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : SweInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                Vec<SweState>,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
        },
    };

    let result = physics_fvm::solve_shallow_water_1d(
        input.h,
        input.hu,
        input.dx,
        input.dt,
        input.steps,
        input.g,
    );

    to_bincode_buffer(&FfiResult::<
        Vec<SweState>,
        String,
    >::ok(result))
}
