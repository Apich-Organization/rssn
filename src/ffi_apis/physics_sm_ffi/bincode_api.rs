//! Bincode-based FFI API for physics SM functions.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sm::AdvectionDiffusionConfig as AdvectionDiffusionConfig2d;
use crate::physics::physics_sm::{
    self,
}; // Alias for clarity if needed

#[derive(Deserialize)]

struct AdvectionDiffusion2DInput {
    initial_condition : Vec<f64>,
    config : physics_sm::AdvectionDiffusionConfig,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sm_solve_advection_2d_bincode(
    buffer : BincodeBuffer
) -> BincodeBuffer {

    let input : AdvectionDiffusion2DInput = match from_bincode_buffer(&buffer) {
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

    let res = physics_sm::solve_advection_diffusion_2d(
        &input.initial_condition,
        &input.config,
    );

    to_bincode_buffer(&FfiResult::<
        Vec<f64>,
        String,
    >::ok(res))
}
