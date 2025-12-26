//! Bincode-based FFI API for physics CNM functions.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_cnm::HeatEquationSolverConfig;
use crate::physics::physics_cnm::{
    self,
};

#[derive(Deserialize)]

struct Heat2DInput {
    initial_condition : Vec<f64>,
    config : HeatEquationSolverConfig,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_cnm_solve_heat_2d_bincode(
    buffer : BincodeBuffer
) -> BincodeBuffer {

    let input : Heat2DInput = match from_bincode_buffer(&buffer) {
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

    let res = physics_cnm::solve_heat_equation_2d_cn_adi(
        &input.initial_condition,
        &input.config,
    );

    to_bincode_buffer(&FfiResult::<
        Vec<f64>,
        String,
    >::ok(res))
}
