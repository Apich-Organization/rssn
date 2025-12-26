//! Bincode-based FFI API for physics sim linear elasticity functions.

use crate::ffi_apis::common::{
    from_bincode_buffer,
    to_bincode_buffer,
    BincodeBuffer,
};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::linear_elasticity::{
    self,
    ElasticityParameters,
};

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_linear_elasticity_run_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let params: ElasticityParameters =
        match from_bincode_buffer(&buffer) {
            | Some(p) => p,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    "Invalid Bincode".to_string(),
                ))
            },
        };

    match linear_elasticity::run_elasticity_simulation(
        &params,
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
