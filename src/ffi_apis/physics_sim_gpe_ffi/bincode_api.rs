//! Bincode-based FFI API for physics sim GPE superfluidity functions.

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::gpe_superfluidity::GpeParameters;
use crate::physics::physics_sim::gpe_superfluidity::{
    self,
};

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_gpe_run_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let params : GpeParameters = match from_bincode_buffer(&buffer) {
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

    match gpe_superfluidity::run_gpe_ground_state_finder(&params) {
        | Ok(res) => {
            to_bincode_buffer(&FfiResult::<
                Vec<f64>,
                String,
            >::ok(
                res.into_raw_vec(),
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
