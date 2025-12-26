//! JSON-based FFI API for physics sim GPE superfluidity functions.

use crate::ffi_apis::common::{
    from_json_string,
    to_c_string,
};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::gpe_superfluidity::{
    self,
    GpeParameters,
};
use std::os::raw::c_char;

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_gpe_run_json(input: *const c_char) -> *mut c_char {

    let params: GpeParameters = match from_json_string(input) {
        | Some(p) => p,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    match gpe_superfluidity::run_gpe_ground_state_finder(&params) {
        | Ok(res) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::ok(
                    res.into_raw_vec(),
                ))
                .unwrap(),
            )
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<f64>,
                    String,
                >::err(
                    e
                ))
                .unwrap(),
            )
        },
    }
}
