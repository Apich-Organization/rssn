//! JSON-based FFI API for physics sim geodesic relativity functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::geodesic_relativity::GeodesicParameters;
use crate::physics::physics_sim::geodesic_relativity::{
    self,
};

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_geodesic_run_json(
    input: *const c_char
) -> *mut c_char {

    let params : GeodesicParameters = match from_json_string(input) {
        | Some(p) => p,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<(f64, f64)>,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let path = geodesic_relativity::run_geodesic_simulation(&params);

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                Vec<(f64, f64)>,
                String,
            >::ok(path),
        )
        .unwrap(),
    )
}
