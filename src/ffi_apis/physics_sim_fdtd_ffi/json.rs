//! JSON-based FFI API for physics sim FDTD electrodynamics functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::fdtd_electrodynamics::FdtdParameters;
use crate::physics::physics_sim::fdtd_electrodynamics::{
    self,
};

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_fdtd_run_json(
    input: *const c_char
) -> *mut c_char {

    let params : FdtdParameters = match from_json_string(input) {
        | Some(p) => p,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<Vec<f64>>,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let snapshots = fdtd_electrodynamics::run_fdtd_simulation(&params);

    // Return only the last Ez field as Vec<Vec<f64>> for simplicity in JSON
    if let Some(final_ez) =
        snapshots.last()
    {

        let mut out =
            Vec::with_capacity(
                params.width,
            );

        for row in final_ez
            .axis_iter(ndarray::Axis(0))
        {

            out.push(row.to_vec());
        }

        to_c_string(
            serde_json::to_string(
                &FfiResult::<
                    Vec<Vec<f64>>,
                    String,
                >::ok(
                    out
                ),
            )
            .unwrap(),
        )
    } else {

        to_c_string(
            serde_json::to_string(
                &FfiResult::<
                    Vec<Vec<f64>>,
                    String,
                >::err(
                    "No snapshots"
                        .to_string(),
                ),
            )
            .unwrap(),
        )
    }
}
