//! JSON-based FFI API for physics sim Ising statistical functions.

use crate::ffi_apis::common::{
    from_json_string,
    to_c_string,
};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::ising_statistical::{
    self,
    IsingParameters,
};
use serde::{
    Deserialize,
    Serialize,
};
use std::os::raw::c_char;

#[derive(Serialize)]

struct IsingOutput {
    pub grid: Vec<i8>,
    pub magnetization: f64,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_ising_run_json(input: *const c_char) -> *mut c_char {

    let params: IsingParameters = match from_json_string(input) {
        Some(p) => p,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    IsingOutput,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        }
    };

    let (grid, magnetization) = ising_statistical::run_ising_simulation(&params);

    let out = IsingOutput {
        grid,
        magnetization,
    };

    to_c_string(
        serde_json::to_string(&FfiResult::<
            IsingOutput,
            String,
        >::ok(out))
        .unwrap(),
    )
}
