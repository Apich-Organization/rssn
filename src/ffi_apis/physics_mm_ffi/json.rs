//! JSON-based FFI API for physics MM functions.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_mm::Particle;
use crate::physics::physics_mm::SPHSystem;
use crate::physics::physics_mm::{
    self,
};

#[derive(Deserialize)]

struct SphInput {
    system : SPHSystem,
    dt : f64,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_mm_sph_update_json(
    input : *const c_char
) -> *mut c_char {

    let mut input : SphInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    SPHSystem,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    input
        .system
        .update(input.dt);

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                SPHSystem,
                String,
            >::ok(
                input.system
            ),
        )
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_mm_simulate_dam_break_json(
) -> *mut c_char {

    let res = physics_mm::simulate_dam_break_2d_scenario();

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                Vec<(f64, f64)>,
                String,
            >::ok(res),
        )
        .unwrap(),
    )
}
