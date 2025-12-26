//! JSON-based FFI API for physics CNM functions.

use crate::ffi_apis::common::{
    from_json_string,
    to_c_string,
};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_cnm::{
    self,
    HeatEquationSolverConfig,
};
use serde::{
    Deserialize,
    Serialize,
};
use std::os::raw::c_char;

#[derive(Deserialize)]

struct Heat2DInput {
    initial_condition: Vec<f64>,
    config: HeatEquationSolverConfig,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_cnm_solve_heat_2d_json(
    input: *const c_char
) -> *mut c_char {

    let input: Heat2DInput = match from_json_string(input) {
        | Some(i) => i,
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

    let res = physics_cnm::solve_heat_equation_2d_cn_adi(
        &input.initial_condition,
        &input.config,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                Vec<f64>,
                String,
            >::ok(res),
        )
        .unwrap(),
    )
}
