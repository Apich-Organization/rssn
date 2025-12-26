//! JSON-based FFI API for physics SM functions.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sm::AdvectionDiffusionConfig;
use crate::physics::physics_sm::{
    self,
};

#[derive(Deserialize)]

struct AdvectionDiffusion1DInput {
    initial_condition : Vec<f64>,
    dx : f64,
    c : f64,
    d : f64,
    dt : f64,
    steps : usize,
}

#[derive(Deserialize)]

struct AdvectionDiffusion2DInput {
    initial_condition : Vec<f64>,
    config : AdvectionDiffusionConfig,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sm_solve_advection_1d_json(
    input : *const c_char
) -> *mut c_char {

    let input: AdvectionDiffusion1DInput =
        match from_json_string(input) {
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

    let res = physics_sm::solve_advection_diffusion_1d(
        &input.initial_condition,
        input.dx,
        input.c,
        input.d,
        input.dt,
        input.steps,
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

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sm_solve_advection_2d_json(
    input : *const c_char
) -> *mut c_char {

    let input: AdvectionDiffusion2DInput =
        match from_json_string(input) {
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

    let res = physics_sm::solve_advection_diffusion_2d(
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
