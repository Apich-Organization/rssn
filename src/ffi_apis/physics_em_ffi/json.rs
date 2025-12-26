//! JSON-based FFI API for physics EM functions.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_em::OrbitalSystem;
use crate::physics::physics_em::{
    self,
};
use crate::physics::physics_rkm::DampedOscillatorSystem;
use crate::physics::physics_rkm::LorenzSystem;

#[derive(Deserialize)]

struct EulerInput {
    system_type : String, // "lorenz", "oscillator", "orbital"
    params : serde_json::Value,
    y0 : Vec<f64>,
    t_span : (f64, f64),
    dt : f64,
    method : String, // "forward", "midpoint", "heun"
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_em_solve_json(input : *const c_char) -> *mut c_char {

    let input : EulerInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<(f64, Vec<f64>)>,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let res = match input
        .system_type
        .as_str()
    {
        | "lorenz" => {

            let sys : LorenzSystem = match serde_json::from_value(input.params) {
                | Ok(s) => s,
                | Err(e) => {
                    return to_c_string(
                        serde_json::to_string(&FfiResult::<
                            Vec<(f64, Vec<f64>)>,
                            String,
                        >::err(
                            e.to_string()
                        ))
                        .unwrap(),
                    )
                },
            };

            solve_with_method(
                &sys,
                &input.y0,
                input.t_span,
                input.dt,
                &input.method,
            )
        },
        | "oscillator" => {

            let sys : DampedOscillatorSystem = match serde_json::from_value(input.params) {
                | Ok(s) => s,
                | Err(e) => {
                    return to_c_string(
                        serde_json::to_string(&FfiResult::<
                            Vec<(f64, Vec<f64>)>,
                            String,
                        >::err(
                            e.to_string()
                        ))
                        .unwrap(),
                    )
                },
            };

            solve_with_method(
                &sys,
                &input.y0,
                input.t_span,
                input.dt,
                &input.method,
            )
        },
        | _ => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<(f64, Vec<f64>)>,
                    String,
                >::err(
                    "Unknown system type".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    to_c_string(
        serde_json::to_string(&FfiResult::<
            Vec<(f64, Vec<f64>)>,
            String,
        >::ok(res))
        .unwrap(),
    )
}

fn solve_with_method<S : crate::physics::physics_rkm::OdeSystem>(
    sys : &S,
    y0 : &[f64],
    t_span : (f64, f64),
    dt : f64,
    method : &str,
) -> Vec<(f64, Vec<f64>)> {

    match method {
        | "midpoint" => physics_em::solve_midpoint_euler(sys, y0, t_span, dt),
        | "heun" => physics_em::solve_heun_euler(sys, y0, t_span, dt),
        | _ => physics_em::solve_forward_euler(sys, y0, t_span, dt),
    }
}
