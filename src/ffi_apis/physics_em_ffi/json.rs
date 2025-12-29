//! JSON-based FFI API for physics EM functions.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_em::{
    self,
};
use crate::physics::physics_rkm::DampedOscillatorSystem;
use crate::physics::physics_rkm::LorenzSystem;

#[derive(Deserialize)]

struct EulerInput {
    system_type: String, /* "lorenz", "oscillator", "orbital" */
    params: serde_json::Value,
    y0: Vec<f64>,
    t_span: (f64, f64),
    dt: f64,
    method: String, /* "forward", "midpoint", "heun" */
}

/// Solves ODE systems using Euler methods (forward, midpoint, or Heun) via JSON serialization.
///
/// Supports various dynamical systems including Lorenz attractor, damped oscillators,
/// and orbital mechanics, using explicit Euler integration methods.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `system_type`: System identifier ("lorenz", "oscillator", "orbital")
///   - `params`: System parameters as JSON object
///   - `y0`: Initial state vector
///   - `t_span`: Time interval [`t_start`, `t_end`]
///   - `dt`: Time step size
///   - `method`: Integration method ("forward", "midpoint", "heun")
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<(f64, Vec<f64>)>, String>` with
/// the trajectory as (time, state) pairs.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_em_solve_json(
    input: *const c_char
) -> *mut c_char {

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
        | _ => return to_c_string(
            serde_json::to_string(
                &FfiResult::<
                    Vec<(
                        f64,
                        Vec<f64>,
                    )>,
                    String,
                >::err(
                    "Unknown system \
                     type"
                        .to_string(),
                ),
            )
            .unwrap(),
        ),
    };

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                Vec<(f64, Vec<f64>)>,
                String,
            >::ok(res),
        )
        .unwrap(),
    )
}

fn solve_with_method<S : crate::physics::physics_rkm::OdeSystem>(
    sys : &S,
    y0 : &[f64],
    t_span : (f64, f64),
    dt : f64,
    method : &str,
) -> Vec<(f64, Vec<f64>)>{

    match method {
        | "midpoint" => physics_em::solve_midpoint_euler(sys, y0, t_span, dt),
        | "heun" => physics_em::solve_heun_euler(sys, y0, t_span, dt),
        | _ => physics_em::solve_forward_euler(sys, y0, t_span, dt),
    }
}
