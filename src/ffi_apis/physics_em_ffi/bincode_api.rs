//! Bincode-based FFI API for physics EM functions.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_em;
use crate::physics::physics_rkm::DampedOscillatorSystem;
use crate::physics::physics_rkm::LorenzSystem;

#[derive(Deserialize)]

struct EulerInput {
    system_type: String,
    params_bincode: Vec<u8>,
    y0: Vec<f64>,
    t_span: (f64, f64),
    dt: f64,
    method: String,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_em_solve_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : EulerInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                Vec<(f64, Vec<f64>)>,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
        },
    };

    let res = match input
        .system_type
        .as_str()
    {
        | "lorenz" => {

            let (sys, _) : (LorenzSystem, usize) = match bincode_next::serde::decode_from_slice(
                &input.params_bincode,
                bincode_next::config::standard(),
            ) {
                | Ok(s) => s,
                | Err(e) => {
                    return to_bincode_buffer(&FfiResult::<
                        Vec<(f64, Vec<f64>)>,
                        String,
                    >::err(
                        e.to_string()
                    ))
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

            let (sys, _) : (
                DampedOscillatorSystem,
                usize,
            ) = match bincode_next::serde::decode_from_slice(
                &input.params_bincode,
                bincode_next::config::standard(),
            ) {
                | Ok(s) => s,
                | Err(e) => {
                    return to_bincode_buffer(&FfiResult::<
                        Vec<(f64, Vec<f64>)>,
                        String,
                    >::err(
                        e.to_string()
                    ))
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
            return to_bincode_buffer(
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
        },
    };

    to_bincode_buffer(&FfiResult::<
        Vec<(f64, Vec<f64>)>,
        String,
    >::ok(res))
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
