//! Bincode-based FFI API for numerical ODE solvers.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::ode::OdeSolverMethod;
use crate::numerical::ode::{
    self,
};
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct OdeInput {
    funcs : Vec<Expr>,
    y0 : Vec<f64>,
    x_range : (f64, f64),
    num_steps : usize,
    method : OdeSolverMethod,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_ode_solve_bincode(buffer : BincodeBuffer) -> BincodeBuffer {

    let input : OdeInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<Vec<f64>>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    let res = ode::solve_ode_system(
        &input.funcs,
        &input.y0,
        input.x_range,
        input.num_steps,
        input.method,
    );

    let ffi_res = match res {
        | Ok(v) => {
            FfiResult {
                ok : Some(v),
                err : None,
            }
        },
        | Err(e) => {
            FfiResult {
                ok : None,
                err : Some(e),
            }
        },
    };

    to_bincode_buffer(&ffi_res)
}
