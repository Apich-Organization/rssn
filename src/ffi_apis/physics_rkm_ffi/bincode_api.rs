//! Bincode-based FFI API for physics RKM (Runge-Kutta Methods) functions.

use crate::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer, BincodeBuffer};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_rkm::{self, DormandPrince54};
use serde::{Deserialize, Serialize};

#[derive(Deserialize)]

struct LorenzInput {
    sigma: f64,
    rho: f64,
    beta: f64,
    y0: Vec<f64>,
    t_span: (f64, f64),
    dt_initial: f64,
    tol: (f64, f64),
}

#[derive(Serialize)]

struct OdeResult {
    time: Vec<f64>,
    states: Vec<Vec<f64>>,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_physics_rkm_lorenz_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: LorenzInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<
                OdeResult,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
        }
    };

    let system = physics_rkm::LorenzSystem {
        sigma: input.sigma,
        rho: input.rho,
        beta: input.beta,
    };

    let solver = DormandPrince54::new();

    let results = solver.solve(
        &system,
        &input.y0,
        input.t_span,
        input.dt_initial,
        input.tol,
    );

    let mut time = Vec::with_capacity(results.len());

    let mut states = Vec::with_capacity(results.len());

    for (t, y) in results {

        time.push(t);

        states.push(y);
    }

    to_bincode_buffer(&FfiResult::<
        OdeResult,
        String,
    >::ok(
        OdeResult { time, states },
    ))
}
