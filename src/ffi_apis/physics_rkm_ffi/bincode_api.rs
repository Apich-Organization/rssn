//! Bincode-based FFI API for physics RKM (Runge-Kutta Methods) functions.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_rkm::DormandPrince54;
use crate::physics::physics_rkm::{
    self,
};

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

/// Solves the Lorenz system using adaptive Dormand-Prince RK5(4) method via bincode serialization.
///
/// The Lorenz system is a chaotic dynamical system defined by:
/// dx/dt = σ(y - x), dy/dt = x(ρ - z) - y, dz/dt = xy - βz.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `LorenzInput` with:
///   - `sigma`: Prandtl number σ
///   - `rho`: Rayleigh number ρ
///   - `beta`: Geometric parameter β
///   - `y0`: Initial state [x₀, y₀, z₀]
///   - `t_span`: Time interval (`t_start`, `t_end`)
///   - `dt_initial`: Initial time step size
///   - `tol`: Error tolerances (absolute, relative)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<OdeResult, String>` with either:
/// - `ok`: Object containing `time` and `states` arrays
/// - `err`: Error message if computation failed
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_rkm_lorenz_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : LorenzInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                OdeResult,
                String,
            >::err(
                "Invalid Bincode".to_string(),
            ))
        },
    };

    let system =
        physics_rkm::LorenzSystem {
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

    let mut time = Vec::with_capacity(
        results.len(),
    );

    let mut states = Vec::with_capacity(
        results.len(),
    );

    for (t, y) in results {

        time.push(t);

        states.push(y);
    }

    to_bincode_buffer(&FfiResult::<
        OdeResult,
        String,
    >::ok(
        OdeResult {
            time,
            states,
        },
    ))
}
