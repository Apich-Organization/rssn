//! JSON-based FFI API for physics sim Ising statistical functions.

use std::os::raw::c_char;

use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::physics::physics_sim::ising_statistical::IsingParameters;
use crate::physics::physics_sim::ising_statistical::{
    self,
};

#[derive(Serialize)]

struct IsingOutput {
    pub grid: Vec<i8>,
    pub magnetization: f64,
}

/// Runs a 2D Ising model Monte Carlo simulation using the Metropolis algorithm via JSON serialization.
///
/// The Ising model with Hamiltonian H = -J∑⟨i,j⟩sᵢsⱼ - h∑ᵢsᵢ describes phase transitions
/// in magnetic systems. The simulation uses Metropolis-Hastings sampling to evolve spin
/// configurations toward thermal equilibrium.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `width`, `height`: Grid dimensions
///   - `temperature`: Temperature T in units of J/k_B
///   - `mc_steps`: Number of Monte Carlo sweeps to perform
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<IsingOutput, String>` with:
/// - `grid`: Final spin configuration (±1 values)
/// - `magnetization`: Average magnetization M = ⟨∑ᵢsᵢ⟩/N
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_physics_sim_ising_run_json(
    input: *const c_char
) -> *mut c_char {

    let params : IsingParameters = match from_json_string(input) {
        | Some(p) => p,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    IsingOutput,
                    String,
                >::err(
                    "Invalid JSON".to_string(),
                ))
                .unwrap(),
            )
        },
    };

    let (grid, magnetization) = ising_statistical::run_ising_simulation(&params);

    let out = IsingOutput {
        grid,
        magnetization,
    };

    to_c_string(
        serde_json::to_string(
            &FfiResult::<
                IsingOutput,
                String,
            >::ok(out),
        )
        .unwrap(),
    )
}
