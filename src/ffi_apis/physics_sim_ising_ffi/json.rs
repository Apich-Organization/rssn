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
///   - `temperature`: Temperature T in units of `J/k_B`
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
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

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
